import re
from pathlib import Path
from zipfile import BadZipFile, ZipFile

import pandas as pd
import requests
from tqdm import tqdm

from protein_graph.data_sources.base import BaseDataSource

#: Matches the concrete release version embedded in the inner archive
#: member name, e.g. "BIOGRID-ALL-5.0.260.tab3.txt" -> "5.0.260". BioGRID's
#: "LATEST" download alias always resolves to one specific numbered release;
#: this regex recovers that concrete version so builds stay reproducible
#: even though the download URL itself never names a version.
_RELEASE_VERSION_PATTERN = re.compile(
    r"BIOGRID-ALL-(?P<version>[0-9][0-9A-Za-z.\-]*)\.tab3\.txt$"
)

#: Every real BioGRID release version looks like "5.0.260": digits, dots,
#: and hyphens only. A user-supplied ``release`` is validated against this
#: same shape before it is interpolated into a URL and a local cache file
#: path - without this, a value such as "../../../etc/passwd" would be
#: interpolated into ``archive_path`` (raw_dir / f"BIOGRID-ALL-{release}...")
#: and, because it contains path separators and "..", would resolve OUTSIDE
#: raw_dir, turning a downloaded archive into an arbitrary local file write.
_SAFE_RELEASE_PATTERN = re.compile(r"^[0-9][0-9A-Za-z.\-]*$")


class BioGRIDDataSource(BaseDataSource):
    """Global BioGRID physical interaction evidence data source.

    Ingests the official global ``BIOGRID-ALL-<version>.tab3.zip`` product,
    which covers every organism BioGRID curates - including cross-species
    and host-pathogen interactions - rather than ``BIOGRID-ORGANISM`` (one
    species) or ``BIOGRID-MV-Physical`` (multi-validated only, which can
    discard a single strong direct-binding experiment because it measures
    corroboration rather than evidence quality).

    The normalized output contains individual BioGRID physical evidence
    records. It does not resolve identifiers into canonical ProteinGraph
    protein IDs, does not aggregate multiple evidence rows into unique
    protein pairs, and does not decide which rows are strong/direct enough
    for the training set - see
    :mod:`protein_graph.data_sources.biogrid_evidence_policy` and
    :class:`protein_graph.pipelines.biogrid_pair_acceptance.BioGRIDPairAcceptancePipeline`
    for those decisions. This class's only job is faithful ingestion.
    """

    #: Default when ``release`` is not pinned: BioGRID's "Latest-Release"
    #: alias, which always resolves to whatever the current numbered
    #: release is. This is an explicitly documented default, not a silent
    #: one - the concrete version actually served is still recovered from
    #: the downloaded archive's inner file name and recorded as
    #: ``release_version`` (and, downstream, in every build report), so a
    #: build is never silently unreproducible just because it used the
    #: default.
    LATEST_RELEASE_URL = (
        "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/"
        "BIOGRID-ALL-LATEST.tab3.zip"
    )

    #: Backwards-compatible alias.
    DOWNLOAD_URL = LATEST_RELEASE_URL

    #: Template for a specific pinned release, e.g. release="5.0.255" ->
    #: .../Release-Archive/BIOGRID-5.0.255/BIOGRID-ALL-5.0.255.tab3.zip
    #: (verified reachable against BioGRID's real archive on 2026-08-24).
    _RELEASE_ARCHIVE_URL_TEMPLATE = (
        "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/"
        "BIOGRID-{release}/BIOGRID-ALL-{release}.tab3.zip"
    )

    def __init__(
        self,
        *,
        interaction_type: str = "physical",
        taxonomy_id: int | None = None,
        species: str | None = None,
        raw_dir: str | Path = "data/raw/biogrid",
        sample_limit: int | None = None,
        release: str | None = None,
    ) -> None:
        """Configure the global BioGRID data source.

        Parameters
        ----------
        interaction_type:
            BioGRID "Experimental System Type" to keep. Defaults to
            "physical", which excludes purely genetic-interaction rows.
            This is a structural filter (physical vs. genetic evidence
            category), not an evidence-strength judgment.
        taxonomy_id:
            Optional development/debug narrowing filter. When set, only
            rows where at least one interactor belongs to this NCBI
            taxonomy ID are kept (the partner may belong to a different
            organism, so cross-species pairs involving this organism are
            still retained). Production default is ``None``: no organism
            restriction, matching the requirement that the dataset stays
            global and does not treat differing taxonomy as ineligibility.
        species:
            Optional human-readable label stored alongside a supplied
            ``taxonomy_id`` for readability. Purely descriptive; it does
            not select a per-species archive member (the global product
            has exactly one data file for all organisms).
        sample_limit:
            Optional development/debug row cap applied while reading the
            raw TAB3 file (``pandas.read_csv(..., nrows=sample_limit)``).
            Production default is ``None``, meaning the full selected
            BioGRID dataset is processed. The same normalization code path
            runs in both modes; sample mode only reduces input row count.
        release:
            Optional explicit BioGRID release to pin, e.g. ``"5.0.255"``
            (must match a version published under BioGRID's
            Release-Archive). When ``None`` (default), the "Latest-Release"
            alias is used and whatever concrete version it currently
            resolves to is recorded after download - see
            ``LATEST_RELEASE_URL``. Pinning a release is how a later build
            reproduces an earlier one exactly instead of always tracking
            "whatever is newest today".
        """
        super().__init__(raw_dir=raw_dir)

        if not interaction_type.strip():
            raise ValueError("interaction_type must not be empty")

        if taxonomy_id is not None and taxonomy_id <= 0:
            raise ValueError("taxonomy_id must be a positive integer")

        if species is not None and not species.strip():
            raise ValueError("species must not be empty when provided")

        if sample_limit is not None and sample_limit <= 0:
            raise ValueError("sample_limit must be a positive integer")

        if release is not None:
            if not release.strip():
                raise ValueError("release must not be empty when provided")

            if not _SAFE_RELEASE_PATTERN.match(release.strip()):
                raise ValueError(
                    f"release {release!r} is not a valid BioGRID release "
                    "identifier (expected digits/dots/hyphens only, e.g. "
                    "'5.0.255')."
                )

        self.interaction_type = interaction_type.strip().casefold()
        self.taxonomy_id = taxonomy_id
        self.species = species.strip() if species else None
        self.sample_limit = sample_limit
        self.release = release.strip() if release else None

        archive_filename = (
            f"BIOGRID-ALL-{self.release}.tab3.zip"
            if self.release is not None
            else "BIOGRID-ALL-LATEST.tab3.zip"
        )
        self.archive_path = self.raw_dir / archive_filename

        #: Concrete release version resolved from the downloaded archive's
        #: inner member name (e.g. "5.0.260"). Populated by ``_parse``.
        #: Equals ``self.release`` when a release was pinned (and is
        #: cross-checked against it there); resolved from "LATEST" only
        #: when ``self.release`` is ``None``.
        self.release_version: str | None = None

    @property
    def download_url(self) -> str:
        """Return the archive URL for the configured release.

        Single source of truth for URL construction - ``_download`` never
        builds a URL itself, so pinned-release and latest-release logic is
        defined exactly once.
        """
        if self.release is not None:
            return self._RELEASE_ARCHIVE_URL_TEMPLATE.format(
                release=self.release
            )

        return self.LATEST_RELEASE_URL

    def _download(self, force: bool = False) -> None:
        """Download the global BioGRID archive.

        Existing cached data is reused unless ``force`` is True. Downloads
        are written to a temporary ``.part`` file and moved into place only
        after the ZIP archive has been validated.
        """
        if (
            self.archive_path.exists()
            and not force
            and self._is_valid_archive(self.archive_path)
        ):
            return

        temporary_path = self.archive_path.with_suffix(
            self.archive_path.suffix + ".part"
        )

        temporary_path.unlink(missing_ok=True)

        try:
            with requests.get(
                self.download_url,
                stream=True,
                timeout=(10, 120),
            ) as response:
                response.raise_for_status()

                total_size_header = response.headers.get("content-length")
                total_size = (
                    int(total_size_header)
                    if total_size_header is not None
                    else None
                )

                with temporary_path.open("wb") as output_file:
                    with tqdm(
                        total=total_size,
                        unit="B",
                        unit_scale=True,
                        desc="Downloading BioGRID-ALL",
                    ) as progress:
                        for chunk in response.iter_content(
                            chunk_size=1024 * 1024
                        ):
                            if not chunk:
                                continue

                            output_file.write(chunk)
                            progress.update(len(chunk))

            if not self._is_valid_archive(temporary_path):
                raise RuntimeError(
                    "Downloaded BioGRID file is not a valid ZIP archive. "
                    "The URL may have changed or returned an error page."
                )

            temporary_path.replace(self.archive_path)

        except Exception:
            temporary_path.unlink(missing_ok=True)
            raise

    def _parse(self) -> pd.DataFrame:
        """Parse the global BioGRID-ALL Tab3 file."""
        if not self.archive_path.exists():
            raise FileNotFoundError(
                f"BioGRID archive does not exist: {self.archive_path}"
            )

        if not self._is_valid_archive(self.archive_path):
            raise RuntimeError(
                f"BioGRID archive is not valid: {self.archive_path}"
            )

        member_name = self._find_all_organism_file()
        self.release_version = self._extract_release_version(member_name)

        if (
            self.release is not None
            and self.release_version != self.release
        ):
            raise RuntimeError(
                f"Pinned BioGRID release {self.release!r} was requested, "
                f"but the downloaded archive resolves to release "
                f"{self.release_version!r}. Refusing to silently use a "
                "different release than the one pinned."
            )

        with ZipFile(self.archive_path) as archive:
            with archive.open(member_name) as source_file:
                data = pd.read_csv(
                    source_file,
                    sep="\t",
                    dtype="string",
                    low_memory=False,
                    nrows=self.sample_limit,
                )

        self._validate_required_columns(data)

        return data

    def _normalize(self, data: pd.DataFrame) -> pd.DataFrame:
        """Normalize BioGRID rows into a global evidence-record table."""
        interaction_type = (
            data["Experimental System Type"]
            .fillna("")
            .str.strip()
            .str.casefold()
        )

        organism_a = pd.to_numeric(
            data["Organism ID Interactor A"],
            errors="coerce",
        )
        organism_b = pd.to_numeric(
            data["Organism ID Interactor B"],
            errors="coerce",
        )

        mask = interaction_type.eq(self.interaction_type)

        if self.taxonomy_id is not None:
            mask &= organism_a.eq(self.taxonomy_id) | organism_b.eq(
                self.taxonomy_id
            )

        filtered = data.loc[mask].copy()

        publication_column = self._get_first_available_column(
            filtered,
            [
                "Publication Source",
                "Pubmed ID",
                "PubMed ID",
            ],
        )

        normalized = pd.DataFrame(
            {
                "source": pd.Series(
                    "BioGRID",
                    index=filtered.index,
                    dtype="string",
                ),
                "biogrid_release_version": pd.Series(
                    self.release_version,
                    index=filtered.index,
                    dtype="string",
                ),
                "organism_id_a": organism_a.loc[filtered.index].astype(
                    "Int64"
                ),
                "organism_id_b": organism_b.loc[filtered.index].astype(
                    "Int64"
                ),
                "organism_name_a": self._get_optional_column(
                    filtered,
                    "Organism Name Interactor A",
                ),
                "organism_name_b": self._get_optional_column(
                    filtered,
                    "Organism Name Interactor B",
                ),
                "interaction_type": filtered["Experimental System Type"],
                "experimental_system": filtered["Experimental System"],
                "biogrid_interaction_id": filtered[
                    "#BioGRID Interaction ID"
                ],
                "gene_symbol_a": filtered["Official Symbol Interactor A"],
                "gene_symbol_b": filtered["Official Symbol Interactor B"],
                "biogrid_id_a": filtered["BioGRID ID Interactor A"],
                "biogrid_id_b": filtered["BioGRID ID Interactor B"],
                "entrez_id_a": filtered["Entrez Gene Interactor A"],
                "entrez_id_b": filtered["Entrez Gene Interactor B"],
                "swissprot_accessions_a": self._get_optional_column(
                    filtered,
                    "SWISS-PROT Accessions Interactor A",
                ),
                "swissprot_accessions_b": self._get_optional_column(
                    filtered,
                    "SWISS-PROT Accessions Interactor B",
                ),
                "trembl_accessions_a": self._get_optional_column(
                    filtered,
                    "TREMBL Accessions Interactor A",
                ),
                "trembl_accessions_b": self._get_optional_column(
                    filtered,
                    "TREMBL Accessions Interactor B",
                ),
                "refseq_accessions_a": self._get_optional_column(
                    filtered,
                    "REFSEQ Accessions Interactor A",
                ),
                "refseq_accessions_b": self._get_optional_column(
                    filtered,
                    "REFSEQ Accessions Interactor B",
                ),
                "publication_id": filtered[publication_column],
                "source_database": filtered["Source Database"],
                "throughput": filtered["Throughput"],
                "score": filtered["Score"],
            }
        )

        normalized["biogrid_pair_id"] = self._build_undirected_pair_id(
            normalized["biogrid_id_a"],
            normalized["biogrid_id_b"],
        )

        normalized = normalized.drop_duplicates(
            subset=["biogrid_interaction_id"],
        )

        return normalized.reset_index(drop=True)

    def _find_all_organism_file(self) -> str:
        """Find the single global TAB3 member inside the archive."""
        with ZipFile(self.archive_path) as archive:
            candidates = [
                name
                for name in archive.namelist()
                if Path(name).name.casefold().endswith(".tab3.txt")
            ]

        if not candidates:
            raise FileNotFoundError(
                "No BIOGRID-ALL TAB3 file was found inside the archive."
            )

        if len(candidates) > 1:
            raise RuntimeError(
                "Multiple TAB3 files were found inside the BIOGRID-ALL "
                f"archive: {candidates}"
            )

        return candidates[0]

    @staticmethod
    def _extract_release_version(member_name: str) -> str:
        """Recover the concrete release version from the inner file name."""
        match = _RELEASE_VERSION_PATTERN.search(Path(member_name).name)

        if match is None:
            raise ValueError(
                "Could not determine the BioGRID release version from "
                f"archive member name {member_name!r}."
            )

        return match.group("version")

    @classmethod
    def _validate_required_columns(
        cls,
        data: pd.DataFrame,
    ) -> None:
        """Validate columns needed for basic normalization."""
        required_columns = {
            "#BioGRID Interaction ID",
            "Experimental System",
            "Experimental System Type",
            "Official Symbol Interactor A",
            "Official Symbol Interactor B",
            "BioGRID ID Interactor A",
            "BioGRID ID Interactor B",
            "Entrez Gene Interactor A",
            "Entrez Gene Interactor B",
            "Organism ID Interactor A",
            "Organism ID Interactor B",
            "Source Database",
            "Throughput",
            "Score",
        }

        missing_columns = required_columns.difference(data.columns)

        if missing_columns:
            raise KeyError(
                "BioGRID data is missing required columns: "
                f"{sorted(missing_columns)}. "
                f"Available columns: {list(data.columns)}"
            )

        cls._get_first_available_column(
            data,
            [
                "Publication Source",
                "Pubmed ID",
                "PubMed ID",
            ],
        )

    @staticmethod
    def _get_first_available_column(
        data: pd.DataFrame,
        candidates: list[str],
    ) -> str:
        """Return the first existing column from ``candidates``."""
        for column in candidates:
            if column in data.columns:
                return column

        raise KeyError(
            "None of the expected columns were found. "
            f"Expected one of: {candidates}. "
            f"Available columns: {list(data.columns)}"
        )

    @staticmethod
    def _get_optional_column(
        data: pd.DataFrame,
        column: str,
    ) -> pd.Series:
        """Return an optional BioGRID column or an empty string Series."""
        if column in data.columns:
            return data[column]

        return pd.Series(
            pd.NA,
            index=data.index,
            dtype="string",
        )

    @staticmethod
    def _build_undirected_pair_id(
        id_a: pd.Series,
        id_b: pd.Series,
    ) -> pd.Series:
        """Create a source-specific undirected BioGRID pair key."""
        a = id_a.astype("string")
        b = id_b.astype("string")

        a_missing = a.isna()
        b_missing = b.isna()

        first = a.where(a <= b, b)
        second = b.where(a <= b, a)

        pair_id = first + "__" + second

        return pair_id.mask(a_missing | b_missing)

    @staticmethod
    def _is_valid_archive(path: Path) -> bool:
        """Return True when ``path`` is a readable, undamaged ZIP file."""
        if not path.is_file():
            return False

        try:
            with ZipFile(path) as archive:
                return archive.testzip() is None
        except (BadZipFile, OSError):
            return False

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"interaction_type={self.interaction_type!r}, "
            f"taxonomy_id={self.taxonomy_id!r}, "
            f"species={self.species!r}, "
            f"sample_limit={self.sample_limit!r}, "
            f"release={self.release!r}, "
            f"raw_dir={str(self.raw_dir)!r})"
        )
