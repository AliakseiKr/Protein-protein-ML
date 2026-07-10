from pathlib import Path
from zipfile import BadZipFile, ZipFile

import pandas as pd
import requests
from tqdm import tqdm

from protein_graph.data_sources.base import BaseDataSource


class BioGRIDDataSource(BaseDataSource):
    """BioGRID interaction evidence data source.

    The normalized output contains individual BioGRID evidence records.
    It does not resolve identifiers into canonical ProteinGraph protein IDs
    and does not aggregate multiple evidence rows into unique protein pairs.
    """

    DOWNLOAD_URL = (
        "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/"
        "BIOGRID-ORGANISM-LATEST.tab3.zip"
    )

    def __init__(
        self,
        species: str = "Rattus norvegicus",
        taxonomy_id: int = 10116,
        interaction_type: str = "physical",
        raw_dir: str | Path = "data/raw/biogrid",
    ) -> None:
        super().__init__(raw_dir=raw_dir)

        if not species.strip():
            raise ValueError("species must not be empty")

        if taxonomy_id <= 0:
            raise ValueError("taxonomy_id must be a positive integer")

        if not interaction_type.strip():
            raise ValueError("interaction_type must not be empty")

        self.species = species.strip()
        self.taxonomy_id = taxonomy_id
        self.interaction_type = interaction_type.strip().casefold()

        self.archive_path = (
            self.raw_dir / "BIOGRID-ORGANISM-LATEST.tab3.zip"
        )

    def _download(self, force: bool = False) -> None:
        """Download the BioGRID organism archive.

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
                self.DOWNLOAD_URL,
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
                        desc="Downloading BioGRID",
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
        """Parse the organism-specific BioGRID Tab3 file."""
        if not self.archive_path.exists():
            raise FileNotFoundError(
                f"BioGRID archive does not exist: {self.archive_path}"
            )

        if not self._is_valid_archive(self.archive_path):
            raise RuntimeError(
                f"BioGRID archive is not valid: {self.archive_path}"
            )

        species_file = self._find_species_file()

        with ZipFile(self.archive_path) as archive:
            with archive.open(species_file) as source_file:
                data = pd.read_csv(
                    source_file,
                    sep="\t",
                    dtype="string",
                    low_memory=False,
                )

        self._validate_required_columns(data)

        return data

    def _normalize(self, data: pd.DataFrame) -> pd.DataFrame:
        """Normalize BioGRID rows into an evidence-record table."""
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

        mask = (
            interaction_type.eq(self.interaction_type)
            & organism_a.eq(self.taxonomy_id)
            & organism_b.eq(self.taxonomy_id)
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
                "species": pd.Series(
                    self.species,
                    index=filtered.index,
                    dtype="string",
                ),
                "taxonomy_id": pd.Series(
                    self.taxonomy_id,
                    index=filtered.index,
                    dtype="Int64",
                ),
                "interaction_type": filtered[
                    "Experimental System Type"
                ],
                "experimental_system": filtered[
                    "Experimental System"
                ],
                "biogrid_interaction_id": filtered[
                    "#BioGRID Interaction ID"
                ],
                "gene_symbol_a": filtered[
                    "Official Symbol Interactor A"
                ],
                "gene_symbol_b": filtered[
                    "Official Symbol Interactor B"
                ],
                "biogrid_id_a": filtered[
                    "BioGRID ID Interactor A"
                ],
                "biogrid_id_b": filtered[
                    "BioGRID ID Interactor B"
                ],
                "entrez_id_a": filtered[
                    "Entrez Gene Interactor A"
                ],
                "entrez_id_b": filtered[
                    "Entrez Gene Interactor B"
                ],
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

    def _find_species_file(self) -> str:
        """Find the archive member matching the selected species."""
        species_token = self.species.replace(" ", "_").casefold()

        with ZipFile(self.archive_path) as archive:
            candidates = [
                name
                for name in archive.namelist()
                if species_token in Path(name).name.casefold()
                and name.casefold().endswith(".txt")
            ]

        if not candidates:
            raise FileNotFoundError(
                "No BioGRID file was found for species "
                f"{self.species!r}."
            )

        if len(candidates) > 1:
            raise RuntimeError(
                "Multiple BioGRID files were found for species "
                f"{self.species!r}: {candidates}"
            )

        return candidates[0]

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
            f"species={self.species!r}, "
            f"taxonomy_id={self.taxonomy_id!r}, "
            f"interaction_type={self.interaction_type!r}, "
            f"raw_dir={str(self.raw_dir)!r})"
        )
