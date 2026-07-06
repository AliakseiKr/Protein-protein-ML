from pathlib import Path
from zipfile import BadZipFile, ZipFile

import pandas as pd
import requests
from tqdm import tqdm

from protein_graph.data_sources.base import BaseDataSource


class BioGRIDDataSource(BaseDataSource):
    """BioGRID interaction data source."""

    DOWNLOAD_URL = (
        "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/"
        "BIOGRID-ORGANISM-LATEST.tab3.zip"
    )

    def __init__(
        self,
        species: str = "Rattus norvegicus",
        interaction_type: str = "physical",
        raw_dir: str | Path = "data/raw/biogrid",
    ) -> None:
        super().__init__(raw_dir=raw_dir)
        self.species = species
        self.interaction_type = interaction_type.lower()
        self.archive_path = self.raw_dir / "BIOGRID-ORGANISM-LATEST.tab3.zip"

    def _download(self) -> None:
        """Download the latest BioGRID organism archive if missing."""
        if self.archive_path.exists():
            return

        response = requests.get(self.DOWNLOAD_URL, stream=True, timeout=60)
        response.raise_for_status()

        total_size = int(response.headers.get("content-length", 0))

        with open(self.archive_path, "wb") as file:
            with tqdm(
                total=total_size,
                unit="B",
                unit_scale=True,
                desc="Downloading BioGRID",
            ) as progress:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        file.write(chunk)
                        progress.update(len(chunk))

        try:
            with ZipFile(self.archive_path):
                pass
        except BadZipFile as error:
            self.archive_path.unlink(missing_ok=True)
            raise RuntimeError(
                "Downloaded BioGRID file is not a valid ZIP archive. "
                "The download URL may have changed or returned an HTML error page."
            ) from error

    def _parse(self) -> pd.DataFrame:
        """Parse the BioGRID Tab3 file for the selected species."""
        species_file = self._find_species_file()

        with ZipFile(self.archive_path) as archive:
            with archive.open(species_file) as file:
                data = pd.read_csv(file, sep="\t", low_memory=False)

        return data

    def _normalize(self, data: pd.DataFrame) -> pd.DataFrame:
        """Normalize BioGRID interactions into ProteinGraph format."""
        filtered = data.copy()

        filtered = filtered[
            filtered["Experimental System Type"].str.lower()
            == self.interaction_type
        ].copy()

        publication_id_column = self._get_first_available_column(
            filtered,
            ["Publication Source", "Pubmed ID", "PubMed ID"],
        )

        normalized = pd.DataFrame(
            {
                "source": "BioGRID",
                "species": self.species,
                "interaction_type": filtered["Experimental System Type"],
                "experimental_system": filtered["Experimental System"],
                "biogrid_interaction_id": filtered["#BioGRID Interaction ID"],
                "gene_a": filtered["Official Symbol Interactor A"],
                "gene_b": filtered["Official Symbol Interactor B"],
                "biogrid_id_a": filtered["BioGRID ID Interactor A"],
                "biogrid_id_b": filtered["BioGRID ID Interactor B"],
                "entrez_id_a": filtered["Entrez Gene Interactor A"],
                "entrez_id_b": filtered["Entrez Gene Interactor B"],
                "publication_id": filtered[publication_id_column],
                "source_database": filtered["Source Database"],
                "throughput": filtered["Throughput"],
                "score": filtered["Score"],
            }
        )

        normalized = normalized.drop_duplicates(
            subset=["gene_a", "gene_b", "experimental_system", "publication_id"]
        ).reset_index(drop=True)

        return normalized

    def _find_species_file(self) -> str:
        """Find the BioGRID file inside the archive matching the target species."""
        species_token = self.species.replace(" ", "_")

        with ZipFile(self.archive_path) as archive:
            candidates = [
                name
                for name in archive.namelist()
                if species_token in name and name.endswith(".txt")
            ]

        if not candidates:
            raise FileNotFoundError(
                f"No BioGRID file found for species: {self.species}"
            )

        if len(candidates) > 1:
            raise RuntimeError(
                f"Multiple BioGRID files found for species {self.species}: "
                f"{candidates}"
            )

        return candidates[0]

    @staticmethod
    def _get_first_available_column(
        data: pd.DataFrame,
        candidates: list[str],
    ) -> str:
        """Return the first available column name from a list of candidates."""
        for column in candidates:
            if column in data.columns:
                return column

        raise KeyError(
            "None of the expected columns were found. "
            f"Expected one of: {candidates}. "
            f"Available columns: {list(data.columns)}"
        )

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"species={self.species!r}, "
            f"interaction_type={self.interaction_type!r}, "
            f"raw_dir={str(self.raw_dir)!r})"
        )
