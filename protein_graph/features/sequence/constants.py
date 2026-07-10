"""Shared amino acid groups for sequence features."""

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

HYDROPHOBIC = set("AVILMFWY")
HYDROPHILIC = set("RNDQEHKSTYC")
CHARGED = set("DEKRH")
POSITIVE = set("KRH")
NEGATIVE = set("DE")
POLAR = set("STNQCY")
AROMATIC = set("FWY")

GLOBAL_GROUPS = {
    "hydrophobic": HYDROPHOBIC,
    "hydrophilic": HYDROPHILIC,
    "charged": CHARGED,
    "positive": POSITIVE,
    "negative": NEGATIVE,
    "polar": POLAR,
    "aromatic": AROMATIC,
}

LOCAL_GROUPS = {
    "hydrophobic": HYDROPHOBIC,
    "charged": CHARGED,
    "positive": POSITIVE,
    "negative": NEGATIVE,
}

WINDOW_SIZES = [5, 10, 20, 30]
RUN_BINS = [(0, 5), (5, 10), (10, 15), (15, 20), (20, 25), (25, 30), (30, 35), (35, 40)]
PERIODS = [2, 3, 4, 7]
