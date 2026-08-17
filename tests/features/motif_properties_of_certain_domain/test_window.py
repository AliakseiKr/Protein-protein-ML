import pytest

from protein_graph.features.motif_properties_of_certain_domain.window import (
    WindowInterval,
)


def test_from_biological_converts_to_zero_based_half_open() -> None:
    assert WindowInterval.from_biological(3, 6) == WindowInterval(2, 6)


def test_overlapping_windows_overlap() -> None:
    assert WindowInterval(2, 6).overlaps(WindowInterval(5, 8)) is True


def test_adjacent_windows_do_not_overlap() -> None:
    assert WindowInterval(2, 6).overlaps(WindowInterval(6, 8)) is False


def test_from_biological_rejects_start_below_one() -> None:
    with pytest.raises(ValueError, match="biological start must be >= 1"):
        WindowInterval.from_biological(0, 5)


def test_rejects_non_positive_length_window() -> None:
    with pytest.raises(ValueError, match="end must be greater than start"):
        WindowInterval(3, 3)
