from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.common import normalize_contig


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("chr1", "1"),
        ("Chr01", "1"),
        ("chr0007", "7"),
        ("1", "1"),
        ("01", "1"),
        ("chrM", "m"),
        (" ChrX ", "x"),
        ("000", "0"),
        ("chr000", "0"),
    ],
)
def test_normalize_contig_handles_common_aliases(raw: str, expected: str) -> None:
    assert normalize_contig(raw) == expected
