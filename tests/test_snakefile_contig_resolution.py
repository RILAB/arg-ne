import ast
import re
from pathlib import Path


def _extract_function_source(text: str, name: str) -> str:
    lines = text.splitlines(keepends=True)
    start = None
    for i, line in enumerate(lines):
        if line.startswith(f"def {name}("):
            start = i
            break
    if start is None:
        raise AssertionError(f"Function {name} not found in Snakefile")

    end = len(lines)
    i = start
    paren_balance = 0
    while i < len(lines):
        paren_balance += lines[i].count("(") - lines[i].count(")")
        if paren_balance <= 0 and lines[i].rstrip().endswith(":"):
            i += 1
            break
        i += 1

    for i in range(i, len(lines)):
        line = lines[i]
        if not line.strip():
            continue
        if line.startswith(" ") or line.startswith("\t"):
            continue
        if line.startswith("#"):
            continue
        if line.startswith("def "):
            end = i
            break
        end = i
        break
    return "".join(lines[start:end])


def _load_snakefile_functions():
    snakefile = Path(__file__).resolve().parents[1] / "Snakefile"
    text = snakefile.read_text(encoding="utf-8")
    snippet = (
        _extract_function_source(text, "_normalize_contig")
        + "\n"
        + _extract_function_source(text, "_resolve_requested_contigs")
    )
    module = ast.parse(snippet, filename=str(snakefile))
    ast.fix_missing_locations(module)
    ns = {"re": re}
    exec(compile(module, str(snakefile), "exec"), ns)
    return ns["_normalize_contig"], ns["_resolve_requested_contigs"]


def test_normalize_contig_handles_chr_prefix_and_zero_padding():
    normalize, _ = _load_snakefile_functions()
    assert normalize("chr01") == "1"
    assert normalize("001") == "1"
    assert normalize("chrX") == "x"
    assert normalize("chr000") == "0"


def test_resolve_requested_contigs_remaps_dedupes_and_drops_ambiguous():
    _, resolve = _load_snakefile_functions()

    kept, dropped, remapped = resolve(["chr1", "01", "2", "chr2", "chr3"], ["1", "2"])
    assert kept == ["1", "2"]
    assert dropped == ["chr3"]
    assert remapped == [("chr1", "1"), ("01", "1"), ("chr2", "2")]

    kept2, dropped2, remapped2 = resolve(["chr1"], ["1", "chr01"])
    assert kept2 == []
    assert dropped2 == ["chr1"]
    assert remapped2 == []
