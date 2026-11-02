"""Which 10x chemistry Cell Ranger found in a pool, and whether HLA naming can work on it.

HLA alleles differ mostly in exons 2 and 3, near the 5' end of the genes. 5'
libraries read that end; 3' libraries read the other end, where the alleles hardly
differ (in a 10x 3' test library, under 1% of the HLA-A reads covered exons 1-3,
against 84-92% in 5' libraries). So with 3' data souporcell still splits the
patients, but most clusters cannot be named from the HLA table.
"""

import glob
import re
from pathlib import Path

THREE_PRIME_NOTE = (
    "HLA naming needs 5' data: 3' reads cover the end of the HLA genes where the alleles "
    "hardly differ, so clusters often stay unmatched. The clusters themselves (souporcell) "
    "are not affected; name them with a manual_assignments file if needed."
)

_CHEMISTRY = re.compile(r'"Chemistry","([^"]+)"')


def gene_expression_chemistry(cellranger_outs):
    """Cell Ranger's name for the gene expression chemistry, e.g. "Single Cell 5' R2-only v3",
    read from its web summary; None when it cannot be found."""
    outs = Path(cellranger_outs)
    paths = sorted(glob.glob(str(outs / "per_sample_outs" / "*" / "web_summary.html")))
    paths.append(str(outs / "web_summary.html"))   # cellranger count
    for path in paths:
        try:
            text = Path(path).read_text(errors="replace")
        except OSError:
            continue
        for name in _CHEMISTRY.findall(text):
            name = name.replace("\\u0027", "'")
            if "V(D)J" not in name:   # the TCR library has its own chemistry
                return name
    return None


def is_five_prime(chemistry):
    """True for 5' gene expression, False for anything else (3', Flex...), None if unknown."""
    if chemistry is None:
        return None
    return "5'" in chemistry
