
from pathlib import Path
import random

random.seed(42)

outdir = Path("benchmark_tfasty")
outdir.mkdir(exist_ok=True)

codons = {
    "A": "GCT", "C": "TGT", "D": "GAT",
    "E": "GAA", "F": "TTT", "G": "GGT",
    "H": "CAT", "I": "ATT", "K": "AAA",
    "L": "CTG", "M": "ATG", "N": "AAT",
    "P": "CCT", "Q": "CAA", "R": "CGT",
    "S": "TCT", "T": "ACT", "V": "GTT",
    "W": "TGG", "Y": "TAT"
}

protein = "M" + "".join(
    random.choices(list(codons), k=199)
)

cds = "".join(codons[aa] for aa in protein) + "TAA"

# Mutation positions (zero-based nucleotide coordinates)
p1 = 180
p2 = 420

variants = {
    "intact": cds,

    "insertion_1nt":
        cds[:p1] + "A" + cds[p1:],

    "deletion_1nt":
        cds[:p1] + cds[p1 + 1:],

    "deletion_3nt":
        cds[:p1] + cds[p1 + 3:],

    "premature_stop":
        cds[:p1] + "TGA" + cds[p1 + 3:],

    "missing_start":
        "GCG" + cds[3:],

    "missing_stop":
        cds[:-3] + "CAA",

    "compensatory_frameshifts":
        cds[:p1] + "A" +
        cds[p1:p2] + cds[p2 + 1:],
}

with (outdir / "reference.faa").open("w") as f:
    f.write(f">reference_protein\n{protein}\n")

with (outdir / "variants.fna").open("w") as f:
    for name, seq in variants.items():
        f.write(f">{name}\n{seq}\n")

for name, seq in variants.items():
    with (outdir / f"{name}.fna").open("w") as f:
        f.write(f">{name}\n{seq}\n")

    print(f"{name:28s} {len(seq):4d} nt")
