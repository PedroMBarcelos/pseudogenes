import pandas as pd

df = pd.read_csv("annotation_audit_156_v2.tsv", sep="\t")

gff = df[df["status"] == "PASS"].copy()
pseudo = gff[gff["pseudogene"] > 0].copy()

print("Genomes with valid GFFs:", len(gff))
print("Genomes with pseudogenes:", len(pseudo))
print("Total annotated pseudogenes:", int(pseudo["pseudogene"].sum()))

comparisons = {
    "CDS masking": "pseudo_by_CDS_any",
    "Protein-coding gene masking": "pseudo_by_protein_gene_any",
}

for label, prefix in comparisons.items():
    full = int(pseudo[f"{prefix}_full"].sum())
    partial = int(pseudo[f"{prefix}_partial"].sum())
    none = int(pseudo[f"{prefix}_none"].sum())

    print(f"\n{label}")
    print("  Fully covered pseudogenes:", full)
    print("  Partially covered:", partial)
    print("  No overlap:", none)
    print(
        "  Genomes with pseudogenes affected:",
        int((
            pseudo[f"{prefix}_full"] +
            pseudo[f"{prefix}_partial"] > 0
        ).sum())
    )

# How effectively would gene-only masking exclude CDS intervals?
print("\nCDS coverage by protein-coding genes:")
print(
    "  Fully covered:",
    int(gff["CDS_by_protein_gene_any_full"].sum())
)
print(
    "  Partially covered:",
    int(gff["CDS_by_protein_gene_any_partial"].sum())
)
print(
    "  Not covered:",
    int(gff["CDS_by_protein_gene_any_none"].sum())
)