
import json
import subprocess
from pathlib import Path

GENOME = "GCA_000008865.2"

INPUT = Path(
    f"results_nf/{GENOME}/pseudogenes/"
    f"{GENOME}.tfasty_test.jsonl"
)

REFERENCE = Path("uniprot_sprot.fasta")

TFASTY = Path(
    "fasta-36.3.8i/bin/tfasty36"
).resolve()

OUTPUT = Path("validation_tfasty") / GENOME
OUTPUT.mkdir(parents=True, exist_ok=True)

SELECTED = {
    "region_3",
    "region_34",
    "region_135",
    "region_136",
    "region_151",
    "region_159",
    "region_210",
    "region_263",
    "region_289",
}

records = []

with INPUT.open() as f:
    for line in f:
        record = json.loads(line)

        if record["region"]["region_id"] in SELECTED:
            records.append(record)

wanted = {
    r["best_protein"]["protein_id"]
    for r in records
}

proteins = {}

# Retrieve only the required Swiss-Prot sequences.
with REFERENCE.open() as f:
    protein_id = None
    chunks = []

    for line in f:
        if line.startswith(">"):

            if protein_id in wanted:
                proteins[protein_id] = "".join(chunks)

            protein_id = line[1:].split()[0]
            chunks = []

        elif protein_id in wanted:
            chunks.append(line.strip())

    if protein_id in wanted:
        proteins[protein_id] = "".join(chunks)

for record in records:

    region = record["region"]
    region_id = region["region_id"]
    protein_id = record["best_protein"]["protein_id"]

    if protein_id not in proteins:
        print(f"Missing reference: {protein_id}")
        continue

    directory = OUTPUT / region_id
    directory.mkdir(exist_ok=True)

    faa = directory / "reference.faa"
    fna = directory / "candidate.fna"

    faa.write_text(
        f">{protein_id}\n"
        f"{proteins[protein_id]}\n"
    )

    fna.write_text(
        f">{region_id}\n"
        f"{record['sequence_nt']}\n"
    )

    (directory / "annotation.json").write_text(
        json.dumps(record, indent=2)
    )

    for mode in ["0", "10", "8CC"]:

        output_file = directory / f"alignment_m{mode}.txt"

        command = [
            str(TFASTY),
            "-3",
            "-m", mode,
            str(faa),
            str(fna),
        ]

        with output_file.open("w") as handle:
            subprocess.run(
                command,
                stdout=handle,
                check=True
            )

    print(
        region_id,
        protein_id,
        region["chromosome"],
        region["start"],
        region["end"],
        region["strand"],
        sep="\t"
    )

print(f"\nValidation files saved in {OUTPUT}")
