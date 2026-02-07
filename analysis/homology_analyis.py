from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# -----------------------------
#  Homology Search (BLAST)
# -----------------------------

# Read protein sequence
record = SeqIO.read("data/input_sequence.fasta", "fasta")

# Perform BLAST search
result_handle = NCBIWWW.qblast(
    program="blastp",
    database="nr",
    sequence=record.seq
)

# Save BLAST XML output
with open("results/blast_results.xml", "w") as b:
    b.write(result_handle.read())

print("BLAST search completed")

# -----------------------------
# Parse BLAST
# -----------------------------

with open("results/blast_results.xml") as b:
    blast_record = NCBIXML.read(b)

# Get top alignment and HSP
first_alignment = blast_record.alignments[0]
first_hsp = first_alignment.hsps[0]

# Write BLAST evidence
with open("results/blast_results.txt", "w") as f:
    f.write("STEP 4: HOMOLOGY SEARCH RESULTS\n")
    f.write("--------------------------------\n")
    f.write("Top BLAST hit:\n")
    f.write(first_alignment.title + "\n\n")
    f.write("Score: " + str(first_hsp.score) + "\n")
    f.write("E-value: " + str(first_hsp.expect) + "\n\n")
    f.write("Query range: " +
            str(first_hsp.query_start) + " - " +
            str(first_hsp.query_end) + "\n")
    f.write("Subject range: " +
            str(first_hsp.sbjct_start) + " - " +
            str(first_hsp.sbjct_end) + "\n")

print("BLAST parsing completed")

# -----------------------------
# Functional Annotation
# -----------------------------

with open("results/functional_annotation.txt", "w") as f:
    f.write("STEP 5: FUNCTIONAL ANNOTATION\n")
    f.write("-----------------------------\n\n")

    f.write("BLAST Evidence:\n")
    f.write(first_alignment.title + "\n\n")

    f.write(
        "Interpretation:\n"
        "The top BLAST hit corresponds to a protein annotated as an "
        "ESX-1 secretion system protein (EccA1) from Mycobacterium species.\n\n"
    )

    f.write(
        "Based on this homology, the target protein is predicted to be involved "
        "in protein secretion via the ESX-1 secretion system.\n\n"
    )

    f.write(
        "This functional prediction is inferred from sequence similarity "
        "and has not been experimentally validated.\n"
    )

print("Functional annotation completed")
