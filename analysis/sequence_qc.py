from Bio import SeqIO

# Read sequence
record = SeqIO.read("data/input_sequence.fasta", "fasta")
sequence = record.seq

# Basic analysis
length = len(sequence)

aa_composition = {}
for aa in sequence:
    if aa in aa_composition:
        aa_composition[aa] += 1
    else:
        aa_composition[aa] = 1

# Validation decision
if length < 100:
    decision = "Sequence rejected (too short for analysis)"
else:
    decision = "Sequence accepted for downstream analysis"

# Print results
print("\nSEQUENCE QUALITY ANALYSIS")
print("------------------------")
print("Sequence ID:", record.id)
print("Sequence Length:", length)
print("\nAmino Acid Composition:")
for aa in aa_composition:
    print(aa, ":", aa_composition[aa])
print("\nDecision:", decision)

# Write results
with open("results/qc_summary.txt", "w") as f:
    f.write("SEQUENCE QUALITY ANALYSIS\n")
    f.write("------------------------\n")
    f.write("Sequence ID: " + record.id + "\n")
    f.write("Sequence Length: " + str(length) + "\n\n")

    f.write("Amino Acid Composition:\n")
    for aa in aa_composition:
        f.write(aa + " : " + str(aa_composition[aa]) + "\n")

    f.write("\nDecision:\n")
    f.write(decision + "\n")

