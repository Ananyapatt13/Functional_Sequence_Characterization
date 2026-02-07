from Bio import SeqIO
from Bio.Blast import NCBIXML
import matplotlib.pyplot as plt
import re

# -----------------------------
# Read protein sequence
# -----------------------------
record = SeqIO.read("data/input_sequence.fasta", "fasta")
seq = str(record.seq)
length = len(seq)

# -----------------------------
# RESULT 1: Motif Analysis
# -----------------------------
walker_a_pattern = r"G....GK[ST]"
walker_b_pattern = r"[ILVFMW][ILVFMW][ILVFMW][ILVFMW]DE"

walker_a_matches = list(re.finditer(walker_a_pattern, seq))
walker_b_matches = list(re.finditer(walker_b_pattern, seq))

with open("results/motif_analysis.txt", "w") as f:
    f.write("MOTIF ANALYSIS RESULTS\n")
    f.write("---------------------\n\n")

    f.write("Walker A motifs:\n")
    for m in walker_a_matches:
        f.write(f"Position {m.start()+1}-{m.end()}: {m.group()}\n")

    f.write("\nWalker B motifs:\n")
    for m in walker_b_matches:
        f.write(f"Position {m.start()+1}-{m.end()}: {m.group()}\n")

print("Motif analysis completed")

# -----------------------------
# Read BLAST XML
# -----------------------------
with open("results/blast_results.xml") as b:
    blast_record = NCBIXML.read(b)

hsp = blast_record.alignments[0].hsps[0]

# -----------------------------
# RESULT 2: Pairwise Identity
# -----------------------------
identity_percent = (hsp.identities / hsp.align_length) * 100

with open("results/identity_analysis.txt", "w") as f:
    f.write("PAIRWISE IDENTITY ANALYSIS\n")
    f.write("--------------------------\n\n")
    f.write(f"Identities: {hsp.identities}\n")
    f.write(f"Alignment length: {hsp.align_length}\n")
    f.write(f"Percent identity: {identity_percent:.2f}%\n")

print("Identity analysis completed")

# -----------------------------
# PLOT 1: Amino Acid Composition
# -----------------------------
aa_counts = {}
for aa in seq:
    aa_counts[aa] = aa_counts.get(aa, 0) + 1

plt.figure()
plt.bar(aa_counts.keys(), aa_counts.values())
plt.xlabel("Amino Acid")
plt.ylabel("Count")
plt.title("Amino Acid Composition of EccA1")
plt.tight_layout()
plt.savefig("results/aa_composition.png")
plt.close()

# -----------------------------
# PLOT 2: Conservation Profile
# -----------------------------
conservation = [1 if c == "|" else 0 for c in hsp.match]

plt.figure()
plt.plot(conservation)
plt.xlabel("Alignment Position")
plt.ylabel("Conservation (1 = conserved)")
plt.title("Conservation Profile of EccA1")
plt.tight_layout()
plt.savefig("results/conservation_plot.png")
plt.close()

# -----------------------------
# PLOT 3: Identity Summary
# -----------------------------
plt.figure()
plt.bar(["EccA1 vs Closest Homolog"], [identity_percent])
plt.ylabel("Percent Identity")
plt.ylim(0, 100)
plt.title("Pairwise Identity with Closest Homolog")
plt.tight_layout()
plt.savefig("results/identity_plot.png")
plt.close()

print("Plots generated successfully")
