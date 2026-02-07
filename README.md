# Functional Sequence Characterization using Biopython

## Project Overview

This project implements a stepwise bioinformatics pipeline for the
functional characterization of a protein sequence using **Biopython**.
The workflow follows a homology-based approach, where protein function
is inferred from sequence similarity rather than experimentally proven.

The target sequence analysed in this project corresponds to the
**ESX-1 secretion system protein EccA1** from *Mycobacterium tuberculosis*.

---

## Project Objectives

- To perform basic sequence quality assessment and validation
- To identify homologous proteins using BLAST
- To infer protein function based on sequence similarity
- To strengthen functional inference using additional computational analyses
- To provide a biologically meaningful interpretation of results

---

## Project Structure

Functional_Sequence_Characterization/
│
├── analysis/
│ ├── sequence_qc.py
│ ├── homology_analysis.py
│ └── additional_analysis.py
│
├── data/
│ └── input_sequence.fasta
│
├── results/
│ ├── qc_summary.txt
│ ├── blast_results.xml
│ ├── blast_results.txt
│ ├── functional_annotation.txt
│ ├── motif_analysis.txt
│ ├── identity_analysis.txt
│ ├── biological_interpretation.txt
│ ├── aa_composition.png
│ ├── conservation_plot.png
│ └── identity_plot.png
│
└── README.md
---

## Pipeline Description

### Step 2: Sequence Quality & Basic Analysis
- Reads the protein sequence using Biopython
- Calculates sequence length
- Computes amino acid composition
- Confirms suitability for downstream analysis

### Step 3: Sequence Filtering & Validation
- Applies a length-based biological criterion
- Accepts or rejects the sequence for further analysis

### Step 4: Homology Search (BLAST)
- Performs similarity search using BLAST (blastp)
- Identifies closest homologs
- Examines conserved regions and evolutionary evidence
- Extracts alignment score, E-value, and alignment ranges

### Step 5: Functional Annotation
- Uses the top BLAST hit as evidence
- Infers protein function from known annotations
- Emphasizes that functional prediction is homology-based

### Step 6: Biological Interpretation (Final Step)
- Integrates results from all analyses
- Interprets biological role and organism relevance
- Clearly states that conclusions are computationally inferred

---

## Additional Computational Analyses

To strengthen functional inference, the following code-based analyses were performed:

- **Motif Analysis**
  - Detection of conserved Walker A and Walker B motifs
  - Supports ATPase activity of the protein

- **Pairwise Identity Analysis**
  - Quantitative measurement of sequence similarity
  - Provides numerical evidence for strong homology

- **Visualization**
  - Amino acid composition bar plot
  - Conservation profile plot from BLAST alignment
  - Percent identity summary plot

These analyses complement the core homology-based pipeline.

---

## Key Findings

- The protein is highly conserved across *Mycobacterium* species
- BLAST analysis shows a full-length alignment with E-value 0.0
- Conserved ATPase motifs were identified
- Strong evidence supports the protein being an ATPase component of the
  ESX-1 (Type VII) secretion system
- The ESX-1 system is known to play a role in virulence of
  *Mycobacterium tuberculosis*

---

## Important Note

All functional conclusions presented in this project are **inferred**
from sequence-based computational analyses and **have not been
experimentally validated**.

---

## Requirements

- Python
- Biopython
- matplotlib

---

## How to Run

From the project root directory:

```bash
python analysis/sequence_qc.py
python analysis/homology_analysis.py
python analysis/additional_analysis.py

