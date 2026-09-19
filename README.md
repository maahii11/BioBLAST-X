# 🧬 BioBLAST-X

## BLAST-Inspired Biological Sequence Similarity Analyzer

BioBLAST-X is a Python-based bioinformatics application designed to analyze DNA sequences using k-mer based similarity screening and sequence comparison.

The application provides an interactive Streamlit interface where users can enter a DNA sequence or upload a FASTA file and perform sequence similarity analysis.

---

## 🎯 Objective

The main objective of BioBLAST-X is to provide a simple and interactive platform for understanding biological sequence similarity.

The application follows a BLAST-inspired workflow:

**Sequence Input → Validation → K-mer Screening → Similarity Search → Sequence Alignment → Results & Visualization**

---

## 🚀 Features

- DNA sequence input through text
- FASTA file upload
- DNA sequence validation
- Sequence length calculation
- GC-content calculation
- Configurable k-mer size
- K-mer based similarity screening
- Sequence similarity comparison
- Similarity ranking
- Local sequence alignment
- Similarity visualization
- Top similarity match identification
- CSV result export
- Interactive Streamlit interface

---

## 🧪 Methodology

### 1. Sequence Input

The user can either:

- Paste a DNA sequence directly into the application
- Upload a FASTA file

### 2. Sequence Validation

The input sequence is checked to ensure that it contains valid DNA nucleotide characters.

### 3. Sequence Statistics

The application calculates basic sequence properties such as:

- Sequence length
- GC content
- Selected k-mer size

### 4. K-mer Screening

The DNA sequence is divided into overlapping k-mers.

These k-mers are compared with sequences stored in the local database to identify potentially similar sequences.

### 5. Similarity Analysis

Candidate sequences are compared with the query sequence and similarity values are calculated.

### 6. Sequence Alignment

The selected similar sequence is aligned with the query sequence to visualize their correspondence.

### 7. Results and Visualization

The application displays:

- Similarity ranking
- Sequence identity
- Shared k-mers
- Top similarity match
- Sequence alignment
- Similarity visualization

Results can also be downloaded as a CSV file.

---

## 🛠️ Technologies Used

- **Python**
- **Streamlit**
- **Bioinformatics**
- **FASTA sequence format**
- **K-mer analysis**
- **Sequence alignment**
- **CSV data processing**

---

## 📂 Project Structure

```text
BioBLAST-X/
│
├── app.py
├── database.py
├── search_engine.py
├── alignment.py
├── requirements.txt
│
└── data/
    └── sequences.fasta