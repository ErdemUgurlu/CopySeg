# CopySeg

**CopySeg** is a high-performance C++ bioinformatics HMM based tool designed for **Copy Number Variation (CNV)** analysis using high-noise k-mer data.

## 🛠 Usage

### 1. Preprocessing (C++)
To run the tool:
```bash
g++ -O3 preprocess_kmers_v5.cpp -o preprocess
./preprocess input.kmers output.clean
