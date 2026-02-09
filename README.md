# CopySeg

**CopySeg** is a high-performance C++ bioinformatics HMM based tool designed for **Copy Number Variation (CNV)** analysis using high-noise k-mer data.

## 🛠 Usage

### 1. Preprocessing (C++)# CopySeg

**CopySeg** is a high-performance C++ bioinformatics HMM based tool designed for **Copy Number Variation (CNV)** analysis using high-noise k-mer data.

## 🛠 Usage



First, compile and run the preprocessing script to clean the k-mer data:

```bash
g++ -O3 preprocess_kmers_v5.cpp -o preprocess
./preprocess input.kmers output.clean

Then, run the 7-state Hidden Markov Model script to detect CNV segments:

python3 segment_cnv_hmm_log_7state.py output.clean
