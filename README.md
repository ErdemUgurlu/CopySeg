# CopySeg

A C++ based k-mer preprocessing and Python (HMM) based segmentation tool.

## 🛠 Usage

### 1. Preprocessing (C++)
To clean and filter k-mer data for noise reduction:
```bash
g++ -O3 preprocess_kmers_v5.cpp -o preprocess
./preprocess input.kmers output.clean
