#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <random>

// === CONFIGURATION ===

// Window size for k-mer aggregation
const int WINDOW_SIZE = 1000;

// Saturation Handling (PacBio HiFi dynamic range limit)
const int SATURATION_THRESHOLD = 65535;
const double SATURATION_FACTOR = 1.7;

// Savitzky-Golay Smoothing Parameters
const int SMOOTH_WINDOW = 5;
const double SG_COEFFS[] = {-3.0/35.0, 12.0/35.0, 17.0/35.0, 12.0/35.0, -3.0/35.0};

// K-mer Filtering & Clamping Thresholds (Option C)
const int MIN_COUNT_THRESHOLD = 10;
const int MAX_COUNT_THRESHOLD_BASE = 350;
const int EXTREME_REPEAT_ABSOLUTE = 1200;

// Window-based decision thresholds
const double UNIFORM_CV_THRESHOLD = 0.50;
const double DUPLICATION_MEAN_FACTOR = 1.8;
const double MAX_BIOLOGICAL_CN_FACTOR = 200.0;

// Visualization Config
const int PLOTS_PER_CATEGORY = 5; // How many samples to save for each category

// === DATA STRUCTURES ===

struct KmerData {
    int position;
    int count;
    double adjusted_count;
};

struct WindowStats {
    double winsorized_mean;
    double iqr;
    double robust_cv;
    double median;
    double mean;
    int num_kmers;
    bool is_uniform_high;
};

// === HELPER FUNCTIONS ===

bool parse_line(const std::string& line, std::string& chrom, int& start, int& end, int& count) {
    std::stringstream ss(line);
    if (!(ss >> chrom >> start >> end >> count)) return false;
    return true;
}

double calculate_percentile_kmer(const std::vector<KmerData>& kmers, double percentile) {
    if (kmers.empty()) return 0.0;
    std::vector<int> counts;
    counts.reserve(kmers.size());
    for (const auto& kmer : kmers) counts.push_back(kmer.count);
    std::sort(counts.begin(), counts.end());
    
    double n = counts.size();
    double h = (n - 1) * (percentile / 100.0);
    int h_floor = static_cast<int>(std::floor(h));
    int h_ceil = static_cast<int>(std::ceil(h));
    
    if (h_floor < 0) return counts[0];
    if (h_ceil >= n) return counts[n-1];
    if (h_floor == h_ceil) return counts[h_floor];
    
    double weight = h - h_floor;
    return counts[h_floor] * (1.0 - weight) + counts[h_ceil] * weight;
}

double calculate_winsorized_mean(const std::vector<KmerData>& kmers, double lower_pct=5.0, double upper_pct=95.0) {
    if (kmers.empty()) return 0.0;
    double p_lower = calculate_percentile_kmer(kmers, lower_pct);
    double p_upper = calculate_percentile_kmer(kmers, upper_pct);
    
    double sum = 0;
    int count = 0;
    for (const auto& kmer : kmers) {
        if (kmer.count >= p_lower && kmer.count <= p_upper) {
            sum += kmer.count;
            count++;
        }
    }
    return (count > 0) ? (sum / count) : 0.0;
}

double calculate_iqr(const std::vector<KmerData>& kmers) {
    if (kmers.empty()) return 0.0;
    double q1 = calculate_percentile_kmer(kmers, 25.0);
    double q3 = calculate_percentile_kmer(kmers, 75.0);
    return q3 - q1;
}

WindowStats calculate_window_stats(const std::vector<KmerData>& kmers, int median) {
    WindowStats stats;
    stats.num_kmers = kmers.size();
    
    if (kmers.empty()) {
        stats.winsorized_mean = 0; stats.iqr = 0; stats.median = 0; stats.robust_cv = 0; stats.mean = 0;
        stats.is_uniform_high = false;
        return stats;
    }

    stats.median = calculate_percentile_kmer(kmers, 50.0);
    stats.winsorized_mean = calculate_winsorized_mean(kmers, 5.0, 95.0);
    stats.iqr = calculate_iqr(kmers);
    stats.robust_cv = (stats.median > 1e-6) ? (stats.iqr / stats.median) : 0.0;
    
    double sum = 0;
    for(const auto& k : kmers) sum += k.count;
    stats.mean = sum / kmers.size();

    bool standard_duplication = (stats.robust_cv < UNIFORM_CV_THRESHOLD) &&
                                (stats.winsorized_mean > median * DUPLICATION_MEAN_FACTOR);

    bool complex_amplification = (stats.winsorized_mean > median * DUPLICATION_MEAN_FACTOR) &&
                                 (stats.median > MIN_COUNT_THRESHOLD);
    
    bool is_biological_valid = (stats.median < median * MAX_BIOLOGICAL_CN_FACTOR);

    stats.is_uniform_high = (standard_duplication || complex_amplification) && is_biological_valid;
    return stats;
}

void export_window_data(std::ofstream& out, const std::string& type, const std::string& chrom, int w_idx, const std::vector<KmerData>& kmers) {
    out << type << "\t" << chrom << "\t" << w_idx * WINDOW_SIZE << "\t";
    for (size_t i = 0; i < kmers.size(); ++i) {
        out << kmers[i].count << (i == kmers.size() - 1 ? "" : ",");
    }
    out << "\n";
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_bed> <output_bed>" << std::endl;
        return 1;
    }

    std::string input_file = argv[1];
    std::string output_file = argv[2];

    std::cout << "=== Clean Preprocessing with Visualization Export ===" << std::endl;

    // 1. Calculate Median (Pass 1)
    std::ifstream infile1(input_file);
    if (!infile1.is_open()) {
        std::cerr << "Error opening input file!" << std::endl;
        return 1;
    }

    std::vector<long long> histogram(65536, 0);
    std::string line, chrom;
    int start, end, count;
    
    while (std::getline(infile1, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (parse_line(line, chrom, start, end, count)) {
            if (count >= 0 && count <= 65535) {
                histogram[count]++;
            }
        }
    }
    infile1.close();

    // Export Global Histogram
    std::ofstream hist_out("global_histogram.dat");
    for (int i = 0; i < 65536; ++i) {
        if (histogram[i] > 0) {
            hist_out << i << "\t" << histogram[i] << "\n";
        }
    }
    hist_out.close();
    std::cout << "Exported global_histogram.dat" << std::endl;

    const int MEDIAN_UPPER_BOUND = 5000;
    long long valid_kmers = 0;
    for (int i = 10; i < MEDIAN_UPPER_BOUND; ++i) valid_kmers += histogram[i];
    
    long long current_sum = 0;
    int median = 51;
    long long target = valid_kmers / 2;
    
    if (valid_kmers > 0) {
        for (int i = 10; i < MEDIAN_UPPER_BOUND; ++i) {
            current_sum += histogram[i];
            if (current_sum >= target) {
                median = i;
                break;
            }
        }
    }
    std::cout << "Calculated Median: " << median << std::endl;

    // 2. Window Processing (Pass 2)
    std::ifstream infile2(input_file);
    std::ofstream outfile(output_file);
    std::ofstream window_debug_out("window_samples.dat"); // For visualization
    
    std::map<std::string, std::map<int, std::vector<KmerData>>> chrom_window_kmers;

    while (std::getline(infile2, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (parse_line(line, chrom, start, end, count)) {
            if (count < MIN_COUNT_THRESHOLD) continue;

            int w_idx = start / WINDOW_SIZE;
            KmerData kmer;
            kmer.position = start;
            kmer.count = count;
            kmer.adjusted_count = (count >= SATURATION_THRESHOLD) ? (count * SATURATION_FACTOR) : count;
            
            chrom_window_kmers[chrom][w_idx].push_back(kmer);
        }
    }
    infile2.close();

    // Counters for visualization samples
    int saved_good = 0;
    int saved_bad = 0;
    int saved_random = 0;
    
    // Random generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // 3. Process Chromosomes
    for (auto& [current_chrom, window_kmers] : chrom_window_kmers) {
        std::vector<std::pair<int, double>> window_cns;

        for (auto& [w_idx, kmers] : window_kmers) {
            WindowStats stats = calculate_window_stats(kmers, median);
            
            // --- VISUALIZATION SAMPLING LOGIC ---
            bool save_this = false;
            std::string type = "";

            // Criteria for "Good" (Clean): Low CV, close to expected CNs (1, 2, etc.)
            if (saved_good < PLOTS_PER_CATEGORY && stats.robust_cv < 0.2 && stats.median > 20) {
                type = "Good_Clean";
                save_this = true;
                saved_good++;
            }
            // Criteria for "Bad" (Noisy): High CV
            else if (saved_bad < PLOTS_PER_CATEGORY && stats.robust_cv > 1.0 && stats.median > 20) {
                type = "Bad_Noisy";
                save_this = true;
                saved_bad++;
            }
            // Criteria for "Random"
            else if (saved_random < PLOTS_PER_CATEGORY && dis(gen) < 0.001) { // 0.1% chance to pick random
                type = "Random";
                save_this = true;
                saved_random++;
            }

            if (save_this) {
                export_window_data(window_debug_out, type, current_chrom, w_idx, kmers);
            }
            // ------------------------------------

            std::vector<double> kept_values;
            int threshold = stats.is_uniform_high ? EXTREME_REPEAT_ABSOLUTE : MAX_COUNT_THRESHOLD_BASE;

            for (const auto& kmer : kmers) {
                double val = kmer.adjusted_count;
                if (val > threshold) val = threshold; // Clamp
                kept_values.push_back(val);
            }

            double mean_count = 0;
            if (!kept_values.empty()) {
                mean_count = std::accumulate(kept_values.begin(), kept_values.end(), 0.0) / kept_values.size();
            }
            
            double cn = mean_count / median;
            window_cns.push_back({w_idx, cn});
        }

        // Smoothing & Output
        for (size_t i = 0; i < window_cns.size(); ++i) {
            int w_idx = window_cns[i].first;
            double smoothed_cn = window_cns[i].second;

            std::vector<double> neighbors;
            int radius = SMOOTH_WINDOW / 2;
            bool enough_neighbors = true;
            
            if (i < radius || i >= window_cns.size() - radius) {
                enough_neighbors = false;
            } else {
                for (int j = -radius; j <= radius; ++j) {
                     if (window_cns[i+j].first != w_idx + j) {
                         enough_neighbors = false;
                         break;
                     }
                     neighbors.push_back(window_cns[i+j].second);
                }
            }

            if (enough_neighbors) {
                smoothed_cn = 0;
                for (size_t k = 0; k < SMOOTH_WINDOW; ++k) {
                    smoothed_cn += neighbors[k] * SG_COEFFS[k];
                }
                if (smoothed_cn < 0) smoothed_cn = 0;
            }
            
            outfile << current_chrom << "\t"
                    << (w_idx * WINDOW_SIZE) << "\t"
                    << ((w_idx + 1) * WINDOW_SIZE) << "\t"
                    << smoothed_cn << "\n";
        }
    }
    
    outfile.close();
    window_debug_out.close();
    std::cout << "Done. Visualization data written to 'global_histogram.dat' and 'window_samples.dat'" << std::endl;
    return 0;
}
