#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// === CONFIGURATION STRUCT ===
struct Config {
  int window_size = 1000;
  int saturation_threshold = 65535;
  double saturation_factor = 1.7;
  double bio_threshold_factor = 150.0;
  int min_kmers_per_window = 5;
  int deletion_floor = 10;
  bool use_local_median = true;
  std::string input_file;
  std::string output_file;
};

// NOISE ELIMINATION PARAMETERS (used by apply_local_median_preprocessing)
const int LOCAL_MEDIAN_WINDOW = 11;
const double OUTLIER_IQR_MULTIPLIER = 3.0;
const int MIN_NEIGHBORS_REQUIRED = 5;

struct KmerData {
  int position;
  int count;
  double adjusted_count;
};

void print_usage(const char *program_name) {
  std::cerr << "Usage: " << program_name
            << " <input.bed> <output.bed> [options]" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Required:" << std::endl;
  std::cerr
      << "  <input.bed>           Input BED file (chrom, start, end, count)"
      << std::endl;
  std::cerr << "  <output.bed>          Output BED file" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Optional:" << std::endl;
  std::cerr << "  --window-size N       Window size in bp (default: 1000)"
            << std::endl;
  std::cerr
      << "  --bio-threshold F     Biological threshold factor (default: 50.0)"
      << std::endl;
  std::cerr << "  --no-local-median     Disable local median preprocessing"
            << std::endl;
  std::cerr << "  --help                Show this help message" << std::endl;
}

Config parse_arguments(int argc, char *argv[]) {
  Config config;

  if (argc < 3) {
    print_usage(argv[0]);
    exit(1);
  }

  // Check for --help first
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
      print_usage(argv[0]);
      exit(0);
    }
  }

  // First two positional arguments
  config.input_file = argv[1];
  config.output_file = argv[2];

  // Parse optional arguments
  for (int i = 3; i < argc; ++i) {
    if (strcmp(argv[i], "--window-size") == 0 && i + 1 < argc) {
      config.window_size = std::atoi(argv[++i]);
      if (config.window_size <= 0) {
        std::cerr << "Error: window-size must be positive" << std::endl;
        exit(1);
      }
    } else if (strcmp(argv[i], "--bio-threshold") == 0 && i + 1 < argc) {
      config.bio_threshold_factor = std::atof(argv[++i]);
      if (config.bio_threshold_factor <= 0) {
        std::cerr << "Error: bio-threshold must be positive" << std::endl;
        exit(1);
      }
    } else if (strcmp(argv[i], "--no-local-median") == 0) {
      config.use_local_median = false;
    } else {
      std::cerr << "Error: Unknown argument: " << argv[i] << std::endl;
      print_usage(argv[0]);
      exit(1);
    }
  }

  return config;
}

bool parse_line(const std::string &line, std::string &chrom, int &start,
                int &end, int &count) {
  std::stringstream ss(line);
  if (!(ss >> chrom >> start >> end >> count))
    return false;
  return true;
}

double calculate_median(const std::vector<double> &values) {
  if (values.empty())
    return 0.0;
  std::vector<double> sorted = values;
  std::sort(sorted.begin(), sorted.end());
  size_t n = sorted.size();
  if (n % 2 == 0) {
    return (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0;
  } else {
    return sorted[n / 2];
  }
}

double calculate_percentile(const std::vector<double> &sorted_values,
                            double percentile) {
  if (sorted_values.empty())
    return 0.0;
  double n = sorted_values.size();
  double h = (n - 1) * (percentile / 100.0);
  int h_floor = static_cast<int>(std::floor(h));
  int h_ceil = static_cast<int>(std::ceil(h));
  if (h_floor < 0)
    return sorted_values[0];
  if (h_ceil >= n)
    return sorted_values[n - 1];
  if (h_floor == h_ceil)
    return sorted_values[h_floor];
  double weight = h - h_floor;
  return sorted_values[h_floor] * (1.0 - weight) +
         sorted_values[h_ceil] * weight;
}

struct GaussianFit {
  double mu;
  double sigma;
  double amplitude;
  double noise_threshold;
};


GaussianFit fit_gaussian_to_main_peak(const std::vector<long long> &histogram) {
  GaussianFit result;
  int peak_pos = 20;
  long long peak_height = 0;
  for (int i = 20; i <= 150; ++i) {
    if (histogram[i] > peak_height) {
      peak_height = histogram[i];
      peak_pos = i;
    }
  }
  result.amplitude = peak_height;
  result.mu = peak_pos;

  long long half_max = peak_height / 2;
  int left_edge = peak_pos;
  for (int i = peak_pos; i >= 10; --i) {
    if (histogram[i] < half_max) {
      left_edge = i;
      break;
    }
  }
  int right_edge = peak_pos;
  for (int i = peak_pos; i <= 200; ++i) {
    if (histogram[i] < half_max) {
      right_edge = i;
      break;
    }
  }

  double fwhm = right_edge - left_edge;
  result.sigma = fwhm / 2.355;
  double poisson_sigma = std::sqrt(result.mu);
  if (result.sigma < poisson_sigma)
    result.sigma = poisson_sigma;

  result.noise_threshold = result.mu - 3.0 * result.sigma;
  if (result.noise_threshold < 10.0)
    result.noise_threshold = 10.0;
  return result;
}

double calculate_robust_percentile(const std::vector<long long> &histogram,
                                   int lower_bound, int upper_bound) {
  long long valid_kmers = 0;
  for (int i = lower_bound; i < upper_bound; ++i)
    valid_kmers += histogram[i];
  if (valid_kmers == 0)
    return 57.0;

  long long current_sum = 0;
  long long target = valid_kmers / 2;
  for (int i = lower_bound; i < upper_bound; ++i) {
    current_sum += histogram[i];
    if (current_sum >= target)
      return static_cast<double>(i);
  }
  return 57.0;
}


int apply_local_median_preprocessing(
    std::map<std::string, std::vector<KmerData>> &chrom_data) {
  int total_corrected = 0;
  for (auto &[chrom, kmers] : chrom_data) {
    if (kmers.empty())
      continue;
    std::sort(kmers.begin(), kmers.end(),
              [](const KmerData &a, const KmerData &b) {
                return a.position < b.position;
              });

    for (size_t i = 0; i < kmers.size(); ++i) {
      std::vector<double> local_values;
      int current_pos = kmers[i].position;

      for (int j = i - 1;
           j >= 0 && local_values.size() < LOCAL_MEDIAN_WINDOW / 2; --j) {
        if (kmers[j].position >= current_pos - 500) {
          local_values.push_back(kmers[j].adjusted_count);
        } else
          break;
      }
      for (size_t j = i + 1;
           j < kmers.size() && local_values.size() < LOCAL_MEDIAN_WINDOW - 1;
           ++j) {
        if (kmers[j].position <= current_pos + 500) {
          local_values.push_back(kmers[j].adjusted_count);
        } else
          break;
      }

      if (local_values.size() >= MIN_NEIGHBORS_REQUIRED) {
        std::sort(local_values.begin(), local_values.end());
        double local_median = calculate_median(local_values);
        double q1 = calculate_percentile(local_values, 25.0);
        double q3 = calculate_percentile(local_values, 75.0);
        double iqr = q3 - q1;
        double deviation = std::abs(kmers[i].adjusted_count - local_median);
        if (deviation > OUTLIER_IQR_MULTIPLIER * iqr) {
          kmers[i].adjusted_count = local_median;
          total_corrected++;
        }
      }
    }
  }
  return total_corrected;
}

double calculate_biological_winsorized_mean(const std::vector<double> &counts,
                                            double biological_threshold,
                                            int &num_filtered) {
  if (counts.empty())
    return 0.0;

  std::vector<double> valid_counts;
  num_filtered = 0;

  // Step 1: Biological filtering
  for (double count : counts) {
    if (count <= biological_threshold) {
      valid_counts.push_back(count);
    } else {
      num_filtered++;
    }
  }

  if (valid_counts.empty())
    return 0.0;

  // Step 2: Winsorized Mean (5-95 percentile) - BASELINE
  std::vector<double> sorted = valid_counts;
  std::sort(sorted.begin(), sorted.end());

  size_t n = sorted.size();
  size_t lower_idx = static_cast<size_t>(n * 0.05);
  size_t upper_idx = static_cast<size_t>(n * 0.95);

  if (upper_idx >= n)
    upper_idx = n - 1;
  if (lower_idx > upper_idx)
    return calculate_median(sorted);

  double sum = 0.0;
  int count = 0;
  for (size_t i = lower_idx; i <= upper_idx; ++i) {
    sum += sorted[i];
    count++;
  }

  return (count > 0) ? (sum / count) : 0.0;
}

int main(int argc, char *argv[]) {
  Config config = parse_arguments(argc, argv);

  std::cout << "========================================" << std::endl;
  std::cout << "KonuSeg v5 - Bio Winsorized (Parametric)" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Input: " << config.input_file << std::endl;
  std::cout << "Output: " << config.output_file << std::endl;
  std::cout << "Window Size: " << config.window_size << " bp" << std::endl;
  std::cout << "Bio Threshold Factor: " << config.bio_threshold_factor << "x"
            << std::endl;
  std::cout << "Local Median: "
            << (config.use_local_median ? "enabled" : "disabled") << std::endl;
  std::cout << std::endl;

  // PASS 1: BUILD HISTOGRAM
  std::cout << "Step 1a: Building genome-wide histogram..." << std::endl;
  std::ifstream infile(config.input_file);
  if (!infile.is_open()) {
    std::cerr << "Error: Cannot open input file" << std::endl;
    return 1;
  }

  std::map<std::string, std::vector<KmerData>> chrom_data;
  std::vector<long long> histogram(65536, 0);
  std::string line, chrom;
  int start, end, count;
  long long total_kmers = 0;

  while (std::getline(infile, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    if (!parse_line(line, chrom, start, end, count))
      continue;
    total_kmers++;
    if (count >= 0 && count <= 65535)
      histogram[count]++;

    double adjusted = count;
    if (count >= config.saturation_threshold)
      adjusted = count * config.saturation_factor;

    KmerData kmer;
    kmer.position = start;
    kmer.count = count;
    kmer.adjusted_count = adjusted;
    chrom_data[chrom].push_back(kmer);
  }
  infile.close();

  std::cout << "  Read " << chrom_data.size() << " chromosomes" << std::endl;
  std::cout << "  Total k-mers: " << total_kmers << std::endl << std::endl;

  // GAUSSIAN FIT MEDIAN
  std::cout << "Step 1b: Calculating genome-wide median..." << std::endl;
  GaussianFit gaussian_fit = fit_gaussian_to_main_peak(histogram);
  double gaussian_median = gaussian_fit.mu;

  const int MEDIAN_UPPER_BOUND = 5000;
  double robust_median =
      calculate_robust_percentile(histogram, 10, MEDIAN_UPPER_BOUND);

  std::cout << "  Gaussian fit (peak detection): mu=" << gaussian_median
            << ", sigma=" << gaussian_fit.sigma << std::endl;
  std::cout << "  Robust percentile (reference): " << robust_median
            << std::endl;

  double GENOME_WIDE_MEDIAN = gaussian_median;
  std::cout << "  Using Gaussian fit: " << GENOME_WIDE_MEDIAN << std::endl;
  std::cout << std::endl;

  // BIOLOGICAL THRESHOLD
  double BIOLOGICAL_THRESHOLD =
      GENOME_WIDE_MEDIAN * config.bio_threshold_factor;
  std::cout << "Step 1c: Biological filtering parameters..." << std::endl;
  std::cout << "  Max Biological CN: " << config.bio_threshold_factor
            << std::endl;
  std::cout << "  Biological Threshold: " << BIOLOGICAL_THRESHOLD << std::endl;
  std::cout << std::endl;

  // NOISE ELIMINATION
  std::cout << "Step 1d: Noise elimination..." << std::endl;

  long long kmers_before = total_kmers;
  for (auto &[chrom_name, kmers] : chrom_data) {
    auto it =
        std::remove_if(kmers.begin(), kmers.end(),
                       [&gaussian_fit, &config](const KmerData &k) {
                         if (k.adjusted_count <= config.deletion_floor)
                           return false;
                         return k.adjusted_count < gaussian_fit.noise_threshold;
                       });
    kmers.erase(it, kmers.end());
  }

  long long kmers_after = 0;
  for (const auto &[chrom_name, kmers] : chrom_data)
    kmers_after += kmers.size();

  double removal_pct =
      (kmers_before > 0) ? (100.0 * (kmers_before - kmers_after) / kmers_before)
                         : 0.0;
  std::cout << "  Threshold filtering: " << kmers_before << " -> "
            << kmers_after << " (" << removal_pct << "% removed)" << std::endl;

  // LOCAL MEDIAN PREPROCESSING (optional)
  int outliers_corrected = 0;
  if (config.use_local_median) {
    outliers_corrected = apply_local_median_preprocessing(chrom_data);
    std::cout << "  Local median correction: " << outliers_corrected
              << " k-mers" << std::endl;
  } else {
    std::cout << "  Local median correction: SKIPPED (--no-local-median)"
              << std::endl;
  }

  std::cout << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "GENOME-WIDE BASELINE: " << GENOME_WIDE_MEDIAN << std::endl;
  std::cout << "WINDOW SIZE: " << config.window_size << " bp" << std::endl;
  std::cout << "AGGREGATION METHOD: Bio Filter (" << config.bio_threshold_factor
            << "x) + Winsorized Mean (10-90%)" << std::endl;
  std::cout << "========================================" << std::endl
            << std::endl;

  // PASS 2: WINDOW AGGREGATION
  std::cout << "Step 2: Window aggregation..." << std::endl;
  std::ofstream outfile(config.output_file);
  if (!outfile.is_open()) {
    std::cerr << "Error: Cannot open output file" << std::endl;
    return 1;
  }

  outfile << "#chrom\tstart\tend\tCN\tmean_count\tlog_ratio\tnum_kmers\tnum_"
             "filtered_repeats"
          << std::endl;

  long long total_windows = 0;
  long long total_filtered_repeats = 0;

  for (auto &[chr, kmers] : chrom_data) {
    if (kmers.empty())
      continue;

    std::sort(kmers.begin(), kmers.end(),
              [](const KmerData &a, const KmerData &b) {
                return a.position < b.position;
              });

    std::cout << "  Processing " << chr << " (" << kmers.size() << " k-mers)..."
              << std::endl;

    for (size_t i = 0; i < kmers.size();) {
      int window_start =
          (kmers[i].position / config.window_size) * config.window_size;
      int window_end = window_start + config.window_size;

      std::vector<double> counts;
      while (i < kmers.size() && kmers[i].position < window_end) {
        counts.push_back(kmers[i].adjusted_count);
        i++;
      }

      if (counts.size() >= static_cast<size_t>(config.min_kmers_per_window)) {
        int num_filtered = 0;
        double mean_count = calculate_biological_winsorized_mean(
            counts, BIOLOGICAL_THRESHOLD, num_filtered);

        total_filtered_repeats += num_filtered;

        int num_valid = counts.size() - num_filtered;
        if (num_valid >= config.min_kmers_per_window && mean_count > 0) {
          double cn = mean_count / GENOME_WIDE_MEDIAN;
          double log_ratio = std::log2(std::max(cn, 0.01));

          outfile << chr << "\t" << window_start << "\t" << window_end << "\t"
                  << cn << "\t" << mean_count << "\t" << log_ratio << "\t"
                  << num_valid << "\t" << num_filtered << std::endl;

          total_windows++;
        }
      }
    }
  }

  outfile.close();

  std::cout << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Output: " << config.output_file << std::endl;
  std::cout << "Total windows: " << total_windows << std::endl;
  std::cout << "Total filtered repeats: " << total_filtered_repeats
            << std::endl;
  std::cout << "CN values are CONTINUOUS (floating-point)" << std::endl;
  std::cout << "Window Size: " << config.window_size << " bp" << std::endl;
  std::cout << "Aggregation: Bio Filter (" << config.bio_threshold_factor
            << "x) + Winsorized Mean (10-90%)" << std::endl;
  std::cout << "========================================" << std::endl;

  return 0;
}
