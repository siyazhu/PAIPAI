#ifndef PAIPAI_PREFAST_H
#define PAIPAI_PREFAST_H

#include <filesystem>
#include <iosfwd>
#include <string>
#include <unordered_map>
#include <vector>

#include "structure.h"

namespace paipai_prefast {

enum class SiteFamily {
    MM = 0,
    MI = 1,
    II = 2,
};

struct PrefastConfig {
    bool enabled = false;
    int candidates_per_slot = 6;
    std::string basis = "ref-dz";
    int nshells = 3;
    double peak_scan_cutoff = 6.0;
    double peak_tol = 0.12;
    double sigma_small = 0.10;
    double sigma_large = 0.30;
    double cutoff_margin = 0.50;
    double learning_rate = 0.01;
    double lms_epsilon = 1.0e-12;
    double weight_decay = 0.0;
};

struct SparseDescriptor {
    std::unordered_map<int, double> values;

    double norm2() const;
};

struct TrialScore {
    SparseDescriptor delta;
    double dE_pred = 0.0;
    double norm_dD = 0.0;
};

class PrefastModel {
public:
    bool build(const Structure& reference, const PrefastConfig& cfg);
    bool enabled() const { return built_ && config_.enabled; }

    TrialScore score(const Structure& before, const Structure& after);
    void update(const SparseDescriptor& delta, double dE_true);
    void append_learning_log(const std::filesystem::path& path,
                             int step,
                             const std::string& trial_id,
                             double dE_pred,
                             double dE_true,
                             bool accepted) const;

    void print_startup_log(std::ostream& out) const;
    void write_basis_log(const std::filesystem::path& path) const;

private:
    struct BasisValue {
        int shell = 0;
        int zeta = 0;
        double value = 0.0;
    };

    struct RefPair {
        int i = -1;
        int j = -1;
        double r0 = 0.0;
        SiteFamily family = SiteFamily::MM;
        std::vector<BasisValue> basis_values;
    };

    PrefastConfig config_;
    bool built_ = false;
    int n_metal_ = 0;
    int n_interstitial_ = 0;
    int n_sites_ = 0;
    std::vector<int> real_atomic_numbers_;
    std::vector<RefPair> pairs_;
    std::vector<std::vector<int>> adjacency_;
    std::vector<std::vector<double>> shell_centers_;
    std::vector<std::vector<int>> shell_counts_;
    std::vector<double> family_cutoffs_;
    std::unordered_map<std::string, int> feature_index_;
    std::vector<std::string> feature_names_;
    std::vector<double> weights_;

    int site_atomic_number(const Structure& s, int site_id) const;
    int feature_index(SiteFamily family, int z1, int z2, int shell, int zeta);
    double predict(const SparseDescriptor& delta) const;
};

} // namespace paipai_prefast

#endif
