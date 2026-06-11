#include "prefast.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>

#include "constants.h"

namespace fs = std::filesystem;

namespace paipai_prefast {

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

std::string family_name(SiteFamily f)
{
    if (f == SiteFamily::MM) return "MM";
    if (f == SiteFamily::MI) return "MI";
    return "II";
}

std::vector<Real> vec3(Real x, Real y, Real z)
{
    return {x, y, z};
}

double dot3(const std::vector<Real>& a, const std::vector<Real>& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

std::vector<Real> matvec_cell(Real f1, Real f2, Real f3,
                              const std::vector<Real>& a,
                              const std::vector<Real>& b,
                              const std::vector<Real>& c)
{
    return vec3(f1 * a[0] + f2 * b[0] + f3 * c[0],
                f1 * a[1] + f2 * b[1] + f3 * c[1],
                f1 * a[2] + f2 * b[2] + f3 * c[2]);
}

bool inverse3x3_columns(const std::vector<Real>& a,
                        const std::vector<Real>& b,
                        const std::vector<Real>& c,
                        Real inv[3][3])
{
    Real m00 = a[0], m01 = b[0], m02 = c[0];
    Real m10 = a[1], m11 = b[1], m12 = c[1];
    Real m20 = a[2], m21 = b[2], m22 = c[2];
    Real det = m00 * (m11 * m22 - m12 * m21)
             - m01 * (m10 * m22 - m12 * m20)
             + m02 * (m10 * m21 - m11 * m20);
    if (std::fabs(det) < 1e-14) return false;
    Real id = 1.0 / det;
    inv[0][0] =  (m11 * m22 - m12 * m21) * id;
    inv[0][1] = -(m01 * m22 - m02 * m21) * id;
    inv[0][2] =  (m01 * m12 - m02 * m11) * id;
    inv[1][0] = -(m10 * m22 - m12 * m20) * id;
    inv[1][1] =  (m00 * m22 - m02 * m20) * id;
    inv[1][2] = -(m00 * m12 - m02 * m10) * id;
    inv[2][0] =  (m10 * m21 - m11 * m20) * id;
    inv[2][1] = -(m00 * m21 - m01 * m20) * id;
    inv[2][2] =  (m00 * m11 - m01 * m10) * id;
    return true;
}

std::vector<Real> cart_to_frac(const std::vector<Real>& r, Real inv[3][3])
{
    return vec3(inv[0][0] * r[0] + inv[0][1] * r[1] + inv[0][2] * r[2],
                inv[1][0] * r[0] + inv[1][1] * r[1] + inv[1][2] * r[2],
                inv[2][0] * r[0] + inv[2][1] * r[1] + inv[2][2] * r[2]);
}

std::vector<Real> minimum_image_delta(const std::vector<Real>& delta,
                                      const std::vector<Real>& a,
                                      const std::vector<Real>& b,
                                      const std::vector<Real>& c,
                                      Real inv[3][3])
{
    std::vector<Real> f = cart_to_frac(delta, inv);
    f[0] -= std::round(f[0]);
    f[1] -= std::round(f[1]);
    f[2] -= std::round(f[2]);
    return matvec_cell(f[0], f[1], f[2], a, b, c);
}

std::vector<Real> site_position(const Structure& s, int site_id)
{
    if (site_id < s.num_metallic_atoms) return s.pos[site_id];
    return s.interstitial_pos[site_id - s.num_metallic_atoms];
}

SiteFamily pair_family(int i, int j, int n_metal)
{
    bool mi = i < n_metal;
    bool mj = j < n_metal;
    if (mi && mj) return SiteFamily::MM;
    if (mi || mj) return SiteFamily::MI;
    return SiteFamily::II;
}

std::vector<std::pair<double, int>> cluster_shells(std::vector<double> distances,
                                                   double tol,
                                                   int nshells)
{
    std::vector<std::pair<double, int>> shells;
    if (distances.empty()) return shells;
    std::sort(distances.begin(), distances.end());

    double sum = 0.0;
    int count = 0;
    double center = distances.front();
    for (double r : distances) {
        if (count > 0 && std::fabs(r - center) > tol) {
            shells.push_back({sum / count, count});
            if ((int)shells.size() >= nshells) return shells;
            sum = 0.0;
            count = 0;
            center = r;
        }
        sum += r;
        ++count;
        center = sum / count;
    }
    if (count > 0 && (int)shells.size() < nshells) {
        shells.push_back({sum / count, count});
    }
    return shells;
}

double cosine_cutoff(double r, double cutoff)
{
    if (cutoff <= 0.0 || r >= cutoff) return 0.0;
    return 0.5 * (std::cos(kPi * r / cutoff) + 1.0);
}

} // namespace

double SparseDescriptor::norm2() const
{
    double n = 0.0;
    for (const auto& kv : values) n += kv.second * kv.second;
    return n;
}

bool PrefastModel::build(const Structure& reference, const PrefastConfig& cfg)
{
    config_ = cfg;
    built_ = false;
    pairs_.clear();
    adjacency_.clear();
    shell_centers_.assign(3, {});
    shell_counts_.assign(3, {});
    family_cutoffs_.assign(3, 0.0);
    feature_index_.clear();
    feature_names_.clear();
    weights_.clear();

    if (!config_.enabled) return true;
    if (config_.basis != "ref-dz") {
        std::cerr << "[prefast] unsupported basis '" << config_.basis
                  << "'; prefast disabled.\n";
        config_.enabled = false;
        return false;
    }

    n_metal_ = reference.num_metallic_atoms;
    n_interstitial_ = reference.num_interstitial;
    n_sites_ = n_metal_ + n_interstitial_;
    adjacency_.assign(n_sites_, {});

    std::set<int> zset;
    for (int z : reference.type) zset.insert(z);
    for (int z : reference.interstitial_type) zset.insert(z);
    real_atomic_numbers_.assign(zset.begin(), zset.end());

    std::vector<Real> a = vec3(reference.cell_x1, reference.cell_y1, reference.cell_z1);
    std::vector<Real> b = vec3(reference.cell_x2, reference.cell_y2, reference.cell_z2);
    std::vector<Real> c = vec3(reference.cell_x3, reference.cell_y3, reference.cell_z3);
    Real inv[3][3];
    bool has_inv = inverse3x3_columns(a, b, c, inv);
    if (!has_inv) {
        std::cerr << "[prefast] warning: singular reference cell; using direct Cartesian distances.\n";
    }

    std::vector<std::vector<double>> distances_by_family(3);
    for (int i = 0; i < n_sites_; ++i) {
        for (int j = i + 1; j < n_sites_; ++j) {
            std::vector<Real> ri = site_position(reference, i);
            std::vector<Real> rj = site_position(reference, j);
            std::vector<Real> d = vec3(rj[0] - ri[0], rj[1] - ri[1], rj[2] - ri[2]);
            if (has_inv) d = minimum_image_delta(d, a, b, c, inv);
            double r0 = std::sqrt(dot3(d, d));
            if (r0 >= config_.peak_scan_cutoff) continue;

            SiteFamily fam = pair_family(i, j, n_metal_);
            RefPair p;
            p.i = i;
            p.j = j;
            p.r0 = r0;
            p.family = fam;
            int fi = static_cast<int>(fam);
            distances_by_family[fi].push_back(r0);
            pairs_.push_back(p);
        }
    }

    for (int fi = 0; fi < 3; ++fi) {
        auto shells = cluster_shells(distances_by_family[fi],
                                     config_.peak_tol,
                                     config_.nshells);
        for (const auto& sh : shells) {
            shell_centers_[fi].push_back(sh.first);
            shell_counts_[fi].push_back(sh.second);
        }
        if (shell_centers_[fi].empty()) {
            std::cerr << "[prefast] warning: no " << family_name((SiteFamily)fi)
                      << " distance shell found; that family contributes no descriptor terms.\n";
            family_cutoffs_[fi] = 0.0;
        } else {
            family_cutoffs_[fi] = shell_centers_[fi].back() + config_.cutoff_margin;
            if ((int)shell_centers_[fi].size() < config_.nshells) {
                std::cerr << "[prefast] warning: only " << shell_centers_[fi].size()
                          << " " << family_name((SiteFamily)fi)
                          << " shells found; requested " << config_.nshells << ".\n";
            }
        }
    }

    for (int pi = 0; pi < (int)pairs_.size(); ++pi) {
        RefPair& p = pairs_[pi];
        int fi = static_cast<int>(p.family);
        double fc = cosine_cutoff(p.r0, family_cutoffs_[fi]);
        for (int k = 0; k < (int)shell_centers_[fi].size(); ++k) {
            double mu = shell_centers_[fi][k];
            double sigmas[2] = {config_.sigma_small, config_.sigma_large};
            for (int zeta = 0; zeta < 2; ++zeta) {
                double sigma = sigmas[zeta];
                if (sigma <= 0.0) continue;
                double x = (p.r0 - mu) / sigma;
                double v = std::exp(-0.5 * x * x) * fc;
                if (v != 0.0) p.basis_values.push_back({k, zeta, v});
            }
        }
        adjacency_[p.i].push_back(pi);
        adjacency_[p.j].push_back(pi);
    }

    for (int fi = 0; fi < 3; ++fi) {
        SiteFamily fam = (SiteFamily)fi;
        for (size_t ai = 0; ai < real_atomic_numbers_.size(); ++ai) {
            for (size_t bi = ai; bi < real_atomic_numbers_.size(); ++bi) {
                for (int k = 0; k < (int)shell_centers_[fi].size(); ++k) {
                    feature_index(fam, real_atomic_numbers_[ai], real_atomic_numbers_[bi], k, 0);
                    feature_index(fam, real_atomic_numbers_[ai], real_atomic_numbers_[bi], k, 1);
                }
            }
        }
    }

    built_ = true;
    return true;
}

int PrefastModel::site_atomic_number(const Structure& s, int site_id) const
{
    if (site_id < s.num_metallic_atoms) {
        int t = s.atomtype[site_id];
        if (t < 0 || t >= (int)s.type.size()) return 0;
        return s.type[t];
    }
    int is = site_id - s.num_metallic_atoms;
    int t = s.interstitial_postype[is];
    if (t < 0) return 0;
    if (t >= (int)s.interstitial_type.size()) return 0;
    return s.interstitial_type[t];
}

int PrefastModel::feature_index(SiteFamily family, int z1, int z2, int shell, int zeta)
{
    if (z2 < z1) std::swap(z1, z2);
    std::ostringstream key;
    key << family_name(family) << ":" << z1 << "-" << z2
        << ":s" << shell << ":z" << zeta;
    std::string k = key.str();
    auto it = feature_index_.find(k);
    if (it != feature_index_.end()) return it->second;

    int idx = (int)feature_names_.size();
    feature_index_[k] = idx;
    feature_names_.push_back(k);
    weights_.push_back(0.0);
    return idx;
}

double PrefastModel::predict(const SparseDescriptor& delta) const
{
    double y = 0.0;
    for (const auto& kv : delta.values) {
        if (kv.first >= 0 && kv.first < (int)weights_.size()) {
            y += weights_[kv.first] * kv.second;
        }
    }
    return y;
}

TrialScore PrefastModel::score(const Structure& before, const Structure& after)
{
    TrialScore ts;
    if (!enabled()) return ts;

    std::vector<char> changed(n_sites_, 0);
    for (int i = 0; i < n_sites_; ++i) {
        if (site_atomic_number(before, i) != site_atomic_number(after, i)) {
            changed[i] = 1;
        }
    }

    std::vector<char> affected_pair(pairs_.size(), 0);
    for (int i = 0; i < n_sites_; ++i) {
        if (!changed[i]) continue;
        for (int pi : adjacency_[i]) affected_pair[pi] = 1;
    }

    for (int pi = 0; pi < (int)pairs_.size(); ++pi) {
        if (!affected_pair[pi]) continue;
        const RefPair& p = pairs_[pi];
        int old_a = site_atomic_number(before, p.i);
        int old_b = site_atomic_number(before, p.j);
        int new_a = site_atomic_number(after, p.i);
        int new_b = site_atomic_number(after, p.j);

        if (old_a > 0 && old_b > 0) {
            for (const auto& bv : p.basis_values) {
                int idx = feature_index(p.family, old_a, old_b, bv.shell, bv.zeta);
                ts.delta.values[idx] -= bv.value;
            }
        }
        if (new_a > 0 && new_b > 0) {
            for (const auto& bv : p.basis_values) {
                int idx = feature_index(p.family, new_a, new_b, bv.shell, bv.zeta);
                ts.delta.values[idx] += bv.value;
            }
        }
    }

    for (auto it = ts.delta.values.begin(); it != ts.delta.values.end(); ) {
        if (std::fabs(it->second) < 1e-15) it = ts.delta.values.erase(it);
        else ++it;
    }

    ts.norm_dD = std::sqrt(ts.delta.norm2());
    ts.dE_pred = predict(ts.delta);
    return ts;
}

void PrefastModel::update(const SparseDescriptor& delta, double dE_true)
{
    if (!enabled()) return;
    double pred = predict(delta);
    double err = dE_true - pred;
    double norm2 = delta.norm2();

    if (config_.weight_decay > 0.0) {
        double factor = std::max(0.0, 1.0 - config_.weight_decay);
        for (double& w : weights_) w *= factor;
    }

    double scale = config_.learning_rate * err / (config_.lms_epsilon + norm2);
    for (const auto& kv : delta.values) {
        if (kv.first >= 0 && kv.first < (int)weights_.size()) {
            weights_[kv.first] += scale * kv.second;
        }
    }
}

void PrefastModel::append_learning_log(const fs::path& path,
                                       int step,
                                       const std::string& trial_id,
                                       double dE_pred,
                                       double dE_true,
                                       bool accepted) const
{
    if (!enabled()) return;
    bool need_header = !fs::exists(path);
    std::ofstream out(path, std::ios::app);
    if (!out) return;
    if (need_header) {
        out << "step\ttrial_id\tdE_pred\tdE_true\terror\tlearning_rate\taccepted\n";
    }
    out << step << "\t" << trial_id
        << "\t" << std::setprecision(12) << dE_pred
        << "\t" << dE_true
        << "\t" << (dE_true - dE_pred)
        << "\t" << config_.learning_rate
        << "\t" << (accepted ? "accepted" : "rejected")
        << "\n";
}

void PrefastModel::print_startup_log(std::ostream& out) const
{
    if (!enabled()) return;
    out << "[prefast] enabled\n";
    out << "[prefast] descriptor = reference-lattice double-zeta radial basis\n";
    out << "[prefast] geometry = reference\n";
    out << "[prefast] peak scan cutoff = " << std::fixed << std::setprecision(3)
        << config_.peak_scan_cutoff << " Angstrom\n";
    out << "[prefast] peak tolerance = " << config_.peak_tol << " Angstrom\n";
    out << "[prefast] nshells = " << config_.nshells << "\n";
    for (int fi = 0; fi < 3; ++fi) {
        out << "[prefast] " << family_name((SiteFamily)fi) << " shell centers:\n";
        for (int k = 0; k < (int)shell_centers_[fi].size(); ++k) {
            out << "    shell " << (k + 1)
                << ": mu = " << shell_centers_[fi][k]
                << " count = " << shell_counts_[fi][k] << "\n";
        }
    }
    out << "[prefast] zeta widths: small = " << config_.sigma_small
        << " Angstrom, large = " << config_.sigma_large << " Angstrom\n";
    out << "[prefast] cutoff function = cosine\n";
    out << "[prefast] descriptor channels = " << feature_names_.size() << "\n";
    out << "[prefast] learning rule = normalized LMS\n";
    out << "[prefast] learning rate = " << config_.learning_rate << "\n";
    out << "[prefast] NOTE: This descriptor is used only for trial ranking/adaptive proposal ordering.\n";
    out << "[prefast] Exact MLIP-relaxed energies are still used for MC acceptance.\n";
}

void PrefastModel::write_basis_log(const fs::path& path) const
{
    if (!enabled()) return;
    std::ofstream out(path);
    if (!out) return;
    out << "family\tshell_index\tmu\tsigma_small\tsigma_large\tcutoff\tnumber_of_pairs_in_peak\n";
    for (int fi = 0; fi < 3; ++fi) {
        for (int k = 0; k < (int)shell_centers_[fi].size(); ++k) {
            out << family_name((SiteFamily)fi)
                << "\t" << (k + 1)
                << "\t" << std::setprecision(12) << shell_centers_[fi][k]
                << "\t" << config_.sigma_small
                << "\t" << config_.sigma_large
                << "\t" << family_cutoffs_[fi]
                << "\t" << shell_counts_[fi][k]
                << "\n";
        }
    }
}

} // namespace paipai_prefast
