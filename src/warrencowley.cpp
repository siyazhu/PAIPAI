#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "analysis_utils.h"

namespace fs = std::filesystem;
using namespace paipai_analysis;

static void usage()
{
    std::cerr << "Usage: warrencowley STATE_DIR "
                 "[--inter-metal-cutoff R] [--metal-metal-cutoff R] "
                 "[--output struc_WC]\n";
}

static std::vector<double> metal_composition(const SaveData& save)
{
    std::vector<double> global(save.metal_species.size(), 0.0);
    for (int t : save.metal_types) {
        if (t >= 0 && t < (int)global.size()) global[t] += 1.0;
    }
    for (double& x : global) {
        if (!save.metal_types.empty()) x /= (double)save.metal_types.size();
    }
    return global;
}

static std::string format_wc(double value)
{
    std::ostringstream ss;
    ss.precision(12);
    ss << value;
    return ss.str();
}

static std::vector<std::pair<std::string, std::string>>
compute_inter_metal_wc(const SaveData& save,
                       const PoscarData& contcar,
                       const std::vector<Vec3>& metal_coords,
                       const std::array<Vec3, 3>& inv,
                       double cutoff)
{
    std::vector<double> global = metal_composition(save);
    std::vector<std::pair<std::string, std::string>> rows;

    for (int it = 0; it < (int)save.interstitial_species.size(); ++it) {
        std::vector<int> sites;
        for (int s = 0; s < (int)save.interstitial_site_types.size(); ++s) {
            if (save.interstitial_site_types[s] == it) sites.push_back(s);
        }

        std::vector<int> neighbor_counts(save.metal_species.size(), 0);
        int total_neighbors = 0;
        for (int s : sites) {
            for (int m = 0; m < (int)metal_coords.size(); ++m) {
                double d = minimum_image_distance(save.interstitial_pos[s], metal_coords[m], contcar.cell, inv);
                if (d <= cutoff) {
                    int mt = save.metal_types[m];
                    if (mt >= 0 && mt < (int)neighbor_counts.size()) {
                        neighbor_counts[mt]++;
                        total_neighbors++;
                    }
                }
            }
        }

        for (int mt = 0; mt < (int)save.metal_species.size(); ++mt) {
            std::string label = "WC_" + save.interstitial_species[it] + "_" + save.metal_species[mt];
            std::string value = "nan";
            if (!sites.empty() && total_neighbors > 0 && global[mt] > 0.0) {
                double p_local = (double)neighbor_counts[mt] / (double)total_neighbors;
                value = format_wc(1.0 - p_local / global[mt]);
            }
            rows.push_back({label, value});
        }
    }
    return rows;
}

static std::vector<std::pair<std::string, std::string>>
compute_metal_metal_wc(const SaveData& save,
                       const PoscarData& contcar,
                       const std::vector<Vec3>& metal_coords,
                       const std::array<Vec3, 3>& inv,
                       double cutoff)
{
    std::vector<double> global = metal_composition(save);
    std::vector<std::pair<std::string, std::string>> rows;

    for (int center_type = 0; center_type < (int)save.metal_species.size(); ++center_type) {
        std::vector<int> neighbor_counts(save.metal_species.size(), 0);
        int total_neighbors = 0;
        int center_count = 0;

        for (int i = 0; i < (int)metal_coords.size(); ++i) {
            if (save.metal_types[i] != center_type) continue;
            center_count++;
            for (int j = 0; j < (int)metal_coords.size(); ++j) {
                if (i == j) continue;
                double d = minimum_image_distance(metal_coords[i], metal_coords[j], contcar.cell, inv);
                if (d <= cutoff) {
                    int nt = save.metal_types[j];
                    if (nt >= 0 && nt < (int)neighbor_counts.size()) {
                        neighbor_counts[nt]++;
                        total_neighbors++;
                    }
                }
            }
        }

        for (int neigh_type = 0; neigh_type < (int)save.metal_species.size(); ++neigh_type) {
            std::string label = "WC_" + save.metal_species[center_type] + "_" + save.metal_species[neigh_type];
            std::string value = "nan";
            if (center_count > 0 && total_neighbors > 0 && global[neigh_type] > 0.0) {
                double p_local = (double)neighbor_counts[neigh_type] / (double)total_neighbors;
                value = format_wc(1.0 - p_local / global[neigh_type]);
            }
            rows.push_back({label, value});
        }
    }
    return rows;
}

static std::vector<std::pair<std::string, std::string>>
compute_wc(const fs::path& state_dir, double inter_metal_cutoff, double metal_metal_cutoff)
{
    SaveData save = read_save(state_dir / "SAVE");
    PoscarData contcar = read_poscar(state_dir / "CONTCAR");
    auto metal_coords = metal_coords_by_save_order(save, contcar);
    auto inv = inverse_cell(contcar.cell);

    std::vector<std::pair<std::string, std::string>> rows;
    if (inter_metal_cutoff > 0.0) {
        auto part = compute_inter_metal_wc(save, contcar, metal_coords, inv, inter_metal_cutoff);
        rows.insert(rows.end(), part.begin(), part.end());
    }
    if (metal_metal_cutoff > 0.0) {
        auto part = compute_metal_metal_wc(save, contcar, metal_coords, inv, metal_metal_cutoff);
        rows.insert(rows.end(), part.begin(), part.end());
    }
    return rows;
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        usage();
        return 2;
    }

    fs::path state_dir = argv[1];
    double inter_metal_cutoff = -1.0;
    double metal_metal_cutoff = -1.0;
    fs::path output = "struc_WC";

    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--inter-metal-cutoff" && i + 1 < argc) {
            inter_metal_cutoff = std::stod(argv[++i]);
        } else if (a == "--metal-metal-cutoff" && i + 1 < argc) {
            metal_metal_cutoff = std::stod(argv[++i]);
        } else if (a == "--output" && i + 1 < argc) {
            output = argv[++i];
        } else if (a == "-h" || a == "--help") {
            usage();
            return 0;
        } else {
            std::cerr << "Unknown argument: " << a << "\n";
            usage();
            return 2;
        }
    }

    if (inter_metal_cutoff <= 0.0 && metal_metal_cutoff <= 0.0) {
        std::cerr << "At least one of --inter-metal-cutoff or --metal-metal-cutoff must be positive.\n";
        return 2;
    }
    if (!output.is_absolute()) output = state_dir / output;

    try {
        auto rows = compute_wc(state_dir, inter_metal_cutoff, metal_metal_cutoff);
        std::ofstream out(output);
        if (!out) throw std::runtime_error("cannot write " + output.string());
        for (const auto& row : rows) {
            out << row.first << "\t" << row.second << "\n";
        }
        std::cout << "Wrote " << output << "\n";
    } catch (const std::exception& e) {
        std::cerr << "warrencowley error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
