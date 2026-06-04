#ifndef PAIPAI_ANALYSIS_UTILS_H
#define PAIPAI_ANALYSIS_UTILS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "json.hpp"

namespace paipai_analysis {

namespace fs = std::filesystem;
using json = nlohmann::json;
using Vec3 = std::array<double, 3>;

struct SaveData {
    std::vector<std::string> metal_species;
    std::vector<int> metal_counts;
    std::vector<std::string> interstitial_species;
    std::vector<int> interstitial_counts;
    int num_interstitial_sites = 0;
    std::vector<int> metal_types;
    std::vector<Vec3> interstitial_pos;
    std::vector<int> interstitial_site_types;
    std::array<Vec3, 3> cell{};
};

struct PoscarData {
    std::vector<std::string> species;
    std::vector<int> counts;
    std::vector<Vec3> coords;
    std::array<Vec3, 3> cell{};
};

inline std::string trim(const std::string& s)
{
    size_t b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

inline std::string strip_comment(const std::string& s)
{
    size_t p = s.find('#');
    return trim(p == std::string::npos ? s : s.substr(0, p));
}

inline std::vector<std::string> split(const std::string& s)
{
    std::stringstream ss(s);
    std::vector<std::string> out;
    std::string x;
    while (ss >> x) out.push_back(x);
    return out;
}

inline std::vector<std::string> clean_lines(const fs::path& path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path.string());
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        line = strip_comment(line);
        if (!line.empty()) lines.push_back(line);
    }
    return lines;
}

inline Vec3 add(Vec3 a, Vec3 b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
inline Vec3 sub(Vec3 a, Vec3 b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
inline Vec3 mul(Vec3 a, double x) { return {a[0] * x, a[1] * x, a[2] * x}; }
inline double norm(Vec3 a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }

inline Vec3 matvec(double f0, double f1, double f2, const std::array<Vec3, 3>& cell)
{
    return {
        f0 * cell[0][0] + f1 * cell[1][0] + f2 * cell[2][0],
        f0 * cell[0][1] + f1 * cell[1][1] + f2 * cell[2][1],
        f0 * cell[0][2] + f1 * cell[1][2] + f2 * cell[2][2],
    };
}

inline std::array<Vec3, 3> inverse_cell(const std::array<Vec3, 3>& cell)
{
    const Vec3& a = cell[0];
    const Vec3& b = cell[1];
    const Vec3& c = cell[2];
    double m00 = a[0], m01 = b[0], m02 = c[0];
    double m10 = a[1], m11 = b[1], m12 = c[1];
    double m20 = a[2], m21 = b[2], m22 = c[2];
    double det = m00 * (m11 * m22 - m12 * m21)
               - m01 * (m10 * m22 - m12 * m20)
               + m02 * (m10 * m21 - m11 * m20);
    if (std::fabs(det) < 1e-14) throw std::runtime_error("singular cell");
    double d = 1.0 / det;
    return {{
        {{ (m11 * m22 - m12 * m21) * d, -(m01 * m22 - m02 * m21) * d,  (m01 * m12 - m02 * m11) * d }},
        {{-(m10 * m22 - m12 * m20) * d,  (m00 * m22 - m02 * m20) * d, -(m00 * m12 - m02 * m10) * d }},
        {{ (m10 * m21 - m11 * m20) * d, -(m00 * m21 - m01 * m20) * d,  (m00 * m11 - m01 * m10) * d }},
    }};
}

inline Vec3 cart_to_frac(Vec3 r, const std::array<Vec3, 3>& inv)
{
    return {
        inv[0][0] * r[0] + inv[0][1] * r[1] + inv[0][2] * r[2],
        inv[1][0] * r[0] + inv[1][1] * r[1] + inv[1][2] * r[2],
        inv[2][0] * r[0] + inv[2][1] * r[1] + inv[2][2] * r[2],
    };
}

inline double minimum_image_distance(Vec3 a, Vec3 b,
                                     const std::array<Vec3, 3>& cell,
                                     const std::array<Vec3, 3>& inv)
{
    Vec3 f = cart_to_frac(sub(b, a), inv);
    f[0] -= std::round(f[0]);
    f[1] -= std::round(f[1]);
    f[2] -= std::round(f[2]);
    return norm(matvec(f[0], f[1], f[2], cell));
}

inline bool line_all_ints(const std::string& line)
{
    auto toks = split(line);
    if (toks.empty()) return false;
    for (const auto& t : toks) {
        char* end = nullptr;
        std::strtol(t.c_str(), &end, 10);
        if (end == t.c_str() || *end != '\0') return false;
    }
    return true;
}

inline SaveData read_save(const fs::path& path)
{
    auto lines = clean_lines(path);
    if (lines.size() < 10) throw std::runtime_error(path.string() + " is too short");
    size_t i = 0;
    ++i; // title
    double scale = std::stod(split(lines.at(i++)).at(0));

    std::vector<double> cell_flat;
    while (i < lines.size() && cell_flat.size() < 9) {
        for (const auto& t : split(lines.at(i))) cell_flat.push_back(std::stod(t) * scale);
        ++i;
    }
    if (cell_flat.size() < 9) throw std::runtime_error("failed to read cell from " + path.string());
    SaveData s;
    s.cell = {{
        {{cell_flat[0], cell_flat[1], cell_flat[2]}},
        {{cell_flat[3], cell_flat[4], cell_flat[5]}},
        {{cell_flat[6], cell_flat[7], cell_flat[8]}},
    }};

    s.metal_species = split(lines.at(i++));
    for (const auto& t : split(lines.at(i++))) s.metal_counts.push_back(std::stoi(t));
    if (s.metal_species.size() != s.metal_counts.size()) throw std::runtime_error("metal species/count mismatch");
    int num_metal = 0;
    for (int n : s.metal_counts) num_metal += n;

    s.interstitial_species = split(lines.at(i++));
    for (const auto& t : split(lines.at(i++))) s.interstitial_counts.push_back(std::stoi(t));
    if (s.interstitial_species.size() != s.interstitial_counts.size()) throw std::runtime_error("interstitial species/count mismatch");
    s.num_interstitial_sites = std::stoi(split(lines.at(i++)).at(0));

    if (i < lines.size()) {
        std::string low = lines[i];
        std::transform(low.begin(), low.end(), low.begin(), [](unsigned char c){ return (char)std::tolower(c); });
        if (low.find("shuffle") != std::string::npos) ++i;
    }

    bool cartesian = true;
    if (i < lines.size()) {
        auto toks = split(lines[i]);
        if (!toks.empty()) {
            char c = (char)std::tolower((unsigned char)toks[0][0]);
            if (c == 'd') cartesian = false;
            if (c == 'c' || c == 'k' || c == 'd') ++i;
        }
    }

    s.metal_types.reserve(num_metal);
    for (int m = 0; m < num_metal; ++m) {
        auto toks = split(lines.at(i++));
        if (toks.size() < 4) throw std::runtime_error("metal coordinate line missing type");
        s.metal_types.push_back(std::stoi(toks[3]));
    }

    s.interstitial_pos.reserve(s.num_interstitial_sites);
    s.interstitial_site_types.reserve(s.num_interstitial_sites);
    for (int site = 0; site < s.num_interstitial_sites; ++site) {
        auto toks = split(lines.at(i++));
        if (toks.size() < 4) throw std::runtime_error("interstitial coordinate line missing occupation");
        Vec3 r{std::stod(toks[0]), std::stod(toks[1]), std::stod(toks[2])};
        if (!cartesian) r = matvec(r[0], r[1], r[2], s.cell);
        s.interstitial_pos.push_back(r);
        s.interstitial_site_types.push_back(std::stoi(toks[3]));
    }
    return s;
}

inline PoscarData read_poscar(const fs::path& path)
{
    auto lines = clean_lines(path);
    if (lines.size() < 8) throw std::runtime_error(path.string() + " is too short");
    size_t i = 0;
    ++i; // title
    double scale = std::stod(split(lines.at(i++)).at(0));
    PoscarData p;
    for (int k = 0; k < 3; ++k) {
        auto toks = split(lines.at(i++));
        p.cell[k] = {{std::stod(toks[0]) * scale, std::stod(toks[1]) * scale, std::stod(toks[2]) * scale}};
    }
    if (line_all_ints(lines.at(i))) {
        auto counts = split(lines.at(i++));
        for (size_t k = 0; k < counts.size(); ++k) p.species.push_back("X" + std::to_string(k + 1));
        for (const auto& t : counts) p.counts.push_back(std::stoi(t));
    } else {
        p.species = split(lines.at(i++));
        for (const auto& t : split(lines.at(i++))) p.counts.push_back(std::stoi(t));
    }
    if (p.species.size() != p.counts.size()) throw std::runtime_error("POSCAR species/count mismatch");
    if (!lines.at(i).empty() && std::tolower((unsigned char)lines.at(i)[0]) == 's') ++i;
    char mode = (char)std::tolower((unsigned char)split(lines.at(i++)).at(0)[0]);
    bool cartesian = (mode == 'c' || mode == 'k');

    int total = 0;
    for (int n : p.counts) total += n;
    p.coords.reserve(total);
    for (int n = 0; n < total; ++n) {
        auto toks = split(lines.at(i++));
        Vec3 r{std::stod(toks[0]), std::stod(toks[1]), std::stod(toks[2])};
        if (cartesian) r = mul(r, scale);
        else r = matvec(r[0], r[1], r[2], p.cell);
        p.coords.push_back(r);
    }
    return p;
}

inline std::vector<Vec3> metal_coords_by_save_order(const SaveData& save, const PoscarData& contcar)
{
    std::map<std::string, int> offsets;
    int start = 0;
    for (size_t i = 0; i < contcar.species.size(); ++i) {
        offsets[contcar.species[i]] = start;
        start += contcar.counts[i];
    }
    std::map<std::string, int> used;
    std::vector<Vec3> coords;
    coords.reserve(save.metal_types.size());
    for (int t : save.metal_types) {
        if (t < 0 || t >= (int)save.metal_species.size()) throw std::runtime_error("metal type out of range");
        const std::string& sp = save.metal_species[t];
        if (!offsets.count(sp)) throw std::runtime_error("CONTCAR missing metallic species " + sp);
        int idx = offsets[sp] + used[sp]++;
        if (idx < 0 || idx >= (int)contcar.coords.size()) throw std::runtime_error("CONTCAR species count mismatch for " + sp);
        coords.push_back(contcar.coords[idx]);
    }
    return coords;
}

inline std::vector<fs::path> numbered_state_dirs(const fs::path& mcprocess)
{
    std::vector<fs::path> dirs;
    for (auto& e : fs::directory_iterator(mcprocess)) {
        if (!e.is_directory()) continue;
        std::string n = e.path().filename().string();
        if (!n.empty() && std::all_of(n.begin(), n.end(), [](unsigned char c){ return std::isdigit(c); })) {
            dirs.push_back(e.path());
        }
    }
    std::sort(dirs.begin(), dirs.end(), [](const fs::path& a, const fs::path& b){
        return std::stoi(a.filename().string()) < std::stoi(b.filename().string());
    });
    return dirs;
}

inline json read_meta(const fs::path& state_dir)
{
    fs::path p = state_dir / "meta.json";
    if (!fs::exists(p)) return json::object();
    try {
        std::ifstream in(p);
        json j;
        in >> j;
        return j;
    } catch (...) {
        return json::object();
    }
}

inline std::string read_energy(const fs::path& state_dir, const json& meta)
{
    fs::path epath = state_dir / "energy";
    if (fs::exists(epath)) {
        std::ifstream in(epath);
        std::string x;
        in >> x;
        if (!x.empty()) return x;
    }
    for (const auto& key : {"energy_final", "energy"}) {
        if (meta.contains(key)) {
            if (meta[key].is_number()) return std::to_string(meta[key].get<double>());
            if (meta[key].is_string()) return meta[key].get<std::string>();
        }
    }
    fs::path info = state_dir / "info.txt";
    if (fs::exists(info)) {
        std::ifstream in(info);
        std::string line;
        while (std::getline(in, line)) {
            if (line.find("energy") == 0) {
                std::replace(line.begin(), line.end(), '=', ' ');
                auto toks = split(line);
                if (!toks.empty()) return toks.back();
            }
        }
    }
    return "nan";
}

inline std::map<std::string, std::string> read_struc_stat_file(const fs::path& path)
{
    std::ifstream in(path);
    std::map<std::string, std::string> out;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        size_t tab = line.find('\t');
        if (tab != std::string::npos) {
            out[trim(line.substr(0, tab))] = trim(line.substr(tab + 1));
        } else {
            auto toks = split(line);
            if (toks.empty()) continue;
            std::string key = toks[0];
            std::string value;
            if (toks.size() > 1) value = toks[1];
            out[key] = value;
        }
    }
    return out;
}

} // namespace paipai_analysis

#endif
