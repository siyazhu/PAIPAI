#include <array>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "analysis_utils.h"

namespace fs = std::filesystem;
using namespace paipai_analysis;

static void usage()
{
    std::cerr << "Usage: packmc "
                 "[--root RUN_ROOT] "
                 "[--output mcprocess_summary.tsv] "
                 "[--contcar-dir mcprocess_CONTCARs] "
                 "[--tar mcprocess_CONTCARs.tar]\n";
}

static bool report_progress(size_t i, size_t total)
{
    return i == 1 || i == total || i % 100 == 0;
}

static void progress(const std::string& stage, size_t i, size_t total, const fs::path& item = {})
{
    if (!report_progress(i, total)) return;
    std::cout << "[" << stage << "] " << i << "/" << total;
    if (!item.empty()) std::cout << " " << item;
    std::cout << "\n";
}

static std::string trial_mode(const json& meta)
{
    for (const auto& key : {"move_type", "source"}) {
        if (meta.contains(key)) {
            if (meta[key].is_string()) return meta[key].get<std::string>();
            return meta[key].dump();
        }
    }
    return "unknown";
}

static void put_octal(char* dst, size_t n, uint64_t value)
{
    std::ostringstream ss;
    ss << std::oct << value;
    std::string s = ss.str();
    std::memset(dst, '0', n);
    if (s.size() + 1 > n) s = s.substr(s.size() + 1 - n);
    std::memcpy(dst + n - s.size() - 1, s.data(), s.size());
    dst[n - 1] = '\0';
}

static void write_tar_header(std::ofstream& out, const std::string& name, uint64_t size, char type)
{
    std::array<char, 512> h{};
    std::string nm = name;
    if (nm.size() > 100) nm = nm.substr(nm.size() - 100);
    std::memcpy(h.data(), nm.data(), nm.size());
    put_octal(h.data() + 100, 8, 0644);
    put_octal(h.data() + 108, 8, 0);
    put_octal(h.data() + 116, 8, 0);
    put_octal(h.data() + 124, 12, size);
    put_octal(h.data() + 136, 12, 0);
    std::memset(h.data() + 148, ' ', 8);
    h[156] = type;
    std::memcpy(h.data() + 257, "ustar", 5);
    std::memcpy(h.data() + 263, "00", 2);

    unsigned int sum = 0;
    for (unsigned char c : h) sum += c;
    put_octal(h.data() + 148, 8, sum);
    h[155] = ' ';
    out.write(h.data(), h.size());
}

static void write_tar(const fs::path& dir, const fs::path& tar_path)
{
    std::ofstream out(tar_path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot write " + tar_path.string());

    std::string root_name = dir.filename().string();
    write_tar_header(out, root_name + "/", 0, '5');

    std::vector<fs::path> files;
    for (auto& e : fs::directory_iterator(dir)) {
        if (e.is_regular_file()) files.push_back(e.path());
    }
    std::sort(files.begin(), files.end());

    std::array<char, 8192> buf{};
    for (size_t fi = 0; fi < files.size(); ++fi) {
        const auto& file = files[fi];
        progress("write tar", fi + 1, files.size(), file.filename());
        uint64_t size = (uint64_t)fs::file_size(file);
        write_tar_header(out, root_name + "/" + file.filename().string(), size, '0');
        std::ifstream in(file, std::ios::binary);
        uint64_t remaining = size;
        while (remaining > 0) {
            size_t chunk = (size_t)std::min<uint64_t>(buf.size(), remaining);
            in.read(buf.data(), chunk);
            out.write(buf.data(), chunk);
            remaining -= chunk;
        }
        size_t pad = (512 - (size % 512)) % 512;
        if (pad) {
            std::array<char, 512> zeros{};
            out.write(zeros.data(), pad);
        }
    }
    std::array<char, 1024> zeros{};
    out.write(zeros.data(), zeros.size());
}

int main(int argc, char** argv)
{
    fs::path root = ".";
    fs::path output = "mcprocess_summary.tsv";
    fs::path contcar_dir = "mcprocess_CONTCARs";
    fs::path tar_path = "mcprocess_CONTCARs.tar";

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--root" && i + 1 < argc) root = argv[++i];
        else if (a == "--output" && i + 1 < argc) output = argv[++i];
        else if (a == "--contcar-dir" && i + 1 < argc) contcar_dir = argv[++i];
        else if (a == "--tar" && i + 1 < argc) tar_path = argv[++i];
        else if (a == "-h" || a == "--help") { usage(); return 0; }
        else {
            std::cerr << "Unknown argument: " << a << "\n";
            usage();
            return 2;
        }
    }

    try {
        root = fs::absolute(root);
        fs::path mcprocess = root / "mcprocess";
        fs::path base = root;
        if (!fs::exists(mcprocess) && fs::is_directory(root)) {
            auto direct_states = numbered_state_dirs(root);
            if (!direct_states.empty()) {
                mcprocess = root;
                base = root.parent_path();
            }
        }
        if (!output.is_absolute()) output = base / output;
        if (!contcar_dir.is_absolute()) contcar_dir = base / contcar_dir;
        if (!tar_path.is_absolute()) tar_path = base / tar_path;

        auto states = numbered_state_dirs(mcprocess);
        if (states.empty()) throw std::runtime_error("no numbered state directories found in " + mcprocess.string());

        if (fs::exists(contcar_dir)) fs::remove_all(contcar_dir);
        fs::create_directories(contcar_dir);
        std::cout << "[collect CONTCAR] collecting " << states.size() << " states\n";
        for (size_t si = 0; si < states.size(); ++si) {
            const auto& state = states[si];
            progress("collect CONTCAR", si + 1, states.size(), state.filename());
            fs::path src = state / "CONTCAR";
            if (fs::exists(src)) fs::copy_file(src, contcar_dir / ("CONTCAR" + state.filename().string()), fs::copy_options::overwrite_existing);
        }
        std::cout << "[write tar] packing CONTCAR files\n";
        write_tar(contcar_dir, tar_path);

        std::vector<std::string> columns;
        std::set<std::string> seen;
        std::map<std::string, std::map<std::string, std::string>> values_by_state;
        std::cout << "[scan stats] reading struc_* files\n";
        for (size_t si = 0; si < states.size(); ++si) {
            const auto& state = states[si];
            progress("scan stats", si + 1, states.size(), state.filename());
            std::map<std::string, std::string> vals;
            std::vector<fs::path> stat_files;
            for (auto& e : fs::directory_iterator(state)) {
                if (e.is_regular_file() && e.path().filename().string().find("struc_") == 0) {
                    stat_files.push_back(e.path());
                }
            }
            std::sort(stat_files.begin(), stat_files.end());
            for (const auto& sf : stat_files) {
                for (const auto& kv : read_struc_stat_file(sf)) {
                    std::string col = sf.filename().string() + ":" + kv.first;
                    vals[col] = kv.second;
                    if (!seen.count(col)) {
                        seen.insert(col);
                        columns.push_back(col);
                    }
                }
            }
            values_by_state[state.filename().string()] = vals;
        }

        std::ofstream table(output);
        if (!table) throw std::runtime_error("cannot write " + output.string());
        table << "step\tenergy\ttrial_mode";
        for (const auto& c : columns) table << "\t" << c;
        table << "\n";

        std::cout << "[write summary] writing TSV rows\n";
        for (size_t si = 0; si < states.size(); ++si) {
            const auto& state = states[si];
            progress("write summary", si + 1, states.size(), state.filename());
            json meta = read_meta(state);
            std::string name = state.filename().string();
            table << name << "\t" << read_energy(state, meta) << "\t" << trial_mode(meta);
            const auto& vals = values_by_state[name];
            for (const auto& c : columns) {
                auto it = vals.find(c);
                table << "\t" << (it == vals.end() ? "nan" : it->second);
            }
            table << "\n";
        }

        std::cout << "Wrote " << output << "\n";
        std::cout << "Collected CONTCAR files in " << contcar_dir << "\n";
        std::cout << "Wrote " << tar_path << "\n";
    } catch (const std::exception& e) {
        std::cerr << "packmc error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
