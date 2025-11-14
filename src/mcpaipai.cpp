#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <thread>
#include <filesystem>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <optional>

#include "structure.h"   // your Structure class
#include "json.hpp"      // nlohmann::json single-header

namespace fs = std::filesystem;
using json = nlohmann::json;

/* ---------- CLI configuration ---------- */
struct Args {
    std::string input_struc;     // initial structure file
    int workers    = 4;          // number of fast worker slots (and python fast_worker processes)
    int steps      = 1000;       // number of MC trial steps (proposals)
    double temp    = 0.001;      // Metropolis temperature factor
    int p_swap_metal  = 70;
    int p_swap_inter  = 30;
    int p_exch_metal  = 0;
    int p_exch_inter  = 0;
};

void print_help(const char* prog){
    std::cout << "Usage: " << prog << " INPUT_STRUCTURE "
              << "[--workers N] [--steps K] [--temp T] "
              << "[--p-swap-metal P] [--p-swap-inter P] "
              << "[--p-exch-metal P] [--p-exch-inter P]\n";
}

Args parse_args(int argc, char** argv){
    if (argc < 2) { print_help(argv[0]); std::exit(2); }
    Args a; a.input_struc = argv[1];
    for (int i=2; i<argc; ++i){
        std::string s = argv[i];
        auto need=[&](const char* name){ if (i+1>=argc){ std::cerr<<"Missing "<<name<<"\n"; std::exit(2);} };
        if (s=="--workers"){ need("--workers"); a.workers=std::max(1, std::stoi(argv[++i])); }
        else if (s=="--steps"){ need("--steps"); a.steps=std::max(1, std::stoi(argv[++i])); }
        else if (s=="--temp"){ need("--temp"); a.temp=std::stod(argv[++i]); }
        else if (s=="--p-swap-metal"){ need("--p-swap-metal"); a.p_swap_metal=std::max(0, std::stoi(argv[++i])); }
        else if (s=="--p-swap-inter"){ need("--p-swap-inter"); a.p_swap_inter=std::max(0, std::stoi(argv[++i])); }
        else if (s=="--p-exch-metal"){ need("--p-exch-metal"); a.p_exch_metal=std::max(0, std::stoi(argv[++i])); }
        else if (s=="--p-exch-inter"){ need("--p-exch-inter"); a.p_exch_inter=std::max(0, std::stoi(argv[++i])); }
        else { std::cerr<<"Unknown arg: "<<s<<"\n"; print_help(argv[0]); std::exit(2); }
    }
    int sum = a.p_swap_metal + a.p_swap_inter + a.p_exch_metal + a.p_exch_inter;
    if (sum <= 0){
        std::cerr << "MC move probabilities are incorrect. Please check input parameters.\n";
        std::exit(2);
    }
    return a;
}

/* ---------- small helpers ---------- */

// Read a plain text energy (one double).
std::optional<double> read_energy_text(const fs::path& f){
    std::ifstream ifs(f);
    if(!ifs) return std::nullopt;
    double e; ifs >> e; if(!ifs) return std::nullopt;
    return e;
}

// Copy small file (overwrite).
void copy_file_overwrite(const fs::path& src, const fs::path& dst) {
    if (!fs::exists(src)) return;
    fs::create_directories(dst.parent_path());
    std::error_code ec;
    fs::copy_file(src, dst, fs::copy_options::overwrite_existing, ec);
}

// Atomic-ish integer counter for mc_count (very simple).
int increment_mc_counter(const fs::path& root) {
    fs::path ctr_dir = root / "counters";
    fs::path mc_ctr  = ctr_dir / "mc_count";
    fs::create_directories(ctr_dir);

    int idx = 0;
    if (fs::exists(mc_ctr)) {
        std::ifstream ifs(mc_ctr);
        if (ifs) ifs >> idx;
    }
    int new_idx = idx + 1;
    {
        std::ofstream ofs(mc_ctr);
        ofs << new_idx << "\n";
    }
    return new_idx;
}

// Archive an accepted state into mcprocess/000001/...
void archive_mc_accept(const fs::path& root,
                       const std::string& task_id,
                       double E_final)
{
    fs::path out_dir = root / "refine_outbox" / task_id;
    fs::path mc_root = root / "mcprocess";
    fs::create_directories(mc_root);

    int new_idx = increment_mc_counter(root);

    std::ostringstream oss;
    oss << std::setw(6) << std::setfill('0') << new_idx;
    fs::path mc_dir = mc_root / oss.str();
    fs::create_directories(mc_dir);

    copy_file_overwrite(out_dir / "CONTCAR",  mc_dir / "CONTCAR");
    copy_file_overwrite(out_dir / "SAVE",     mc_dir / "SAVE");
    copy_file_overwrite(out_dir / "meta.json",mc_dir / "meta.json");

    std::ofstream info(mc_dir / "info.txt");
    if (info) {
        info << "task_id = " << task_id << "\n";
        info << "E_final = " << std::setprecision(12) << E_final << "\n";
    }

    std::cout << "[MC] accepted task " << task_id
              << ", archived to " << mc_dir << "\n";
}

// Metropolis acceptance
bool Accept(double Eold, double Enew, double temp,
            std::mt19937_64& rng)
{
    if (Enew <= Eold) return true;
    double p = std::exp(-(Enew - Eold)/temp);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double u = dist(rng);
    return (u < p);
}

// Generate one candidate for a given fast slot: uses current SAVE.
void generate_candidate_for_slot(int slot,
                                 const Args& cfg,
                                 const fs::path& root,
                                 Structure& struc,
                                 int SUMP)
{
    const int P1 = cfg.p_swap_metal;
    const int P2 = cfg.p_swap_inter;
    const int P3 = cfg.p_exch_metal;
    const int P4 = cfg.p_exch_inter;

    // 1) Load current accepted state from SAVE.
    struc.readstruc("SAVE");

    // 2) Random MC move.
    int r = std::rand() % SUMP;
    if (r < P1) {
        // swapMetal
        int a = std::rand() % struc.num_metallic_atoms;
        int b = std::rand() % struc.num_metallic_atoms;
        while (struc.atomtype[a] == struc.atomtype[b]) {
            a = std::rand() % struc.num_metallic_atoms;
            b = std::rand() % struc.num_metallic_atoms;
        }
        struc.swapMetal(a,b);
    } else if ((r -= P1) < P2) {
        // swapInterstitial
        int a = std::rand() % struc.num_interstitial;
        while (struc.interstitial_postype[a] == -1)
            a = (a+1) % struc.num_interstitial;
        int b = std::rand() % struc.num_interstitial;
        while (struc.interstitial_postype[b] == struc.interstitial_postype[a])
            b = std::rand() % struc.num_interstitial;
        struc.swapInterstitial(a,b);
    } else if ((r -= P2) < P3) {
        // TODO: exchangeMetal
        // struc.exchangeMetal(...);
    } else {
        // TODO: exchangeInterstitial
        // struc.exchangeInterstitial(...);
    }

    // 3) Dump candidate directly into fast/POSCARk, fast/SAVEk
    fs::path fast_dir = root / "fast";
    fs::create_directories(fast_dir);

    fs::path pos = fast_dir / ("POSCAR" + std::to_string(slot));
    fs::path sav = fast_dir / ("SAVE"   + std::to_string(slot));

    struc.outputvasp(pos.string().c_str());
    struc.outputsave(sav.string().c_str());

    // 4) Trigger fast worker by touching .go_k
    fs::path gof = fast_dir / (".go_" + std::to_string(slot));
    std::ofstream ofs(gof); // touch
}

/* ---------- process a single slow report ---------- */
bool process_report_file(const fs::path& root,
                         const fs::path& rep_path,
                         double& current_E,
                         bool& have_state,
                         int& mc_steps,
                         int& accept_count,
                         double temp,
                         std::mt19937_64& rng,
                         std::ofstream& log)
{
    // Parse JSON
    json j;
    {
        std::ifstream ifs(rep_path);
        if (!ifs) {
            std::cerr << "[WARN] cannot open report " << rep_path << "\n";
            fs::remove(rep_path);
            return false;
        }
        try {
            ifs >> j;
        } catch (...) {
            std::cerr << "[WARN] bad JSON in " << rep_path << "\n";
            fs::remove(rep_path);
            return false;
        }
    }

    std::string status = j.value("status", "");
    if (status == "error") {
        std::cerr << "[slow] ERROR in report " << rep_path.filename()
                  << ": " << j.value("error", std::string("<no_msg>")) << "\n";
        fs::remove(rep_path);
        return false;
    }

    std::string task_id = j.value("task_id", rep_path.stem().string());
    double E_final = j.value("energy_final", std::numeric_limits<double>::infinity());

    if (!std::isfinite(E_final)) {
        std::cerr << "[WARN] invalid energy_final in report " << rep_path << "\n";
        fs::remove(rep_path);
        return false;
    }

    fs::path out_dir = root / "refine_outbox" / task_id;

    if (!have_state) {
        // First-ever refined structure: treat as initial state.
        have_state = true;
        current_E  = E_final;
        accept_count = 0;
        mc_steps = 0;

        // Update global SAVE, CONTCAR
        copy_file_overwrite(out_dir / "SAVE",    root / "SAVE");
        copy_file_overwrite(out_dir / "CONTCAR", root / "CONTCAR");

        log << "INITIAL_STATE task_id=" << task_id
            << " E = " << std::setprecision(12) << E_final << "\n";
        log.flush();

        fs::remove(rep_path);
        return true;
    }

    // Normal MC proposal
    ++mc_steps;
    bool accept = Accept(current_E, E_final, temp, rng);
    log << "STEP " << mc_steps
        << " proposal task_id=" << task_id
        << " E_new=" << std::setprecision(12) << E_final
        << " E_old=" << current_E
        << " -> "   << (accept ? "ACCEPT" : "REJECT")
        << "\n";
    log.flush();

    if (accept) {
        ++accept_count;
        current_E = E_final;
        // Update global SAVE & CONTCAR
        copy_file_overwrite(out_dir / "SAVE",  root / "SAVE");
        copy_file_overwrite(out_dir / "CONTCAR", root /"CONTCAR");
        // Archive
        archive_mc_accept(root, task_id, E_final);
    }

    fs::remove(rep_path);
    return true;
}

/* ---------- main ---------- */

int main(int argc, char** argv){
    std::ios::sync_with_stdio(false);
    std::srand((unsigned)std::time(nullptr));
    std::mt19937_64 rng((unsigned)std::time(nullptr));

    Args cfg = parse_args(argc, argv);
    const int SUMP = cfg.p_swap_metal + cfg.p_swap_inter +
                     cfg.p_exch_metal + cfg.p_exch_inter;

    fs::path ROOT = ".";

    // Directories expected by workers
    fs::create_directories(ROOT / "fast");
    fs::create_directories(ROOT / "reports");
    fs::create_directories(ROOT / "refine_outbox");
    fs::create_directories(ROOT / "waiting_pool");
    fs::create_directories(ROOT / "waiting_work");
    fs::create_directories(ROOT / "counters");
    fs::create_directories(ROOT / "mcprocess");

    // Load initial structure and output an initial SAVE (used as MC seed)
    Structure struc;
    struc.readstruc(cfg.input_struc.c_str());
    struc.outputsave("SAVE");

    std::ofstream log("mc.log");
    log << "# MC with waiting_pool, fast=" << cfg.workers
        << " steps=" << cfg.steps
        << " temp=" << cfg.temp << "\n";

    bool have_state = false;
    double current_E = 0.0;
    int mc_steps = 0;
    int accept_count = 0;

    // Main loop: keep feeding fast workers and consuming slow reports
    while (mc_steps < cfg.steps) {
        // 1) For each fast slot, if .go_k does not exist, schedule a new candidate.
        for (int k = 1; k <= cfg.workers; ++k) {
            fs::path gof = ROOT / "fast" / (".go_" + std::to_string(k));
            if (fs::exists(gof)) continue;  // slot is busy

            // generate a new candidate and trigger fast worker k
            generate_candidate_for_slot(k, cfg, ROOT, struc, SUMP);
        }

        // 2) Poll reports directory for new slow results
        bool processed_any = false;
        fs::path reports_dir = ROOT / "reports";
        for (auto& entry : fs::directory_iterator(reports_dir)) {
            if (!entry.is_regular_file()) continue;
            fs::path rep_path = entry.path();
            // Only process .json files
            if (rep_path.extension() != ".json") continue;

            bool ok = process_report_file(
                ROOT, rep_path,
                current_E, have_state,
                mc_steps, accept_count,
                cfg.temp, rng, log
            );
            if (ok) processed_any = true;
            if (mc_steps >= cfg.steps) break;
        }

        if (mc_steps >= cfg.steps) break;

        if (!processed_any) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }

    log << "# Finished. MC steps = " << mc_steps
        << ", accepted = " << accept_count << "\n";
    log.close();

    std::cout << "MC finished: steps=" << mc_steps
              << " accepted=" << accept_count << "\n";
    return 0;
}
