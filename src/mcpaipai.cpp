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
#include <algorithm>
#include <limits>
#include <ctime>

#include "structure.h"   // your Structure class
#include "json.hpp"      // nlohmann::json single-header
#include "prefast.h"

namespace fs = std::filesystem;
using json = nlohmann::json;
using paipai_prefast::PrefastConfig;
using paipai_prefast::PrefastModel;
using paipai_prefast::SparseDescriptor;
using paipai_prefast::TrialScore;

/* ---------- CLI configuration ---------- */
struct Args {
    std::string input_struc;     // initial structure file
    int workers    = 4;          // number of fast worker slots (and python fast_worker processes)
    int steps      = 1000;       // number of MC trial steps (proposals)
    double temp    = 0.001;      // Metropolis temperature factor
    int p_swap_metal  = 70;
    int p_swap_inter  = 30;
    int p_hop_inter = 0;
    int p_cluster_inter = 0;
    int p_exch_metal  = 0;
    int p_exch_inter  = 0;
    double intsite_neighbor_cutoff = 3.5;  // cutoff for interstitial-site metal-neighbor map (Angstrom)
    double intsite_hop_cutoff = 3.0;       // cutoff for interstitial-site local hop graph (Angstrom)
    std::string mode = "search";          // "search" = fast+slow pool; "finiteT" = one direct slow task at a time
    double interstitial_site_cutoff = 1.5; // max relaxed interstitial distance to assigned site (Angstrom)
    std::string resume_state_dir;          // existing SAVE/CONTCAR/energy seed
    bool continue_run = false;             // continue from ./mcprocess latest accepted state
    PrefastConfig prefast;
};

void print_help(const char* prog){
    std::cout << "Usage: " << prog << " [INPUT_STRUCTURE] "
              << "[--workers N] [--steps K] [--temp T] "
              << "[--p-swap-metal P] [--p-swap-inter P] "
              << "[--p-hop-inter P] [--p-cluster-inter P] "
              << "[--p-exch-metal P] [--p-exch-inter P] "
              << "[--intsite-neighbor-cutoff R] [--intsite-hop-cutoff R] "
              << "[--mode search|finiteT] [--finiteT] "
              << "[--interstitial-site-cutoff R] "
              << "[--resume-state DIR] [--continue] "
              << "[--prefast on|off] [--prefast-candidates-per-slot N] "
              << "[--prefast-warmup-steps N] "
              << "[--prefast-basis ref-dz] [--prefast-nshells N] "
              << "[--prefast-peak-scan-cutoff R] [--prefast-peak-tol R] "
              << "[--prefast-sigma-small R] [--prefast-sigma-large R] "
              << "[--prefast-cutoff-margin R] [--prefast-learning-rate LR] "
              << "[--prefast-lms-epsilon EPS] [--prefast-weight-decay R] "
              << "[--prefast-diagnostics off|summary|updates|full]\n";
}

Args parse_args(int argc, char** argv){
    if (argc < 2) { print_help(argv[0]); std::exit(2); }
    Args a;
    int i = 1;
    if (i < argc && std::string(argv[i]).rfind("-", 0) != 0) {
        a.input_struc = argv[i++];
    }
    for (; i<argc; ++i){
        std::string s = argv[i];
        auto need=[&](const char* name){ if (i+1>=argc){ std::cerr<<"Missing "<<name<<"\n"; std::exit(2);} };
        if (s=="--workers"){ need("--workers"); a.workers=std::max(1, std::stoi(argv[++i])); }
        else if (s=="--steps"){ need("--steps"); a.steps=std::max(1, std::stoi(argv[++i])); }
        else if (s=="--temp"){ need("--temp"); a.temp=std::stod(argv[++i]); }
        else if (s=="--p-swap-metal"){ need("--p-swap-metal"); a.p_swap_metal=std::max(0, std::stoi(argv[++i])); }
        else if (s=="--p-swap-inter"){ need("--p-swap-inter"); a.p_swap_inter=std::max(0, std::stoi(argv[++i])); }
        else if (s=="--p-hop-inter"){ need("--p-hop-inter"); a.p_hop_inter=std::max(0, std::stoi(argv[++i])); }
        else if (s=="--p-cluster-inter"){ need("--p-cluster-inter"); a.p_cluster_inter=std::max(0, std::stoi(argv[++i])); }
        else if (s=="--p-exch-metal"){ need("--p-exch-metal"); a.p_exch_metal=std::max(0, std::stoi(argv[++i])); }
        else if (s=="--p-exch-inter"){ need("--p-exch-inter"); a.p_exch_inter=std::max(0, std::stoi(argv[++i])); }
        else if (s=="--intsite-neighbor-cutoff"){ need("--intsite-neighbor-cutoff"); a.intsite_neighbor_cutoff=std::stod(argv[++i]); }
        else if (s=="--intsite-hop-cutoff"){ need("--intsite-hop-cutoff"); a.intsite_hop_cutoff=std::stod(argv[++i]); }
        else if (s=="--mode"){ need("--mode"); a.mode=argv[++i]; }
        else if (s=="--finiteT"){ a.mode="finiteT"; }
        else if (s=="--interstitial-site-cutoff"){ need("--interstitial-site-cutoff"); a.interstitial_site_cutoff=std::stod(argv[++i]); }
        else if (s=="--resume-state"){ need("--resume-state"); a.resume_state_dir=argv[++i]; }
        else if (s=="--continue"){ a.continue_run=true; }
        else if (s=="--prefast"){
            need("--prefast");
            std::string v = argv[++i];
            if (v == "on" || v == "true" || v == "1") a.prefast.enabled = true;
            else if (v == "off" || v == "false" || v == "0") a.prefast.enabled = false;
            else { std::cerr << "--prefast must be on or off\n"; std::exit(2); }
        }
        else if (s=="--prefast-candidates-per-slot"){ need("--prefast-candidates-per-slot"); a.prefast.candidates_per_slot=std::max(1, std::stoi(argv[++i])); }
        else if (s=="--prefast-warmup-steps"){ need("--prefast-warmup-steps"); a.prefast.warmup_steps=std::max(0, std::stoi(argv[++i])); }
        else if (s=="--prefast-basis"){ need("--prefast-basis"); a.prefast.basis=argv[++i]; }
        else if (s=="--prefast-nshells"){ need("--prefast-nshells"); a.prefast.nshells=std::max(1, std::stoi(argv[++i])); }
        else if (s=="--prefast-peak-scan-cutoff"){ need("--prefast-peak-scan-cutoff"); a.prefast.peak_scan_cutoff=std::stod(argv[++i]); }
        else if (s=="--prefast-peak-tol"){ need("--prefast-peak-tol"); a.prefast.peak_tol=std::stod(argv[++i]); }
        else if (s=="--prefast-sigma-small"){ need("--prefast-sigma-small"); a.prefast.sigma_small=std::stod(argv[++i]); }
        else if (s=="--prefast-sigma-large"){ need("--prefast-sigma-large"); a.prefast.sigma_large=std::stod(argv[++i]); }
        else if (s=="--prefast-cutoff-margin"){ need("--prefast-cutoff-margin"); a.prefast.cutoff_margin=std::stod(argv[++i]); }
        else if (s=="--prefast-learning-rate"){ need("--prefast-learning-rate"); a.prefast.learning_rate=std::stod(argv[++i]); }
        else if (s=="--prefast-lms-epsilon"){ need("--prefast-lms-epsilon"); a.prefast.lms_epsilon=std::stod(argv[++i]); }
        else if (s=="--prefast-weight-decay"){ need("--prefast-weight-decay"); a.prefast.weight_decay=std::stod(argv[++i]); }
        else if (s=="--prefast-diagnostics"){
            need("--prefast-diagnostics");
            a.prefast.diagnostics=argv[++i];
        }
        else { std::cerr<<"Unknown arg: "<<s<<"\n"; print_help(argv[0]); std::exit(2); }
    }
    int sum = a.p_swap_metal + a.p_swap_inter + a.p_hop_inter + a.p_cluster_inter +
              a.p_exch_metal + a.p_exch_inter;
    if (sum <= 0){
        std::cerr << "MC move probabilities are incorrect. Please check input parameters.\n";
        std::exit(2);
    }
    if (a.mode != "search" && a.mode != "finiteT") {
        std::cerr << "Unknown --mode: " << a.mode
                  << " (valid: search, finiteT)\n";
        std::exit(2);
    }
    if (a.interstitial_site_cutoff <= 0.0) {
        std::cerr << "--interstitial-site-cutoff must be positive.\n";
        std::exit(2);
    }
    if (a.intsite_hop_cutoff <= 0.0) {
        std::cerr << "--intsite-hop-cutoff must be positive.\n";
        std::exit(2);
    }
    if (a.continue_run && !a.resume_state_dir.empty()) {
        std::cerr << "--continue and --resume-state cannot be used together.\n";
        std::exit(2);
    }
    if (a.prefast.peak_scan_cutoff <= 0.0 ||
        a.prefast.peak_tol <= 0.0 ||
        a.prefast.sigma_small <= 0.0 ||
        a.prefast.sigma_large <= 0.0 ||
        a.prefast.cutoff_margin < 0.0 ||
        a.prefast.learning_rate < 0.0 ||
        a.prefast.lms_epsilon <= 0.0 ||
        a.prefast.weight_decay < 0.0) {
        std::cerr << "Invalid prefast parameter value.\n";
        std::exit(2);
    }
    if (a.prefast.diagnostics != "off" &&
        a.prefast.diagnostics != "summary" &&
        a.prefast.diagnostics != "updates" &&
        a.prefast.diagnostics != "full") {
        std::cerr << "--prefast-diagnostics must be off, summary, updates, or full.\n";
        std::exit(2);
    }
    if (a.resume_state_dir.empty() && !a.continue_run && a.input_struc.empty()) {
        std::cerr << "INPUT_STRUCTURE is required unless --resume-state or --continue is used.\n";
        print_help(argv[0]);
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

std::optional<double> read_energy_from_meta(const fs::path& f) {
    if (!fs::exists(f)) return std::nullopt;
    try {
        std::ifstream ifs(f);
        json j;
        ifs >> j;
        if (j.contains("energy_final")) {
            return j.at("energy_final").get<double>();
        }
        if (j.contains("E_final")) {
            return j.at("E_final").get<double>();
        }
        if (j.contains("energy")) {
            return j.at("energy").get<double>();
        }
    } catch (...) {
        return std::nullopt;
    }
    return std::nullopt;
}

std::optional<double> read_resume_energy(const fs::path& resume_dir) {
    if (auto e = read_energy_text(resume_dir / "energy")) return e;
    if (auto e = read_energy_from_meta(resume_dir / "meta.json")) return e;
    return std::nullopt;
}

std::optional<int> read_last_logged_mc_step(const fs::path& log_path) {
    std::ifstream in(log_path);
    if (!in) return std::nullopt;
    int last = 0;
    std::string word;
    while (in >> word) {
        if (word != "STEP") continue;
        int n = 0;
        if (in >> n) last = std::max(last, n);
    }
    return last > 0 ? std::optional<int>(last) : std::nullopt;
}

bool looks_like_state_dir(const fs::path& dir) {
    return fs::exists(dir / "SAVE") &&
           fs::exists(dir / "CONTCAR") &&
           (fs::exists(dir / "energy") || fs::exists(dir / "meta.json"));
}

bool is_all_digits(const std::string& s) {
    if (s.empty()) return false;
    for (char ch : s) {
        if (ch < '0' || ch > '9') return false;
    }
    return true;
}

fs::path resolve_resume_state_dir(const fs::path& requested) {
    if (looks_like_state_dir(requested)) {
        return requested;
    }

    fs::path best;
    std::string best_name;
    if (!fs::exists(requested) || !fs::is_directory(requested)) {
        return requested;
    }

    std::error_code ec;
    for (auto& entry : fs::directory_iterator(requested, ec)) {
        if (ec) break;
        if (!entry.is_directory()) continue;
        std::string name = entry.path().filename().string();
        if (!is_all_digits(name)) continue;
        if (!looks_like_state_dir(entry.path())) continue;
        if (best.empty() || name > best_name) {
            best = entry.path();
            best_name = name;
        }
    }

    return best.empty() ? requested : best;
}

fs::path latest_mcprocess_state_dir(const fs::path& root) {
    return resolve_resume_state_dir(root / "mcprocess");
}

// Copy small file (overwrite).
void copy_file_overwrite(const fs::path& src, const fs::path& dst) {
    if (!fs::exists(src)) return;
    fs::create_directories(dst.parent_path());
    std::error_code ec;
    fs::copy_file(src, dst, fs::copy_options::overwrite_existing, ec);
}

bool initialize_from_resume_state(const fs::path& root,
                                  const fs::path& resume_dir,
                                  double& current_E,
                                  bool& have_state,
                                  int& mc_steps,
                                  int& accept_count,
                                  std::ofstream& log)
{
    fs::path state_dir = resolve_resume_state_dir(resume_dir);
    fs::path save = state_dir / "SAVE";
    fs::path contcar = state_dir / "CONTCAR";
    auto energy = read_resume_energy(state_dir);

    if (!fs::exists(save)) {
        std::cerr << "[resume] missing SAVE in " << state_dir << "\n";
        return false;
    }
    if (!fs::exists(contcar)) {
        std::cerr << "[resume] missing CONTCAR in " << state_dir << "\n";
        return false;
    }
    if (!energy || !std::isfinite(*energy)) {
        std::cerr << "[resume] missing valid energy in " << state_dir
                  << " (expected energy or meta.json)\n";
        return false;
    }

    copy_file_overwrite(save, root / "SAVE");
    copy_file_overwrite(contcar, root / "CONTCAR");
    copy_file_overwrite(state_dir / "REFERENCE_SAVE", root / "REFERENCE_SAVE");

    have_state = true;
    current_E = *energy;
    mc_steps = 0;
    accept_count = 0;

    log << "RESUME_STATE requested=" << resume_dir.string()
        << " resolved=" << state_dir.string()
        << " E = " << std::setprecision(12) << current_E << "\n";
    log.flush();

    std::cout << "[resume] initialized from " << state_dir
              << " E_current = " << std::setprecision(12) << current_E << "\n";
    return true;
}


bool is_ignorable_runtime_name(const std::string& name) {
    if (name.empty()) return true;
    if (name == "." || name == "..") return true;
    if (name.rfind(".tmp_", 0) == 0) return true;
    if (name.size() >= 5 && name.substr(name.size() - 5) == ".lock") return true;
    return false;
}

bool dir_has_non_tmp_entries(const fs::path& p) {
    if (!fs::exists(p)) return false;
    std::error_code ec;
    for (auto& e : fs::directory_iterator(p, ec)) {
        if (ec) break;
        std::string name = e.path().filename().string();
        if (is_ignorable_runtime_name(name)) continue;
        return true;
    }
    return false;
}

// In finiteT mode we want a strict one-proposal-at-a-time chain.  A task is
// pending if it is in the slow queue, claimed by a worker, finished but not yet
// consumed, or has a report waiting for the master.
bool has_pending_slow_task(const fs::path& root) {
    if (dir_has_non_tmp_entries(root / "waiting_pool")) return true;
    if (dir_has_non_tmp_entries(root / "waiting_work")) return true;
    if (dir_has_non_tmp_entries(root / "refine_outbox")) return true;
    if (dir_has_non_tmp_entries(root / "reports")) return true;
    return false;
}

// Remove files/directories associated with a processed slow-worker task.
// After the master consumes reports/<task_id>.json, the corresponding
// refine_outbox/<task_id> directory is no longer needed. Accepted states are
// copied to SAVE/CONTCAR and archived to mcprocess before this cleanup.
void cleanup_processed_task(const fs::path& root,
                            const std::string& task_id,
                            const fs::path& rep_path,
                            bool remove_outbox = true)
{
    std::error_code ec;
    fs::remove(rep_path, ec);

    if (remove_outbox && !task_id.empty()) {
        fs::path out_dir = root / "refine_outbox" / task_id;
        if (fs::exists(out_dir)) {
            fs::remove_all(out_dir, ec);
        }
    }
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

int read_counter_value(const fs::path& root, const std::string& name) {
    fs::path ctr = root / "counters" / name;
    int value = 0;
    if (fs::exists(ctr)) {
        std::ifstream ifs(ctr);
        if (ifs) ifs >> value;
    }
    return value;
}

int count_busy_fast_slots(const fs::path& root, int workers) {
    int busy = 0;
    for (int k = 1; k <= workers; ++k) {
        fs::path gof = root / "fast" / (".go_" + std::to_string(k));
        if (fs::exists(gof)) ++busy;
    }
    return busy;
}

int run_progress_steps(const fs::path& root,
                       const std::string& mode,
                       int mc_steps)
{
    if (mode == "search") {
        return read_counter_value(root, "fast_count");
    }
    return mc_steps;
}

bool reached_step_budget(const fs::path& root,
                         const Args& cfg,
                         int mc_steps)
{
    return run_progress_steps(root, cfg.mode, mc_steps) >= cfg.steps;
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
    copy_file_overwrite(out_dir / "REFERENCE_SAVE", mc_dir / "REFERENCE_SAVE");
    if (fs::exists(out_dir / "BASE_SAVE")) {
        copy_file_overwrite(out_dir / "BASE_SAVE", mc_dir / "BASE_SAVE");
    }
    copy_file_overwrite(out_dir / "meta.json",mc_dir / "meta.json");

    std::ofstream info(mc_dir / "info.txt");
    if (info) {
        info << "task_id = " << task_id << "\n";
        info << "E_final = " << std::setprecision(12) << E_final << "\n";
    }

    std::cout << "[MC] accepted task " << task_id
              << ", archived to " << mc_dir << "\n";
}


Structure make_relaxed_seed_for_trial(const fs::path& root,
                                      const Structure& current_ref,
                                      const Structure& trial_ref)
{
    Structure trial_init = trial_ref;
    fs::path contcar_path = root / "CONTCAR";
    fs::path neigh_path = root / "intsite_metal_neighbors.dat";

    if (!fs::exists(contcar_path) || !fs::exists(neigh_path)) {
        return trial_init;
    }

    if (!trial_init.seedTrialCoordinatesFromContcar(current_ref,
                                                    contcar_path.string().c_str(),
                                                    neigh_path.string().c_str())) {
        std::cerr << "[WARN] failed to seed trial POSCAR from current CONTCAR; "
                  << "falling back to reference coordinates.\n";
        trial_init = trial_ref;
    }
    return trial_init;
}

// Metropolis acceptance
bool Accept(double Eold, double Enew, double temp,
            std::mt19937_64& rng,
            double hastings_ratio = 1.0)
{
    if (hastings_ratio <= 0.0) return false;
    double logp = -(Enew - Eold)/(temp * k_Boltzmann) + std::log(hastings_ratio);
    if (logp >= 0.0) return true;
    double p = std::exp(logp);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double u = dist(rng);
    return (u < p);
}

struct MoveRecord {
    std::string type = "none";
    int a = -1;
    int b = -1;
    int forward_choices = 0;
    int reverse_choices = 0;
    double hastings_ratio = 1.0;
};

json sparse_descriptor_to_json(const SparseDescriptor& d)
{
    json j = json::object();
    for (const auto& kv : d.values) {
        j[std::to_string(kv.first)] = kv.second;
    }
    return j;
}

SparseDescriptor sparse_descriptor_from_json(const json& j)
{
    SparseDescriptor d;
    if (!j.is_object()) return d;
    for (auto it = j.begin(); it != j.end(); ++it) {
        try {
            int idx = std::stoi(it.key());
            double value = it.value().get<double>();
            if (std::isfinite(value) && std::fabs(value) > 0.0) {
                d.values[idx] = value;
            }
        } catch (...) {
        }
    }
    return d;
}

bool pick_metal_swap_sites(const Structure& struc, int& a, int& b)
{
    a = -1;
    b = -1;
    if (struc.num_metallic_atoms < 2) return false;

    for (int attempt = 0; attempt < 1000; ++attempt) {
        int i = std::rand() % struc.num_metallic_atoms;
        int j = std::rand() % struc.num_metallic_atoms;
        if (i != j && struc.atomtype[i] != struc.atomtype[j]) {
            a = i;
            b = j;
            return true;
        }
    }

    for (int i = 0; i < struc.num_metallic_atoms; ++i) {
        for (int j = i + 1; j < struc.num_metallic_atoms; ++j) {
            if (struc.atomtype[i] != struc.atomtype[j]) {
                a = i;
                b = j;
                return true;
            }
        }
    }
    return false;
}

bool pick_interstitial_swap_sites(const Structure& struc, int& a, int& b)
{
    a = -1;
    b = -1;
    if (struc.num_interstitial < 2) return false;

    std::vector<int> occupied;
    occupied.reserve(struc.num_interstitial);
    for (int i = 0; i < struc.num_interstitial; ++i) {
        if (struc.interstitial_postype[i] != -1) occupied.push_back(i);
    }
    if (occupied.empty()) return false;

    for (int attempt = 0; attempt < 1000; ++attempt) {
        int i = occupied[std::rand() % occupied.size()];
        int j = std::rand() % struc.num_interstitial;
        if (i != j && struc.interstitial_postype[i] != struc.interstitial_postype[j]) {
            a = i;
            b = j;
            return true;
        }
    }

    for (int i : occupied) {
        for (int j = 0; j < struc.num_interstitial; ++j) {
            if (i != j && struc.interstitial_postype[i] != struc.interstitial_postype[j]) {
                a = i;
                b = j;
                return true;
            }
        }
    }
    return false;
}

// Apply one MC proposal move to a Structure.  This is shared by the original
// fast+slow search mode and the finiteT direct-to-slow mode.
MoveRecord apply_random_mc_move(Structure& struc,
                                const Args& cfg,
                                const fs::path& root,
                                int SUMP)
{
    MoveRecord rec;
    const int P1 = cfg.p_swap_metal;
    const int P2 = cfg.p_swap_inter;
    const int P3 = cfg.p_hop_inter;
    const int P4 = cfg.p_cluster_inter;
    const int P5 = cfg.p_exch_metal;

    int r = std::rand() % SUMP;
    if (r < P1) {
        // swapMeta
        int a = -1, b = -1;
        if (!pick_metal_swap_sites(struc, a, b)) {
            std::cerr << "[WARN] metal swap skipped because no unlike metal pair is available.\n";
            return rec;
        }
        struc.swapMetal(a,b);
        rec = {"swap_metal", a, b};
    } else if ((r -= P1) < P2) {
        // swapInterstitial: pick one occupied site and one site with a different
        // occupation state/species.  This preserves the total numbers of B/O.
        int a = -1, b = -1;
        if (!pick_interstitial_swap_sites(struc, a, b)) {
            std::cerr << "[WARN] interstitial swap skipped because no occupied differing site pair is available.\n";
            return rec;
        }
        struc.swapInterstitial(a,b);
        rec = {"swap_interstitial", a, b};
    } else if ((r -= P2) < P3) {
        fs::path hop_path = root / "intsite_hop_neighbors.dat";
        if (!struc.readIntsiteHopNeighborMap(hop_path.string().c_str())) {
            std::cerr << "[WARN] local interstitial hop skipped because hop neighbor map could not be read.\n";
            return rec;
        }
        int a = -1, b = -1, forward = 0, reverse = 0;
        int ret = struc.localHopInterstitialRandom(a, b, forward, reverse);
        if (ret != 1 || forward <= 0 || reverse <= 0) {
            std::cerr << "[WARN] local interstitial hop skipped because no valid hop edge was available.\n";
            return rec;
        }
        rec = {"hop_interstitial", a, b, forward, reverse, (double)forward / (double)reverse};
    } else if ((r -= P3) < P4) {
        // clusterInterstitialSwap: move an interstitial occupation and shuffle
        // the metal atom types in the union of the two neighboring cages.
        fs::path neigh_path = root / "intsite_metal_neighbors.dat";
        if (!struc.readIntsiteMetalNeighborMap(neigh_path.string().c_str())) {
            std::cerr << "[WARN] cluster interstitial move skipped because neighbor map could not be read.\n";
            return rec;
        }
        int a = -1, b = -1;
        if (!pick_interstitial_swap_sites(struc, a, b)) {
            std::cerr << "[WARN] cluster interstitial move skipped because no occupied differing site pair is available.\n";
            return rec;
        }
        struc.clusterSwapInterstitial(a,b);
        rec = {"cluster_swap_interstitial", a, b};
    } else if ((r -= P4) < P5) {
        // TODO: exchangeMetal
        // struc.exchangeMetal(...);
    } else {
        // TODO: exchangeInterstitial
        // struc.exchangeInterstitial(...);
    }
    return rec;
}

struct CandidateTrial {
    Structure trial_ref;
    Structure trial_init;
    MoveRecord move;
    TrialScore prefast_score;
    int generated_index = 0;
};

bool generate_one_trial_from_current(const Structure& current_ref,
                                     const Args& cfg,
                                     const fs::path& root,
                                     int SUMP,
                                     CandidateTrial& out)
{
    Structure trial_ref = current_ref;
    MoveRecord move;
    for (int attempt = 0; attempt < 100; ++attempt) {
        trial_ref = current_ref;
        move = apply_random_mc_move(trial_ref, cfg, root, SUMP);
        if (move.type != "none") break;
    }
    if (move.type == "none") return false;

    out.trial_ref = trial_ref;
    out.trial_init = make_relaxed_seed_for_trial(root, current_ref, trial_ref);
    out.move = move;
    return true;
}

// Generate one candidate for a given fast slot: uses current SAVE.
bool generate_candidate_for_slot(int slot,
                                 const Args& cfg,
                                 const fs::path& root,
                                 Structure& struc,
                                 int SUMP,
                                 PrefastModel* prefast,
                                 double current_E,
                                 int prefast_learning_steps)
{
    // 1) Load current discrete reference state from SAVE.
    struc.readstruc("SAVE");
    Structure current_ref = struc;

    // 2) Generate either one ordinary trial, or a small prefast-ranked batch.
    bool prefast_enabled = (prefast && prefast->enabled());
    bool prefast_warmup_active =
        prefast_enabled &&
        cfg.prefast.warmup_steps > 0 &&
        prefast_learning_steps < cfg.prefast.warmup_steps;
    int n_candidates = (prefast_enabled && !prefast_warmup_active)
                     ? std::max(1, cfg.prefast.candidates_per_slot)
                     : 1;
    std::vector<CandidateTrial> candidates;
    candidates.reserve(n_candidates);
    for (int ci = 0; ci < n_candidates; ++ci) {
        CandidateTrial cand;
        cand.generated_index = ci;
        if (!generate_one_trial_from_current(current_ref, cfg, root, SUMP, cand)) {
            continue;
        }
        if (prefast_enabled) {
            cand.prefast_score = prefast->score(current_ref, cand.trial_ref);
        }
        candidates.push_back(cand);
    }

    if (candidates.empty()) {
        std::cerr << "[WARN] no valid MC move could be generated for fast slot "
                  << slot << ".\n";
        return false;
    }

    int best = 0;
    if (prefast_enabled && !prefast_warmup_active) {
        for (int i = 1; i < (int)candidates.size(); ++i) {
            if (candidates[i].prefast_score.dE_pred <
                candidates[best].prefast_score.dE_pred) {
                best = i;
            }
        }
    }
    CandidateTrial& chosen = candidates[best];
    struc = chosen.trial_ref;

    // 3) Dump candidate directly into fast/POSCARk, fast/SAVEk
    fs::path fast_dir = root / "fast";
    fs::create_directories(fast_dir);

    fs::path pos = fast_dir / ("POSCAR" + std::to_string(slot));
    fs::path sav = fast_dir / ("SAVE"   + std::to_string(slot));
    fs::path base_sav = fast_dir / ("BASE_SAVE" + std::to_string(slot));
    fs::path met = fast_dir / ("META"   + std::to_string(slot));

    chosen.trial_init.outputvasp(pos.string().c_str());
    current_ref.outputsave(base_sav.string().c_str());
    chosen.trial_ref.outputsave(sav.string().c_str());
    json meta;
    meta["source"] = "search_fast_screen";
    meta["mode"] = "search";
    meta["source_slot"] = slot;
    meta["move_type"] = chosen.move.type;
    meta["move_site_a"] = chosen.move.a;
    meta["move_site_b"] = chosen.move.b;
    meta["move_forward_choices"] = chosen.move.forward_choices;
    meta["move_reverse_choices"] = chosen.move.reverse_choices;
    meta["hastings_ratio"] = chosen.move.hastings_ratio;
    if (prefast_enabled) {
        meta["prefast_enabled"] = true;
        meta["prefast_basis"] = cfg.prefast.basis;
        meta["prefast_requested_candidates_per_slot"] = cfg.prefast.candidates_per_slot;
        meta["prefast_candidates_per_slot"] = n_candidates;
        meta["prefast_candidates_generated"] = (int)candidates.size();
        meta["prefast_warmup_steps"] = cfg.prefast.warmup_steps;
        meta["prefast_learning_step_at_proposal"] = prefast_learning_steps;
        meta["prefast_warmup_active"] = prefast_warmup_active;
        meta["prefast_ranking_active"] = !prefast_warmup_active && n_candidates > 1;
        meta["prefast_selected_rank"] = 1;
        meta["prefast_selected_generated_index"] = chosen.generated_index;
        meta["prefast_dE_pred"] = chosen.prefast_score.dE_pred;
        meta["prefast_score"] = chosen.prefast_score.dE_pred;
        meta["prefast_norm_dD"] = chosen.prefast_score.norm_dD;
        meta["prefast_base_energy"] = current_E;
        meta["prefast_delta"] = sparse_descriptor_to_json(chosen.prefast_score.delta);
    } else {
        meta["prefast_enabled"] = false;
    }
    meta["stamp"] = "created_by_cpp_master";
    {
        std::ofstream ofs(met);
        ofs << std::setw(2) << meta << "\n";
    }

    // 4) Trigger fast worker by touching .go_k
    fs::path gof = fast_dir / (".go_" + std::to_string(slot));
    std::ofstream ofs(gof); // touch
    return true;
}

// Submit the initial unmodified structure directly to the slow-worker queue.
// This creates a waiting_pool task without going through any fast worker.
// The first report from this task will be accepted unconditionally as the
// initial current state by process_report_file() when have_state == false.
std::string submit_initial_relax_task(const fs::path& root, Structure& struc)
{
    fs::path pool = root / "waiting_pool";
    fs::create_directories(pool);

    auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    std::string task_id = "initial_relax_" + std::to_string(now_ns);

    fs::path tmpd  = pool / (".tmp_" + task_id);
    fs::path final = pool / task_id;

    fs::create_directories(tmpd);

    struc.outputvasp((tmpd / "POSCAR").string().c_str());
    struc.outputsave((tmpd / "SAVE").string().c_str());

    json meta;
    meta["task_id"] = task_id;
    meta["source"] = "initial_relax";
    // Sentinel value: not a physical energy, only used so slow_worker picks
    // the initial relaxation before any ordinary screened candidates.
    meta["energy_screen"] = -1.0e99;
    meta["note"] = "sentinel energy_screen for initial relaxation; not a physical energy";
    meta["stamp"] = "created_by_cpp_master";

    {
        std::ofstream ofs(tmpd / "meta.json");
        ofs << std::setw(2) << meta << "\n";
    }

    std::error_code ec;
    fs::rename(tmpd, final, ec);
    if (ec) {
        std::cerr << "[ERROR] failed to submit initial relax task: "
                  << ec.message() << "\n";
        return "";
    }

    std::cout << "[INIT] submitted initial structure to slow-worker queue: "
              << task_id << "\n";
    return task_id;
}



// Submit one finite-temperature MC proposal directly to the slow-worker queue.
// This bypasses fast_worker and waiting-pool screening.  The master will submit
// a new task only after the previous one has been processed, so the Markov
// chain remains sequential.
std::string submit_finiteT_trial_task(const fs::path& root,
                                      const Args& cfg,
                                      Structure& struc,
                                      int SUMP,
                                      int next_step)
{
    // Load the current accepted discrete reference state.
    if (!struc.readstruc("SAVE")) {
        std::cerr << "[finiteT] failed to read current SAVE before trial submission.\n";
        return "";
    }
    Structure current_ref = struc;

    MoveRecord move;
    for (int attempt = 0; attempt < 100; ++attempt) {
        struc = current_ref;
        move = apply_random_mc_move(struc, cfg, root, SUMP);
        if (move.type != "none") break;
    }
    if (move.type == "none") {
        std::cerr << "[finiteT] no valid MC move was available for trial submission.\n";
        return "";
    }
    Structure trial_init = make_relaxed_seed_for_trial(root, current_ref, struc);

    fs::path pool = root / "waiting_pool";
    fs::create_directories(pool);

    auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    std::string task_id = "finiteT_step_" + std::to_string(next_step)
                        + "_" + std::to_string(now_ns);

    fs::path tmpd  = pool / (".tmp_" + task_id);
    fs::path final = pool / task_id;
    fs::create_directories(tmpd);

    trial_init.outputvasp((tmpd / "POSCAR").string().c_str());
    struc.outputsave((tmpd / "SAVE").string().c_str());

    json meta;
    meta["task_id"] = task_id;
    meta["source"] = "finiteT_direct";
    meta["mode"] = "finiteT";
    meta["mc_step_proposal"] = next_step;
    meta["move_type"] = move.type;
    meta["move_site_a"] = move.a;
    meta["move_site_b"] = move.b;
    meta["move_forward_choices"] = move.forward_choices;
    meta["move_reverse_choices"] = move.reverse_choices;
    meta["hastings_ratio"] = move.hastings_ratio;
    // There is only one pending finiteT task, so this is not used for screening.
    meta["energy_screen"] = 0.0;
    meta["stamp"] = "created_by_cpp_master";

    {
        std::ofstream ofs(tmpd / "meta.json");
        ofs << std::setw(2) << meta << "\n";
    }

    std::error_code ec;
    fs::rename(tmpd, final, ec);
    if (ec) {
        std::cerr << "[finiteT] failed to submit trial task: "
                  << ec.message() << "\n";
        return "";
    }

    std::cout << "[finiteT] submitted trial step " << next_step
              << " to slow-worker queue: " << task_id << "\n";
    return task_id;
}

/* ---------- process a single slow report ---------- */
bool process_report_file(const fs::path& root,
                         const fs::path& rep_path,
                         double& current_E,
                         bool& have_state,
                         int& mc_steps,
                         int& accept_count,
                         const std::string& mode,
                         double interstitial_site_cutoff,
                         double temp,
                         std::mt19937_64& rng,
                         std::ofstream& log,
                         PrefastModel* prefast = nullptr,
                         bool* accepted_state_changed = nullptr)
{
    if (accepted_state_changed) *accepted_state_changed = false;

    // Parse JSON
    json j;
    {
        std::ifstream ifs(rep_path);
        if (!ifs) {
            std::cerr << "[WARN] cannot open report " << rep_path << "\n";
            cleanup_processed_task(root, "", rep_path, false);
            return false;
        }
        try {
            ifs >> j;
        } catch (...) {
            std::cerr << "[WARN] bad JSON in " << rep_path << "\n";
            cleanup_processed_task(root, "", rep_path, false);
            return false;
        }
    }

    std::string status = j.value("status", "");
    if (status == "error") {
        std::cerr << "[slow] ERROR in report " << rep_path.filename()
                  << ": " << j.value("error", std::string("<no_msg>")) << "\n";
        std::string task_id = j.value("task_id", rep_path.stem().string());
        cleanup_processed_task(root, task_id, rep_path, true);
        return false;
    }

    std::string task_id = j.value("task_id", rep_path.stem().string());
    double E_final = j.value("energy_final", std::numeric_limits<double>::infinity());

    if (!std::isfinite(E_final)) {
        std::cerr << "[WARN] invalid energy_final in report " << rep_path << "\n";
        cleanup_processed_task(root, task_id, rep_path, true);
        return false;
    }

    fs::path out_dir = root / "refine_outbox" / task_id;
    fs::path contcar_path = out_dir / "CONTCAR";
    fs::path save_path = out_dir / "SAVE";
    fs::path neigh_path = root / "intsite_metal_neighbors.dat";
    fs::path meta_path = out_dir / "meta.json";

    if (!have_state) {
        // First-ever refined structure: treat as initial state.
        Structure initial_ref;
        int n_reassigned = 0;
        if (!initial_ref.readstruc(save_path.string().c_str()) ||
            !initial_ref.reconcileInterstitialSitesFromContcar(
                contcar_path.string().c_str(),
                neigh_path.string().c_str(),
                1,
                interstitial_site_cutoff,
                contcar_path.string().c_str(),
                n_reassigned)) {
            std::cerr << "[WARN] failed to reconcile initial interstitial sites for "
                      << task_id << "\n";
            cleanup_processed_task(root, task_id, rep_path, false);
            return false;
        }
        if (n_reassigned > 0) {
            initial_ref.outputsave(save_path.string().c_str());
            log << "INITIAL_REASSIGN task_id=" << task_id
                << " n_reassigned=" << n_reassigned << "\n";
        }

        have_state = true;
        current_E  = E_final;
        accept_count = 0;
        mc_steps = 0;

        // Keep SAVE as the discrete reference state.  CONTCAR carries the
        // relaxed physical coordinates used to seed later trial POSCAR files.
        copy_file_overwrite(out_dir / "SAVE",    root / "SAVE");
        copy_file_overwrite(out_dir / "CONTCAR", root / "CONTCAR");

        log << "INITIAL_STATE task_id=" << task_id
            << " E = " << std::setprecision(12) << E_final << "\n";
        log.flush();

        archive_mc_accept(root, task_id, E_final);
        if (accepted_state_changed) *accepted_state_changed = true;
        cleanup_processed_task(root, task_id, rep_path, true);
        return true;
    }

    bool allow_interstitial_reassign = (mode == "search");
    Structure trial_ref;
    int n_reassigned = 0;
    bool reconcile_ok =
        trial_ref.readstruc(save_path.string().c_str()) &&
        trial_ref.reconcileInterstitialSitesFromContcar(
            contcar_path.string().c_str(),
            neigh_path.string().c_str(),
            allow_interstitial_reassign ? 1 : 0,
            interstitial_site_cutoff,
            allow_interstitial_reassign ? contcar_path.string().c_str() : "",
            n_reassigned);

    if (!reconcile_ok) {
        ++mc_steps;
        log << "STEP " << mc_steps
            << " proposal task_id=" << task_id
            << " E_new=" << std::setprecision(12) << E_final
            << " E_old=" << current_E
            << " -> REJECT_INTERSTITIAL_HOP\n";
        log.flush();
        cleanup_processed_task(root, task_id, rep_path, true);
        return true;
    }

    if (allow_interstitial_reassign && n_reassigned > 0) {
        trial_ref.outputsave(save_path.string().c_str());
        log << "SEARCH_REASSIGN task_id=" << task_id
            << " n_reassigned=" << n_reassigned << "\n";
        log.flush();
    }

    double hastings_ratio = 1.0;
    std::string move_type = "<unknown>";
    bool meta_prefast_enabled = false;
    double prefast_dE_pred = 0.0;
    double prefast_base_energy = current_E;
    SparseDescriptor prefast_delta;
    SparseDescriptor prefast_learning_delta;
    std::string prefast_learning_delta_source = "none";
    bool have_prefast_learning_delta = false;
    bool computed_prefast_learning_delta = false;
    bool have_base_ref = false;
    bool final_same_as_base = false;
    try {
        std::ifstream meta_ifs(meta_path);
        if (meta_ifs) {
            json meta;
            meta_ifs >> meta;
            hastings_ratio = meta.value("hastings_ratio", 1.0);
            move_type = meta.value("move_type", move_type);
            meta_prefast_enabled = meta.value("prefast_enabled", false);
            prefast_dE_pred = meta.value("prefast_dE_pred", 0.0);
            prefast_base_energy = meta.value("prefast_base_energy", current_E);
            if (meta.contains("prefast_delta")) {
                prefast_delta = sparse_descriptor_from_json(meta["prefast_delta"]);
            }
        }
    } catch (...) {
        hastings_ratio = 1.0;
    }

    Structure prefast_base_ref;
    fs::path prefast_base_save_path = out_dir / "BASE_SAVE";
    if (fs::exists(prefast_base_save_path) &&
        prefast_base_ref.readstruc(prefast_base_save_path.string().c_str())) {
        have_base_ref = true;
        final_same_as_base =
            prefast_base_ref.atomtype == trial_ref.atomtype &&
            prefast_base_ref.interstitial_postype == trial_ref.interstitial_postype;
        if (prefast && prefast->enabled()) {
            TrialScore learning_score = prefast->score(prefast_base_ref, trial_ref);
            prefast_learning_delta = learning_score.delta;
            prefast_learning_delta_source =
                (allow_interstitial_reassign && n_reassigned > 0)
                    ? "final_reassigned"
                    : "final";
            computed_prefast_learning_delta = true;
            have_prefast_learning_delta = !prefast_learning_delta.values.empty();
        }
    } else if (!prefast_delta.values.empty()) {
        prefast_learning_delta = prefast_delta;
        prefast_learning_delta_source = "trial_meta";
        have_prefast_learning_delta = true;
    }

    if (have_base_ref &&
        (final_same_as_base ||
         (computed_prefast_learning_delta && prefast_learning_delta.values.empty()))) {
        log << "DISCARD_NO_CHANGE task_id=" << task_id
            << " move=" << move_type
            << " delta_source=" << prefast_learning_delta_source
            << " reason="
            << (final_same_as_base ? "final_same_as_base" : "empty_descriptor_delta")
            << " n_reassigned=" << n_reassigned
            << "\n";
        log.flush();
        cleanup_processed_task(root, task_id, rep_path, true);
        return true;
    }

    // Normal MC proposal
    ++mc_steps;
    bool accept = Accept(current_E, E_final, temp, rng, hastings_ratio);
    if (prefast && prefast->enabled() && meta_prefast_enabled && have_prefast_learning_delta) {
        double dE_true_for_learning = E_final - prefast_base_energy;
        auto update_stats = prefast->update(prefast_learning_delta, dE_true_for_learning);
        if (prefast->log_summary()) {
            prefast->append_learning_log(root / "prefast_learning.log",
                                         mc_steps,
                                         task_id,
                                         prefast_dE_pred,
                                         update_stats,
                                         accept,
                                         prefast_learning_delta_source,
                                         n_reassigned);
        }
        if (prefast->log_weight_updates()) {
            prefast->append_weight_update_log(root / "prefast_weight_updates.log",
                                              mc_steps,
                                              task_id,
                                              update_stats);
        }
        if (prefast->log_weight_snapshots()) {
            prefast->append_weight_snapshot_log(root / "prefast_weights.log",
                                                mc_steps,
                                                task_id);
        }
    }
    log << "STEP " << mc_steps
        << " proposal task_id=" << task_id
        << " move=" << move_type
        << " E_new=" << std::setprecision(12) << E_final
        << " E_old=" << current_E
        << " hastings_ratio=" << hastings_ratio
        << " -> "   << (accept ? "ACCEPT" : "REJECT")
        << "\n";
    log.flush();

    if (accept) {
        ++accept_count;
        current_E = E_final;
        // Accept the trial's discrete reference state.  Do not overwrite SAVE
        // with relaxed coordinates; keep CONTCAR as the relaxed physical state.
        copy_file_overwrite(out_dir / "SAVE",  root / "SAVE");
        copy_file_overwrite(out_dir / "CONTCAR", root /"CONTCAR");
        archive_mc_accept(root, task_id, E_final);
        if (accepted_state_changed) *accepted_state_changed = true;
    }

    cleanup_processed_task(root, task_id, rep_path, true);
    return true;
}

// Wait until the initial slow-worker refinement finishes and initializes SAVE/CONTCAR.
void wait_for_initial_state(const fs::path& root,
                            double& current_E,
                            bool& have_state,
                            int& mc_steps,
                            int& accept_count,
                            const std::string& mode,
                            double interstitial_site_cutoff,
                            double temp,
                            std::mt19937_64& rng,
                            std::ofstream& log)
{
    fs::path reports_dir = root / "reports";
    std::cout << "[INIT] waiting for initial relaxed state...\n";

    while (!have_state) {
        bool processed_any = false;
        for (auto& entry : fs::directory_iterator(reports_dir)) {
            if (!entry.is_regular_file()) continue;
            fs::path rep_path = entry.path();
            if (rep_path.extension() != ".json") continue;

            bool ok = process_report_file(
                root, rep_path,
                current_E, have_state,
                mc_steps, accept_count,
                mode, interstitial_site_cutoff,
                temp, rng, log
            );
            if (ok) processed_any = true;
            if (have_state) break;
        }

        if (!processed_any) {
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
        }
    }

    std::cout << "[INIT] initial relaxed state is ready. E_current = "
              << std::setprecision(12) << current_E << "\n";
}


/* ---------- main ---------- */

int main(int argc, char** argv){
    std::ios::sync_with_stdio(false);
    std::srand((unsigned)std::time(nullptr));
    std::mt19937_64 rng((unsigned)std::time(nullptr));

    Args cfg = parse_args(argc, argv);
    const int SUMP = cfg.p_swap_metal + cfg.p_swap_inter +
                     cfg.p_hop_inter +
                     cfg.p_cluster_inter +
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

    // Load the reference structure used for MC moves.  In resume mode, the
    // accepted seed SAVE is the authoritative reference state.
    Structure struc;
    fs::path resolved_resume_dir;
    std::string seed_struc;
    if (cfg.continue_run) {
        resolved_resume_dir = latest_mcprocess_state_dir(ROOT);
        if (!looks_like_state_dir(resolved_resume_dir)) {
            std::cerr << "[ERROR] --continue could not find a valid latest state under "
                      << (ROOT / "mcprocess") << "\n";
            return 1;
        }
        seed_struc = (resolved_resume_dir / "SAVE").string();
    } else if (!cfg.resume_state_dir.empty()) {
        resolved_resume_dir = resolve_resume_state_dir(cfg.resume_state_dir);
        seed_struc = (resolved_resume_dir / "SAVE").string();
    } else {
        seed_struc = cfg.input_struc;
    }
    if (!struc.readstruc(seed_struc.c_str())) {
        std::cerr << "[ERROR] failed to read structure seed: "
                  << seed_struc << "\n";
        return 1;
    }
    struc.buildIntsiteMetalNeighborMap(cfg.intsite_neighbor_cutoff);
    struc.outputIntsiteMetalNeighborMap("intsite_metal_neighbors.dat");
    struc.buildIntsiteHopNeighborMap(cfg.intsite_hop_cutoff);
    struc.outputIntsiteHopNeighborMap("intsite_hop_neighbors.dat");
    struc.outputsave("SAVE");

    PrefastModel prefast;
    if (cfg.prefast.enabled) {
        prefast.build(struc, cfg.prefast);
        if (prefast.log_summary()) {
            prefast.print_startup_log(std::cout);
            prefast.write_basis_log(ROOT / "prefast_basis.log");
        } else if (prefast.enabled()) {
            std::cout << "[prefast] enabled; diagnostics = off\n";
        }
    }

    std::ofstream log;
    if (cfg.continue_run) {
        log.open("mc.log", std::ios::app);
        log << "\n";
    } else {
        log.open("mc.log");
    }
    log << "# MC mode=" << cfg.mode
        << " fast=" << cfg.workers
        << " steps=" << cfg.steps
        << " step_counter=" << (cfg.mode == "search" ? "fast_count" : "mc_steps")
        << " temp=" << cfg.temp
        << " p_swap_metal=" << cfg.p_swap_metal
        << " p_swap_inter=" << cfg.p_swap_inter
        << " p_hop_inter=" << cfg.p_hop_inter
        << " p_cluster_inter=" << cfg.p_cluster_inter
        << " p_exch_metal=" << cfg.p_exch_metal
        << " p_exch_inter=" << cfg.p_exch_inter
        << " intsite_neighbor_cutoff=" << cfg.intsite_neighbor_cutoff
        << " intsite_hop_cutoff=" << cfg.intsite_hop_cutoff
        << " interstitial_site_cutoff=" << cfg.interstitial_site_cutoff
        << " resume_state=" << (cfg.resume_state_dir.empty() ? "<none>" : cfg.resume_state_dir)
        << " continue=" << (cfg.continue_run ? "true" : "false")
        << " prefast=" << (cfg.prefast.enabled ? "on" : "off")
        << " prefast_candidates_per_slot=" << cfg.prefast.candidates_per_slot
        << " prefast_warmup_steps=" << cfg.prefast.warmup_steps
        << " prefast_basis=" << cfg.prefast.basis
        << " prefast_nshells=" << cfg.prefast.nshells
        << " prefast_peak_scan_cutoff=" << cfg.prefast.peak_scan_cutoff
        << " prefast_peak_tol=" << cfg.prefast.peak_tol
        << " prefast_sigma_small=" << cfg.prefast.sigma_small
        << " prefast_sigma_large=" << cfg.prefast.sigma_large
        << " prefast_cutoff_margin=" << cfg.prefast.cutoff_margin
        << " prefast_learning_rate=" << cfg.prefast.learning_rate
        << " prefast_lms_epsilon=" << cfg.prefast.lms_epsilon
        << " prefast_weight_decay=" << cfg.prefast.weight_decay
        << " prefast_diagnostics=" << cfg.prefast.diagnostics
        << "\n";

    bool have_state = false;
    double current_E = 0.0;
    int mc_steps = 0;
    int accept_count = 0;

    if (cfg.continue_run) {
        if (!initialize_from_resume_state(ROOT, resolved_resume_dir,
                                          current_E, have_state,
                                          mc_steps, accept_count, log)) {
            std::cerr << "[ERROR] failed to initialize from --continue state.\n";
            return 1;
        }
        if (auto last_step = read_last_logged_mc_step(ROOT / "mc.log")) {
            mc_steps = *last_step;
        }
        log << "CONTINUE_RUN from=" << resolved_resume_dir.string()
            << " mc_steps=" << mc_steps
            << " fast_screened=" << read_counter_value(ROOT, "fast_count")
            << " mc_count=" << read_counter_value(ROOT, "mc_count")
            << "\n";
        log.flush();
        std::cout << "[continue] continuing from " << resolved_resume_dir
                  << " mc_steps=" << mc_steps
                  << " fast_screened=" << read_counter_value(ROOT, "fast_count")
                  << "\n";
    } else if (!cfg.resume_state_dir.empty()) {
        if (!initialize_from_resume_state(ROOT, resolved_resume_dir,
                                          current_E, have_state,
                                          mc_steps, accept_count, log)) {
            std::cerr << "[ERROR] failed to initialize from --resume-state.\n";
            return 1;
        }
    } else {
        // First, refine the initial structure directly with the slow-worker path.
        // This establishes a well-defined relaxed current state and E_current
        // before any MC candidates are generated.
        std::string init_task_id = submit_initial_relax_task(ROOT, struc);
        if (init_task_id.empty()) {
            std::cerr << "[ERROR] failed to submit initial relaxation task.\n";
            return 1;
        }
        wait_for_initial_state(ROOT, current_E, have_state, mc_steps,
                               accept_count, cfg.mode, cfg.interstitial_site_cutoff,
                               cfg.temp, rng, log);
    }

    // Main loop.
    //   search  mode: stop when fast workers have screened cfg.steps trials.
    //   finiteT mode: strict sequential chain; stop after cfg.steps processed
    //                 MC proposals.
    while (!reached_step_budget(ROOT, cfg, mc_steps)) {
        // 1) Poll reports directory for new slow results first.  In finiteT mode
        // this usually consumes the single pending trial and advances the chain.
        bool processed_any = false;
        fs::path reports_dir = ROOT / "reports";
        for (auto& entry : fs::directory_iterator(reports_dir)) {
            if (!entry.is_regular_file()) continue;
            fs::path rep_path = entry.path();
            if (rep_path.extension() != ".json") continue;

            bool accepted_state_changed = false;
            bool ok = process_report_file(
                ROOT, rep_path,
                current_E, have_state,
                mc_steps, accept_count,
                cfg.mode, cfg.interstitial_site_cutoff,
                cfg.temp, rng, log,
                &prefast,
                &accepted_state_changed
            );
            if (ok) processed_any = true;
            if (reached_step_budget(ROOT, cfg, mc_steps)) break;
        }

        if (reached_step_budget(ROOT, cfg, mc_steps)) break;

        if (cfg.mode == "finiteT") {
            // Only one outstanding proposal is allowed.  This keeps the chain:
            // current -> one proposal -> slow refine -> accept/reject -> next.
            if (!has_pending_slow_task(ROOT)) {
                std::string tid = submit_finiteT_trial_task(
                    ROOT, cfg, struc, SUMP, mc_steps + 1
                );
                if (!tid.empty()) processed_any = true;
            }
        } else {
            // Search mode: keep fast-worker slots filled.
            for (int k = 1; k <= cfg.workers; ++k) {
                int screened = read_counter_value(ROOT, "fast_count");
                int busy_fast = count_busy_fast_slots(ROOT, cfg.workers);
                if (screened + busy_fast >= cfg.steps) break;

                fs::path gof = ROOT / "fast" / (".go_" + std::to_string(k));
                if (fs::exists(gof)) continue;  // slot is busy
                if (generate_candidate_for_slot(k, cfg, ROOT, struc, SUMP, &prefast, current_E, mc_steps)) {
                    processed_any = true;
                    cout << "generating candidate for worker #" << k << " successfully" << endl;
                }
            }
        }

        if (!processed_any) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }

    int fast_screened = read_counter_value(ROOT, "fast_count");
    log << "# Finished. MC steps = " << mc_steps
        << ", fast screened = " << fast_screened
        << ", accepted = " << accept_count << "\n";
    log.close();

    std::cout << "MC finished: steps=" << mc_steps
              << " fast_screened=" << fast_screened
              << " accepted=" << accept_count << "\n";
    return 0;
}
