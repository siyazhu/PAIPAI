#include "findinter.h"
#include "element.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>

using namespace std;

namespace {

struct Vec3 {
    Real x = 0.0;
    Real y = 0.0;
    Real z = 0.0;
};

struct MetalStructure {
    string title;
    Vec3 a, b, c;
    vector<int> types;
    vector<int> counts;
    vector<int> atom_type;
    vector<Vec3> cart;
};

struct Candidate {
    int ix = 0;
    int iy = 0;
    int iz = 0;
    Vec3 frac;
    Vec3 cart;
    Real clearance = 0.0;
    Real nearest_distance = 0.0;
};

constexpr Real DEAD_CLEARANCE = -1.0e30;

Vec3 add(const Vec3& u, const Vec3& v) {
    return {u.x + v.x, u.y + v.y, u.z + v.z};
}

Vec3 sub(const Vec3& u, const Vec3& v) {
    return {u.x - v.x, u.y - v.y, u.z - v.z};
}

Vec3 scale(const Vec3& u, Real s) {
    return {u.x * s, u.y * s, u.z * s};
}

Real dot(const Vec3& u, const Vec3& v) {
    return u.x*v.x + u.y*v.y + u.z*v.z;
}

Vec3 frac_to_cart(const Vec3& f, const Vec3& a, const Vec3& b, const Vec3& c) {
    return add(add(scale(a, f.x), scale(b, f.y)), scale(c, f.z));
}

bool inverse3x3_columns(const Vec3& a, const Vec3& b, const Vec3& c, Real inv[3][3]) {
    Real m00=a.x, m01=b.x, m02=c.x;
    Real m10=a.y, m11=b.y, m12=c.y;
    Real m20=a.z, m21=b.z, m22=c.z;
    Real det = m00*(m11*m22 - m12*m21)
             - m01*(m10*m22 - m12*m20)
             + m02*(m10*m21 - m11*m20);
    if (fabs(det) < 1e-14) return false;
    Real id = 1.0/det;
    inv[0][0] =  (m11*m22 - m12*m21)*id;
    inv[0][1] = -(m01*m22 - m02*m21)*id;
    inv[0][2] =  (m01*m12 - m02*m11)*id;
    inv[1][0] = -(m10*m22 - m12*m20)*id;
    inv[1][1] =  (m00*m22 - m02*m20)*id;
    inv[1][2] = -(m00*m12 - m02*m10)*id;
    inv[2][0] =  (m10*m21 - m11*m20)*id;
    inv[2][1] = -(m00*m21 - m01*m20)*id;
    inv[2][2] =  (m00*m11 - m01*m10)*id;
    return true;
}

Vec3 cart_to_frac(const Vec3& r, Real inv[3][3]) {
    return {
        inv[0][0]*r.x + inv[0][1]*r.y + inv[0][2]*r.z,
        inv[1][0]*r.x + inv[1][1]*r.y + inv[1][2]*r.z,
        inv[2][0]*r.x + inv[2][1]*r.y + inv[2][2]*r.z
    };
}

Vec3 wrap_frac(Vec3 f) {
    f.x -= floor(f.x);
    f.y -= floor(f.y);
    f.z -= floor(f.z);
    return f;
}

Vec3 minimum_image_delta(const Vec3& delta, const Vec3& a, const Vec3& b, const Vec3& c, Real inv[3][3]) {
    Vec3 f = cart_to_frac(delta, inv);
    f.x -= round(f.x);
    f.y -= round(f.y);
    f.z -= round(f.z);
    return frac_to_cart(f, a, b, c);
}

Real norm(const Vec3& u) {
    return sqrt(dot(u,u));
}

bool read_nonempty_line(ifstream& input, string& line) {
    while (getline(input, line)) {
        bool only_space = true;
        for (char ch : line) {
            if (!isspace(static_cast<unsigned char>(ch))) {
                only_space = false;
                break;
            }
        }
        if (!only_space) return true;
    }
    return false;
}

bool line_all_ints(const string& line) {
    stringstream ss(line);
    string tok;
    bool any = false;
    while (ss >> tok) {
        any = true;
        char* endptr = nullptr;
        strtol(tok.c_str(), &endptr, 10);
        if (endptr == tok.c_str() || *endptr != '\0') return false;
    }
    return any;
}

bool read_metal_poscar(const string& filename, MetalStructure& s) {
    ifstream input(filename);
    if (!input) {
        cerr << "Cannot open POSCAR file: " << filename << endl;
        return false;
    }

    if (!getline(input, s.title)) return false;
    Real scale_factor = 1.0;
    input >> scale_factor;
    input >> s.a.x >> s.a.y >> s.a.z;
    input >> s.b.x >> s.b.y >> s.b.z;
    input >> s.c.x >> s.c.y >> s.c.z;
    s.a = scale(s.a, scale_factor);
    s.b = scale(s.b, scale_factor);
    s.c = scale(s.c, scale_factor);
    string line;
    getline(input, line);

    string species_line;
    if (!read_nonempty_line(input, species_line)) return false;
    string counts_line;
    vector<string> species;
    if (line_all_ints(species_line)) {
        counts_line = species_line;
    } else {
        stringstream ss(species_line);
        string name;
        while (ss >> name) species.push_back(name);
        if (!read_nonempty_line(input, counts_line)) return false;
    }

    stringstream cs(counts_line);
    int count = 0;
    while (cs >> count) s.counts.push_back(count);
    if (s.counts.empty()) {
        cerr << "Failed to parse POSCAR atom counts." << endl;
        return false;
    }

    if (species.empty()) {
        cerr << "VASP4-style POSCAR without species names is not supported by findinter." << endl;
        return false;
    }
    if (species.size() != s.counts.size()) {
        cerr << "Species/count mismatch in POSCAR." << endl;
        return false;
    }
    for (const string& name : species) {
        int z = elementnum(name);
        if (z == 0) {
            cerr << "Unknown element in POSCAR: " << name << endl;
            return false;
        }
        s.types.push_back(z);
    }

    string coordtype;
    if (!read_nonempty_line(input, coordtype)) return false;
    if (!coordtype.empty() && (coordtype[0] == 'S' || coordtype[0] == 's')) {
        if (!read_nonempty_line(input, coordtype)) return false;
    }
    bool cartesian = (!coordtype.empty() &&
                      (coordtype[0] == 'C' || coordtype[0] == 'c' ||
                       coordtype[0] == 'K' || coordtype[0] == 'k'));

    int total = 0;
    for (int n : s.counts) total += n;
    s.cart.reserve(total);
    s.atom_type.reserve(total);
    for (int kind = 0; kind < (int)s.counts.size(); kind++) {
        for (int j = 0; j < s.counts[kind]; j++) {
            if (!read_nonempty_line(input, line)) {
                cerr << "Unexpected EOF while reading POSCAR coordinates." << endl;
                return false;
            }
            stringstream ls(line);
            Real x,y,z;
            if (!(ls >> x >> y >> z)) {
                cerr << "Cannot parse coordinate line: " << line << endl;
                return false;
            }
            Vec3 r;
            if (cartesian) r = {x * scale_factor, y * scale_factor, z * scale_factor};
            else r = frac_to_cart({x,y,z}, s.a, s.b, s.c);
            s.cart.push_back(r);
            s.atom_type.push_back(kind);
        }
    }
    return true;
}

Real builtin_radius(int atomic_number) {
    if (atomic_number <= 0 || atomic_number > 86) return 1.0;
    return ELEMENT_RADII[atomic_number - 1];
}

bool write_default_radii(const string& filename,
                         const MetalStructure& structure,
                         const vector<int>& inter_types) {
    Real min_metal = numeric_limits<Real>::infinity();
    Real min_inter = numeric_limits<Real>::infinity();
    for (int z : structure.types) min_metal = min(min_metal, builtin_radius(z));
    for (int z : inter_types) min_inter = min(min_inter, builtin_radius(z));

    ofstream out(filename);
    if (!out) {
        cerr << "Cannot write radii file: " << filename << endl;
        return false;
    }
    out << "# PAIPAI findinter hard-sphere radii in Angstrom\n";
    out << "# Edit this file if the default radii are not appropriate, then rerun the same findinter command.\n";
    out << "# Metal elements are initialized to the minimum built-in radius among the metal species.\n";
    out << "# Interstitial is initialized to the minimum built-in radius among requested interstitial species.\n";
    out << "# Symbol Radius\n";
    for (int z : structure.types) {
        out << PERIODICTABLE[z-1] << " " << setprecision(8) << min_metal << "\n";
    }
    out << "Interstitial " << setprecision(8) << min_inter << "\n";
    return true;
}

bool read_radii(const string& filename, map<string, Real>& radii) {
    ifstream input(filename);
    if (!input) return false;
    string line;
    while (getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        string key;
        Real value;
        stringstream ss(line);
        if (ss >> key >> value) radii[key] = value;
    }
    return true;
}

bool ensure_radii(const FindInterOptions& options,
                  const MetalStructure& structure,
                  map<string, Real>& radii) {
    if (!read_radii(options.radii_file, radii)) {
        if (!write_default_radii(options.radii_file, structure, options.interstitial_types)) {
            return false;
        }

        cout << "Generated " << options.radii_file << " from built-in radii.\n";
        cout << "Please inspect/edit the file, then rerun the same findinter command.\n\n";
        cout << left << setw(16) << "Symbol" << "Radius_A\n";
        cout << left << setw(16) << "------" << "--------\n";
        map<string, Real> generated;
        read_radii(options.radii_file, generated);
        for (int z : structure.types) {
            string key = PERIODICTABLE[z-1];
            cout << left << setw(16) << key << generated[key] << "\n";
        }
        cout << left << setw(16) << "Interstitial" << generated["Interstitial"] << "\n";
        return false;
    }

    for (int z : structure.types) {
        string key = PERIODICTABLE[z-1];
        if (radii.find(key) == radii.end()) {
            cerr << "Missing radius for metal element " << key
                 << " in " << options.radii_file << endl;
            return false;
        }
    }
    if (radii.find("Interstitial") == radii.end()) {
        cerr << "Missing Interstitial radius in " << options.radii_file << endl;
        return false;
    }
    return true;
}

int grid_index(int ix, int iy, int iz, int n) {
    return (ix * n + iy) * n + iz;
}

bool evaluate_grid_point(const MetalStructure& structure,
                         const map<string, Real>& radii,
                         const Vec3& frac,
                         Real inter_radius,
                         Real min_void_factor,
                         Real max_void_factor,
                         Real inv[3][3],
                         Candidate& cand) {
    Vec3 cart = frac_to_cart(frac, structure.a, structure.b, structure.c);
    Real best_clearance = numeric_limits<Real>::infinity();
    Real nearest_distance = numeric_limits<Real>::infinity();

    for (int i = 0; i < (int)structure.cart.size(); i++) {
        int kind = structure.atom_type[i];
        string symbol = PERIODICTABLE[structure.types[kind] - 1];
        Real metal_radius = radii.at(symbol);
        Vec3 delta = minimum_image_delta(sub(cart, structure.cart[i]),
                                         structure.a, structure.b, structure.c, inv);
        Real d = norm(delta);
        Real clearance = d - (metal_radius + inter_radius);
        if (clearance < best_clearance) {
            best_clearance = clearance;
            nearest_distance = d;
        }
    }

    Real nearest_radius_sum = nearest_distance - best_clearance;
    if (nearest_distance < min_void_factor * nearest_radius_sum) return false;
    if (nearest_distance > max_void_factor * nearest_radius_sum) return false;

    cand.frac = frac;
    cand.cart = cart;
    cand.clearance = best_clearance;
    cand.nearest_distance = nearest_distance;
    return true;
}

vector<Candidate> find_grid_maxima(const FindInterOptions& options,
                                   const MetalStructure& structure,
                                   const map<string, Real>& radii) {
    int n = options.grid;
    int total = n*n*n;
    vector<Real> clearance(total, DEAD_CLEARANCE);
    vector<Candidate> grid_candidates(total);
    Real inv[3][3];
    if (!inverse3x3_columns(structure.a, structure.b, structure.c, inv)) {
        cerr << "Cell matrix is singular." << endl;
        return {};
    }

    Real inter_radius = radii.at("Interstitial");
    for (int ix = 0; ix < n; ix++) {
        for (int iy = 0; iy < n; iy++) {
            for (int iz = 0; iz < n; iz++) {
                Vec3 frac = {ix / (Real)n,
                             iy / (Real)n,
                             iz / (Real)n};
                Candidate cand;
                cand.ix = ix; cand.iy = iy; cand.iz = iz;
                int idx = grid_index(ix, iy, iz, n);
                if (evaluate_grid_point(structure, radii, frac, inter_radius,
                                        options.min_void_factor,
                                        options.max_void_factor, inv, cand)) {
                    clearance[idx] = cand.clearance;
                    grid_candidates[idx] = cand;
                }
            }
        }
    }

    vector<Candidate> selected;
    Real merge_distance = options.merge_distance;
    Real exclusion_distance = 2.0 * inter_radius * options.min_void_factor;

    while ((int)selected.size() < options.max_sites) {
        Candidate picked;
        bool found_pick = false;
        Real best_pick_clearance = DEAD_CLEARANCE;

        for (int ix = 0; ix < n; ix++) {
            for (int iy = 0; iy < n; iy++) {
                for (int iz = 0; iz < n; iz++) {
                    int idx = grid_index(ix, iy, iz, n);
                    if (clearance[idx] <= DEAD_CLEARANCE * 0.5) continue;
                    bool is_max = true;
                    for (int dx = -1; dx <= 1 && is_max; dx++) {
                        for (int dy = -1; dy <= 1 && is_max; dy++) {
                            for (int dz = -1; dz <= 1; dz++) {
                                if (dx == 0 && dy == 0 && dz == 0) continue;
                                int jx = (ix + dx + n) % n;
                                int jy = (iy + dy + n) % n;
                                int jz = (iz + dz + n) % n;
                                int jdx = grid_index(jx, jy, jz, n);
                                if (clearance[jdx] > clearance[idx]) {
                                    is_max = false;
                                    break;
                                }
                            }
                        }
                    }
                    if (!is_max) continue;

                    const Candidate& cand = grid_candidates[idx];
                    bool too_close = false;
                    if (merge_distance > 0.0) {
                        for (const Candidate& existing : selected) {
                            Vec3 delta = minimum_image_delta(sub(cand.cart, existing.cart),
                                                             structure.a, structure.b, structure.c, inv);
                            if (norm(delta) < merge_distance) {
                                too_close = true;
                                break;
                            }
                        }
                    }
                    if (too_close) continue;

                    if (!found_pick || cand.clearance > best_pick_clearance) {
                        picked = cand;
                        best_pick_clearance = cand.clearance;
                        found_pick = true;
                    }
                }
            }
        }

        if (!found_pick) break;
        selected.push_back(picked);

        for (int idx = 0; idx < total; idx++) {
            if (clearance[idx] <= DEAD_CLEARANCE * 0.5) continue;
            Vec3 delta = minimum_image_delta(sub(grid_candidates[idx].cart, picked.cart),
                                             structure.a, structure.b, structure.c, inv);
            if (norm(delta) < exclusion_distance) {
                clearance[idx] = DEAD_CLEARANCE;
            }
        }
    }

    return selected;
}

bool write_struc_in(const FindInterOptions& options,
                    const MetalStructure& structure,
                    const vector<Candidate>& sites) {
    int requested_interstitials = 0;
    for (int n : options.interstitial_counts) requested_interstitials += n;
    if ((int)sites.size() < requested_interstitials) {
        cerr << "Found only " << sites.size()
             << " interstitial sites, but requested "
             << requested_interstitials << " occupied interstitial atoms." << endl;
        return false;
    }

    ofstream out(options.output_struc);
    if (!out) {
        cerr << "Cannot write output file: " << options.output_struc << endl;
        return false;
    }

    out << "PAIPAI findinter generated structure\n";
    out << "1.0\n";
    out << structure.a.x << " " << structure.a.y << " " << structure.a.z << "\n";
    out << structure.b.x << " " << structure.b.y << " " << structure.b.z << "\n";
    out << structure.c.x << " " << structure.c.y << " " << structure.c.z << "\n";

    for (int i = 0; i < (int)structure.types.size(); i++) {
        if (i) out << " ";
        out << PERIODICTABLE[structure.types[i]-1];
    }
    out << "\n";
    for (int i = 0; i < (int)structure.counts.size(); i++) {
        if (i) out << " ";
        out << structure.counts[i];
    }
    out << "\n";

    for (int i = 0; i < (int)options.interstitial_types.size(); i++) {
        if (i) out << " ";
        out << PERIODICTABLE[options.interstitial_types[i]-1];
    }
    out << "\n";
    for (int i = 0; i < (int)options.interstitial_counts.size(); i++) {
        if (i) out << " ";
        out << options.interstitial_counts[i];
    }
    out << "\n";

    out << sites.size() << "\n";
    out << "No Shuffle\n";
    out << "Cartesian\n";

    for (int i = 0; i < (int)structure.cart.size(); i++) {
        out << structure.cart[i].x << " " << structure.cart[i].y << " "
            << structure.cart[i].z << " " << structure.atom_type[i] << "\n";
    }

    int occupied = 0;
    for (int species = 0; species < (int)options.interstitial_counts.size(); species++) {
        for (int j = 0; j < options.interstitial_counts[species]; j++) {
            const Candidate& site = sites[occupied++];
            out << site.cart.x << " " << site.cart.y << " " << site.cart.z
                << " " << species << "\n";
        }
    }
    for (int i = occupied; i < (int)sites.size(); i++) {
        out << sites[i].cart.x << " " << sites[i].cart.y << " " << sites[i].cart.z
            << " -1\n";
    }
    return true;
}

bool write_site_poscar(const string& filename,
                       const MetalStructure& structure,
                       const vector<Candidate>& sites) {
    if (filename.empty()) return true;

    ofstream out(filename);
    if (!out) {
        cerr << "Cannot write site POSCAR file: " << filename << endl;
        return false;
    }

    out << "PAIPAI findinter site visualization\n";
    out << "1.0\n";
    out << structure.a.x << " " << structure.a.y << " " << structure.a.z << "\n";
    out << structure.b.x << " " << structure.b.y << " " << structure.b.z << "\n";
    out << structure.c.x << " " << structure.c.y << " " << structure.c.z << "\n";

    for (int i = 0; i < (int)structure.types.size(); i++) {
        out << PERIODICTABLE[structure.types[i]-1] << " ";
    }
    out << "H\n";

    for (int i = 0; i < (int)structure.counts.size(); i++) {
        out << structure.counts[i] << " ";
    }
    out << sites.size() << "\n";

    out << "Cartesian\n";
    for (const Vec3& r : structure.cart) {
        out << r.x << " " << r.y << " " << r.z << "\n";
    }
    for (const Candidate& site : sites) {
        out << site.cart.x << " " << site.cart.y << " " << site.cart.z << "\n";
    }
    return true;
}

} // namespace

int run_findinter(const FindInterOptions& options) {
    MetalStructure structure;
    if (!read_metal_poscar(options.input_poscar, structure)) return 1;

    map<string, Real> radii;
    if (!ensure_radii(options, structure, radii)) return 2;

    cout << "Running findinter grid search...\n";
    cout << "  input = " << options.input_poscar << "\n";
    cout << "  output = " << options.output_struc << "\n";
    cout << "  grid = " << options.grid << "^3\n";
    cout << "  min_void_factor = " << options.min_void_factor << "\n";
    cout << "  max_void_factor = " << options.max_void_factor << "\n";

    vector<Candidate> sites = find_grid_maxima(options, structure, radii);
    cout << "Found " << sites.size() << " candidate interstitial sites.\n";

    if (!write_struc_in(options, structure, sites)) return 3;
    cout << "Wrote " << options.output_struc << "\n";
    if (!options.site_poscar.empty()) {
        if (!write_site_poscar(options.site_poscar, structure, sites)) return 4;
        cout << "Wrote site visualization POSCAR " << options.site_poscar << "\n";
    }
    return 0;
}
