#include "semigrand.h"
#include <cassert>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

static Structure toy() {
    Structure s;
    s.num_elements = 3;
    s.type = {22, 13, 14}; // Ti Al Si
    s.num_atoms = {2, 2, 0};
    s.num_metallic_atoms = 4;
    s.atomtype = {0, 0, 1, 1};
    s.metal_site_sublattice = {0, 0, 1, 1};
    s.num_interstitial_elements = 1;
    s.interstitial_type = {6}; // C
    s.num_interstitial_atoms = {0};
    s.num_interstitial = 2;
    s.interstitial_postype = {-1, -1};
    s.cell_x1 = s.cell_y2 = s.cell_z3 = 8.0;
    s.pos = {{0,0,0},{4,4,0},{4,0,4},{0,4,4}};
    s.interstitial_pos = {{2,2,2},{6,6,6}};
    return s;
}

static void check_counts(const Structure& s) {
    int sum = 0;
    std::vector<int> actual(s.num_elements, 0);
    for (int x : s.atomtype) { assert(x >= 0 && x < s.num_elements); ++actual[x]; }
    for (int i = 0; i < s.num_elements; ++i) {
        sum += s.num_atoms[i];
        assert(s.num_atoms[i] == actual[i]);
    }
    assert(sum == s.num_metallic_atoms);
    std::vector<int> actual_inter(s.num_interstitial_elements, 0);
    for (int x : s.interstitial_postype)
        if (x >= 0) { assert(x < s.num_interstitial_elements); ++actual_inter[x]; }
    for (int i = 0; i < s.num_interstitial_elements; ++i)
        assert(actual_inter[i] == s.num_interstitial_atoms[i]);
}

int main() {
    std::srand(7);
    Structure s = toy();
    SpeciesPolicy fixed;
    fixed.configure(s, {});
    int site = -1, target = -1;
    assert(!fixed.exchange_metal(s, site, target));
    assert(!fixed.exchange_interstitial(s, site, target));

    SpeciesPolicy open;
    open.configure(s, {"Ti=mu:-4.0", "Al=mu:-3.0",
                       "Si=mu:-5.0", "C=mu:-6.0"});
    assert(std::abs(open.chemical_term(s) - (-14.0)) < 1e-10);
    int first_si = 0;
    for (int j = 0; j < 300; ++j) {
        assert(open.exchange_metal(s, site, target));
        assert(open.exchange_interstitial(s, site, target));
        check_counts(s);
        if (s.num_atoms[2] > 0) ++first_si;
    }
    assert(first_si > 0); // The initially absent Si can be inserted on metal sites.
    assert(s.num_interstitial_atoms[0] <= s.num_interstitial);

    // A fixed Ti population survives both swaps and all open-species moves.
    s = toy();
    SpeciesPolicy mixed;
    mixed.configure(s, {"Ti=fixed", "Al=mu:-3", "Si=mu:-5", "C=mu:-6"});
    for (int j = 0; j < 100; ++j) {
        assert(mixed.exchange_metal(s, site, target));
        check_counts(s);
        assert(s.num_atoms[0] == 2);
    }
    assert(std::abs(mixed.chemical_term(s) -
           (-3.0*s.num_atoms[1]-5.0*s.num_atoms[2]-6.0*s.num_interstitial_atoms[0])) < 1e-10);

    // Only one open substitutional element cannot change its number while
    // the remaining fully occupied lattice species are all canonical.
    s = toy();
    SpeciesPolicy impossible;
    impossible.configure(s, {"Si=mu:-5"});
    assert(!impossible.exchange_metal(s, site, target));

    bool threw = false;
    try { impossible.configure(s, {"Xx=mu:-1"}); }
    catch (const std::runtime_error&) { threw = true; }
    assert(threw);

    // Writing a substitutional-only SAVE/POSCAR must not index element [-1].
    s = toy();
    s.num_interstitial_elements = 0;
    s.interstitial_type.clear();
    s.num_interstitial_atoms.clear();
    s.interstitial_postype.clear();
    s.interstitial_pos.clear();
    s.num_interstitial = 0;
    s.outputsave("v22_test_SAVE");
    s.outputvasp("v22_test_POSCAR");
    Structure recovered;
    assert(recovered.readstruc("v22_test_SAVE") == 1);
    assert(recovered.num_interstitial_elements == 0);
    assert(recovered.num_atoms == s.num_atoms);
    assert(recovered.metal_site_sublattice == s.metal_site_sublattice);
    std::remove("v22_test_SAVE");
    std::remove("v22_test_POSCAR");

    std::cout << "v2.2 species-policy checks passed\n";
    return 0;
}
