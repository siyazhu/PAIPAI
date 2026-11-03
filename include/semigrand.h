#ifndef PAIPAI_SEMIGRAND_H
#define PAIPAI_SEMIGRAND_H

// PAIPAI v2.2 mixed canonical / semi-grand-canonical species controls.
// A fully occupied substitutional lattice has a fixed total site count;
// interstitial sites may additionally be empty.  Energies are in eV.
#include "structure.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct SpeciesPolicy {
    std::vector<bool> metal_open, inter_open;
    std::vector<double> metal_mu, inter_mu;

    static std::string species_name(int atomic_number) {
        return std::string(PERIODICTABLE[atomic_number - 1]);
    }

    void configure(const Structure& s, const std::vector<std::string>& specs) {
        metal_open.assign(s.num_elements, false);
        inter_open.assign(s.num_interstitial_elements, false);
        metal_mu.assign(s.num_elements, 0.0);
        inter_mu.assign(s.num_interstitial_elements, 0.0);
        std::vector<bool> metal_seen(s.num_elements, false);
        std::vector<bool> inter_seen(s.num_interstitial_elements, false);

        for (const std::string& spec : specs) {
            const size_t eq = spec.find('=');
            if (eq == std::string::npos || eq == 0 || eq + 1 == spec.size())
                throw std::runtime_error("Invalid --species rule: " + spec +
                                         " (expected Ti=fixed or Si=mu:-4.5)");
            const std::string name = spec.substr(0, eq);
            const std::string value = spec.substr(eq + 1);
            bool open = false;
            double mu = 0.0;
            if (value == "fixed") {
                open = false;
            } else if (value.rfind("mu:", 0) == 0) {
                size_t n = 0;
                try {
                    mu = std::stod(value.substr(3), &n);
                } catch (...) {
                    throw std::runtime_error("Invalid chemical potential in " + spec);
                }
                if (n != value.size() - 3 || !std::isfinite(mu))
                    throw std::runtime_error("Invalid chemical potential in " + spec);
                open = true;
            } else {
                throw std::runtime_error("Invalid --species policy: " + value);
            }

            bool found = false;
            for (int i = 0; i < s.num_elements; ++i) {
                if (species_name(s.type[i]) != name) continue;
                if (metal_seen[i]) throw std::runtime_error("Duplicate species rule: " + name);
                metal_seen[i] = true;
                metal_open[i] = open;
                metal_mu[i] = mu;
                found = true;
            }
            for (int i = 0; i < s.num_interstitial_elements; ++i) {
                if (species_name(s.interstitial_type[i]) != name) continue;
                if (inter_seen[i]) throw std::runtime_error("Duplicate species rule: " + name);
                inter_seen[i] = true;
                inter_open[i] = open;
                inter_mu[i] = mu;
                found = true;
            }
            if (!found) {
                throw std::runtime_error("Unknown --species element " + name +
                                         "; add it to the struc.in species header (count may be zero)");
            }
        }
    }

    double chemical_term(const Structure& s) const {
        double term = 0.0;
        for (int i = 0; i < s.num_elements; ++i)
            if (metal_open[i]) term += metal_mu[i] * s.num_atoms[i];
        for (int i = 0; i < s.num_interstitial_elements; ++i)
            if (inter_open[i]) term += inter_mu[i] * s.num_interstitial_atoms[i];
        return term;
    }

    // Pick uniformly among all currently open-occupied substitutional sites,
    // then uniformly among the other open species.  The number of open sites
    // is invariant under these transmutations. Consequently the forward and
    // reverse proposal probabilities are identical, even at variable composition.
    // Fixed species can still move by regular swapMetal (which preserves counts).
    bool exchange_metal(Structure& s, int& site, int& to_type) const {
        std::vector<int> types, sites;
        for (int i = 0; i < s.num_elements; ++i)
            if (metal_open[i]) types.push_back(i);
        if (types.size() < 2) return false;
        for (int i = 0; i < s.num_metallic_atoms; ++i)
            if (metal_open[s.atomtype[i]]) sites.push_back(i);
        if (sites.empty()) return false;
        site = sites[std::rand() % sites.size()];
        std::vector<int> alternatives;
        for (int t : types) if (t != s.atomtype[site]) alternatives.push_back(t);
        to_type = alternatives[std::rand() % alternatives.size()];
        return s.exchangeMetal(site, to_type) == 1;
    }

    // Empty (-1) is an always-available state at interstitial sites.  Fixed
    // interstitial species cannot be deleted; C=mu enables empty <-> C moves.
    // Uniform over mutable sites and alternative open states is reversible.
    bool exchange_interstitial(Structure& s, int& site, int& to_type) const {
        std::vector<int> states = {-1};
        for (int i = 0; i < s.num_interstitial_elements; ++i)
            if (inter_open[i]) states.push_back(i);
        if (states.size() < 2) return false;
        std::vector<int> sites;
        for (int i = 0; i < s.num_interstitial; ++i) {
            int t = s.interstitial_postype[i];
            if (t == -1 || (t >= 0 && inter_open[t])) sites.push_back(i);
        }
        if (sites.empty()) return false;
        site = sites[std::rand() % sites.size()];
        std::vector<int> alternatives;
        for (int t : states) if (t != s.interstitial_postype[site]) alternatives.push_back(t);
        to_type = alternatives[std::rand() % alternatives.size()];
        return s.exchangeInterstitial(site, to_type) == 1;
    }

    std::string counts(const Structure& s) const {
        std::ostringstream out;
        for (int i = 0; i < s.num_elements; ++i)
            out << (i ? "," : "") << species_name(s.type[i]) << "=" << s.num_atoms[i];
        for (int i = 0; i < s.num_interstitial_elements; ++i)
            out << "," << species_name(s.interstitial_type[i]) << "(inter)="
                << s.num_interstitial_atoms[i];
        return out.str();
    }

    bool has_open_species() const {
        for (bool x : metal_open) if (x) return true;
        for (bool x : inter_open) if (x) return true;
        return false;
    }
};

#endif // PAIPAI_SEMIGRAND_H
