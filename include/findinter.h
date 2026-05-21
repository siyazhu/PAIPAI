#ifndef FINDINTER_H
#define FINDINTER_H

#include <map>
#include <string>
#include <vector>
#include "constants.h"

struct FindInterOptions {
    std::string input_poscar = "POSCAR";
    std::string output_struc = "struc.in";
    std::string site_poscar;
    std::string radii_file = "radii.dat";
    std::vector<int> interstitial_types;
    std::vector<int> interstitial_counts;
    int grid_nx = 0;
    int grid_ny = 0;
    int grid_nz = 0;
    Real grid_step = 0.0;
    int max_sites = 0;
    Real min_void_factor = 1.0;
    Real max_void_factor = 2.0;
    Real merge_distance = 0.0;
};

int run_findinter(const FindInterOptions& options);

#endif
