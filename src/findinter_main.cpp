#include "findinter.h"
#include "element.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace {

void print_help(const char* prog) {
    cout << "Usage: " << prog
         << " --input POSCAR --inter B,O --internum 4,5 [options]\n\n"
         << "Options:\n"
         << "  --output FILE                Output PAIPAI struc.in file (default: struc.in)\n"
         << "  --site-poscar FILE           Optional POSCAR with H atoms marking every found interstitial site\n"
         << "  --radii FILE                 Radii file (default: radii.dat)\n"
         << "  --grid N                     Fractional grid per axis (default: 50)\n"
         << "  --max-sites N                Maximum sites to write (default: 2000)\n"
         << "  --min-void-factor X          Require d_nearest >= X*(r_metal+r_inter); default: 1.0\n"
         << "  --max-void-factor X          Reject sites farther than X*(r_metal+r_inter) from nearest metal (default: 2.0)\n"
         << "  --merge-distance X           Optional extra merge distance for nearly duplicate selected sites (default: 0)\n"
         << "  -h, --help                   Show this help message\n";
}

vector<string> split_csv(const string& text) {
    vector<string> out;
    string item;
    stringstream ss(text);
    while (getline(ss, item, ',')) {
        if (!item.empty()) out.push_back(item);
    }
    return out;
}

string value_after_equal_or_next(const string& arg, int& i, int argc, char** argv) {
    size_t pos = arg.find('=');
    if (pos != string::npos) return arg.substr(pos + 1);
    if (i + 1 >= argc) {
        cerr << "Missing value for " << arg << endl;
        exit(2);
    }
    return argv[++i];
}

} // namespace

int main(int argc, char** argv) {
    FindInterOptions options;
    string inter_text;
    string internum_text;

    if (argc == 1) {
        print_help(argv[0]);
        return 2;
    }

    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_help(argv[0]);
            return 0;
        } else if (arg == "--input" || arg.rfind("--input=", 0) == 0) {
            options.input_poscar = value_after_equal_or_next(arg, i, argc, argv);
        } else if (arg == "--output" || arg.rfind("--output=", 0) == 0) {
            options.output_struc = value_after_equal_or_next(arg, i, argc, argv);
        } else if (arg == "--site-poscar" || arg.rfind("--site-poscar=", 0) == 0) {
            options.site_poscar = value_after_equal_or_next(arg, i, argc, argv);
        } else if (arg == "--radii" || arg.rfind("--radii=", 0) == 0) {
            options.radii_file = value_after_equal_or_next(arg, i, argc, argv);
        } else if (arg == "--inter" || arg.rfind("--inter=", 0) == 0) {
            inter_text = value_after_equal_or_next(arg, i, argc, argv);
        } else if (arg == "--internum" || arg.rfind("--internum=", 0) == 0) {
            internum_text = value_after_equal_or_next(arg, i, argc, argv);
        } else if (arg == "--grid" || arg.rfind("--grid=", 0) == 0) {
            options.grid = max(3, stoi(value_after_equal_or_next(arg, i, argc, argv)));
        } else if (arg == "--max-sites" || arg.rfind("--max-sites=", 0) == 0) {
            options.max_sites = max(1, stoi(value_after_equal_or_next(arg, i, argc, argv)));
        } else if (arg == "--min-void-factor" || arg.rfind("--min-void-factor=", 0) == 0) {
            options.min_void_factor = stod(value_after_equal_or_next(arg, i, argc, argv));
        } else if (arg == "--max-void-factor" || arg.rfind("--max-void-factor=", 0) == 0) {
            options.max_void_factor = stod(value_after_equal_or_next(arg, i, argc, argv));
        } else if (arg == "--merge-distance" || arg.rfind("--merge-distance=", 0) == 0) {
            options.merge_distance = max(0.0, stod(value_after_equal_or_next(arg, i, argc, argv)));
        } else {
            cerr << "Unknown option: " << arg << endl;
            print_help(argv[0]);
            return 2;
        }
    }

    if (inter_text.empty() || internum_text.empty()) {
        cerr << "--inter and --internum are required." << endl;
        print_help(argv[0]);
        return 2;
    }

    vector<string> inter_names = split_csv(inter_text);
    vector<string> count_texts = split_csv(internum_text);
    if (inter_names.size() != count_texts.size()) {
        cerr << "--inter and --internum must have the same number of entries." << endl;
        return 2;
    }

    for (const string& name : inter_names) {
        int z = elementnum(name);
        if (z == 0) {
            cerr << "Unknown interstitial element: " << name << endl;
            return 2;
        }
        options.interstitial_types.push_back(z);
    }
    for (const string& count : count_texts) {
        options.interstitial_counts.push_back(max(0, stoi(count)));
    }

    if (options.min_void_factor <= 0.0) {
        cerr << "--min-void-factor must be > 0." << endl;
        return 2;
    }
    if (options.max_void_factor <= 0.0) {
        cerr << "--max-void-factor must be > 0." << endl;
        return 2;
    }
    if (options.min_void_factor > options.max_void_factor) {
        cerr << "--min-void-factor must be <= --max-void-factor." << endl;
        return 2;
    }

    return run_findinter(options);
}
