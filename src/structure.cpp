#include "structure.h"
#include <string>
#include <fstream>
#include <vector>
#include <numeric>
#include <iostream>
#include <time.h>
#include <math.h>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cctype>
using namespace std;


static bool read_nonempty_line(ifstream& inputfile, string& line){
    while (getline(inputfile, line)){
        if (line.size() == 0) continue;
        bool only_space = true;
        for (char c : line){
            if (!isspace(static_cast<unsigned char>(c))){
                only_space = false;
                break;
            }
        }
        if (!only_space) return true;
    }
    return false;
}

static bool parse_coord_occupation_line(const string& line, Real& x, Real& y, Real& z, int& occ, bool& has_occ){
    stringstream ss(line);
    if (!(ss >> x >> y >> z)) return false;
    if (ss >> occ) has_occ = true;
    else {
        occ = 0;
        has_occ = false;
    }
    return true;
}

static vector<Real> make_vec3(Real x, Real y, Real z){
    vector<Real> v(3);
    v[0] = x; v[1] = y; v[2] = z;
    return v;
}

static Real dot3(const vector<Real>& a, const vector<Real>& b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static vector<Real> matvec_cell(Real f1, Real f2, Real f3,
                                const vector<Real>& a,
                                const vector<Real>& b,
                                const vector<Real>& c){
    return make_vec3(f1*a[0] + f2*b[0] + f3*c[0],
                     f1*a[1] + f2*b[1] + f3*c[1],
                     f1*a[2] + f2*b[2] + f3*c[2]);
}

static bool inverse3x3_columns(const vector<Real>& a, const vector<Real>& b, const vector<Real>& c,
                               Real inv[3][3]){
    // Matrix columns are a, b, c.
    Real m00=a[0], m01=b[0], m02=c[0];
    Real m10=a[1], m11=b[1], m12=c[1];
    Real m20=a[2], m21=b[2], m22=c[2];
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

static vector<Real> cart_to_frac(const vector<Real>& r, Real inv[3][3]){
    return make_vec3(inv[0][0]*r[0] + inv[0][1]*r[1] + inv[0][2]*r[2],
                     inv[1][0]*r[0] + inv[1][1]*r[1] + inv[1][2]*r[2],
                     inv[2][0]*r[0] + inv[2][1]*r[1] + inv[2][2]*r[2]);
}

static vector<Real> minimum_image_delta(const vector<Real>& delta,
                                        const vector<Real>& a,
                                        const vector<Real>& b,
                                        const vector<Real>& c,
                                        Real inv[3][3]){
    vector<Real> f = cart_to_frac(delta, inv);
    f[0] -= round(f[0]);
    f[1] -= round(f[1]);
    f[2] -= round(f[2]);
    return matvec_cell(f[0], f[1], f[2], a, b, c);
}

static bool line_all_ints(const string& line){
    stringstream ss(line);
    string tok;
    bool any = false;
    while (ss >> tok){
        any = true;
        char* endptr = nullptr;
        strtol(tok.c_str(), &endptr, 10);
        if (endptr == tok.c_str() || *endptr != '\0') return false;
    }
    return any;
}

static bool read_poscar_cart_coords(const char* filename, vector<vector<Real>>& coords,
                                      vector<Real>* out_a = nullptr,
                                      vector<Real>* out_b = nullptr,
                                      vector<Real>* out_c = nullptr){
    coords.clear();
    ifstream inputfile(filename);
    if (!inputfile){
        cerr << "Cannot read POSCAR/CONTCAR file: " << filename << endl;
        return false;
    }
    string line;
    if (!getline(inputfile, line)) return false; // title
    Real scale = 1.0;
    inputfile >> scale;
    vector<Real> a(3), b(3), c(3);
    inputfile >> a[0] >> a[1] >> a[2];
    inputfile >> b[0] >> b[1] >> b[2];
    inputfile >> c[0] >> c[1] >> c[2];
    for (int i=0;i<3;i++){ a[i]*=scale; b[i]*=scale; c[i]*=scale; }
    if (out_a) *out_a = a;
    if (out_b) *out_b = b;
    if (out_c) *out_c = c;
    getline(inputfile, line); // consume end of lattice line

    string species_or_counts;
    if (!read_nonempty_line(inputfile, species_or_counts)) return false;
    string counts_line;
    if (line_all_ints(species_or_counts)){
        counts_line = species_or_counts; // old VASP4 style
    } else {
        if (!read_nonempty_line(inputfile, counts_line)) return false;
    }

    stringstream css(counts_line);
    vector<int> counts;
    int count_tmp;
    int total = 0;
    while (css >> count_tmp){
        counts.push_back(count_tmp);
        total += count_tmp;
    }
    if (total <= 0){
        cerr << "Failed to parse atom counts in " << filename << endl;
        return false;
    }

    string coordtype;
    if (!read_nonempty_line(inputfile, coordtype)) return false;
    if (coordtype.size() > 0 && (coordtype[0]=='S' || coordtype[0]=='s')){
        // Selective dynamics line; next line is coordinate type.
        if (!read_nonempty_line(inputfile, coordtype)) return false;
    }
    bool cartesian = (coordtype.size() > 0 &&
                     (coordtype[0]=='C' || coordtype[0]=='c' || coordtype[0]=='K' || coordtype[0]=='k'));

    coords.reserve(total);
    for (int i=0; i<total; i++){
        if (!read_nonempty_line(inputfile, line)){
            cerr << "Unexpected EOF while reading coordinates in " << filename << endl;
            return false;
        }
        stringstream ss(line);
        Real x,y,z;
        if (!(ss >> x >> y >> z)){
            cerr << "Failed to parse coordinate line in " << filename << ": " << line << endl;
            return false;
        }
        vector<Real> r(3);
        if (cartesian){
            r[0] = x*scale;
            r[1] = y*scale;
            r[2] = z*scale;
        } else {
            r = matvec_cell(x, y, z, a, b, c);
        }
        coords.push_back(r);
    }
    return true;
}

Structure::Structure(){
	cell_x1 = 0;
	cell_x2 = 0;
	cell_x3 = 0;
	cell_y1 = 0;
	cell_y2 = 0;
	cell_y3 = 0;
	cell_z1 = 0;
	cell_z2 = 0;
	cell_z3 = 0;
	num_metallic_atoms = 0;
	num_interstitial = 0;
	num_elements = 0;
	num_interstitial_elements = 0;
	vector<vector<Real>>().swap(pos);      // metallic atom position, in Cartesian with scalingfactor = 1 (real position)
        vector<int>().swap(atomtype);       // element type of each atom (as the order in the Species Names, from 0 to num_elements-1)
	vector<int>().swap(num_atoms);		//number of atoms of each species
	vector<int>().swap(type);		//the atomic number of each species
	vector<vector<Real>>().swap(interstitial_pos); //interstitial positions, in Cartesian with scalingfactor = 1 (real position)
	vector<int>().swap(interstitial_type);	// the atomic number of each interstitial species
	vector<int>().swap(interstitial_postype);	// the type of each interstitial site occupation. -1 for not occupied, 0 to num_interstitial_elements-1 for the occupation of the corresponding type
	vector<int>().swap(num_interstitial_atoms);	//number of atoms of each interstitial species
	vector<vector<int>>().swap(intsite_metal_neighbors);
}


int Structure::readstruc(const char* filename){
	std::ifstream inputfile(filename);
        string elements;
        if (!inputfile) {
                cerr << "Read structure wrong, file not exist!";
                return(0);
        }
	if (!getline(inputfile, elements)) {
                cerr << "File .str is empty!";
                inputfile.close();
		return(0);
        };
	Real scalingfactor;

	vector<vector<Real>>().swap(pos);      // metallic atom position, in Cartesian with scalingfactor = 1 (real position)
        vector<int>().swap(atomtype);       // element type of each atom/site (0 to num_elements-1)
        vector<int>().swap(num_atoms);      // number of atoms of each metallic species
        vector<int>().swap(type);           // atomic number of each metallic species
        vector<vector<Real>>().swap(interstitial_pos); // interstitial site positions
        vector<int>().swap(interstitial_type);         // atomic number of each interstitial species
        vector<int>().swap(interstitial_postype);      // occupation of each interstitial site: -1 empty, 0..num_interstitial_elements-1 occupied
        vector<int>().swap(num_interstitial_atoms);    // number of atoms of each interstitial species
        vector<vector<int>>().swap(intsite_metal_neighbors);

	inputfile >> scalingfactor; 
	
	inputfile >> cell_x1 >> cell_y1 >> cell_z1 >> cell_x2 >> cell_y2 >> cell_z2 >> cell_x3 >> cell_y3 >> cell_z3;
	cell_x1 = cell_x1 * scalingfactor;
	cell_x2 = cell_x2 * scalingfactor;
	cell_x3 = cell_x3 * scalingfactor;
	cell_y1 = cell_y1 * scalingfactor;
	cell_y2 = cell_y2 * scalingfactor;
	cell_y3 = cell_y3 * scalingfactor;
	cell_z1 = cell_z1 * scalingfactor;
	cell_z2 = cell_z2 * scalingfactor;
	cell_z3 = cell_z3 * scalingfactor;

	// Deal with metallic atoms type and number
	int pointer = 0;
        string tempstring;
        getline(inputfile, tempstring); // consume end of cell line
        getline(inputfile, tempstring);
        string element_name;
        num_elements = 0;
        int atomicnum;
	vector<int>().swap(type);       //clear
        while (pointer != string::npos) {
                pointer = tempstring.find(" ");
                element_name = tempstring.substr(0, pointer);
                if (element_name.size() == 0){
                    tempstring.erase(0, pointer + 1);
                    continue;
                }
                cout << "ELEMENT: " << element_name << endl;
		atomicnum = elementnum(element_name);
                if (atomicnum == 0) {
                        cerr << "Element name illegal!";
                        inputfile.close();
			return 0;
                }
                else {
                        type.push_back(atomicnum);
                        num_elements++;
                }
                tempstring.erase(0, pointer + 1);
        }

	int tempnum;
	vector<int>().swap(num_atoms);      //clear
	int i,j;
	num_metallic_atoms = 0;
	for (i = 0; i < num_elements; i++){
		inputfile >> tempnum;
		num_atoms.push_back(tempnum);
		num_metallic_atoms += tempnum;
	}

	// Read interstitial species
        pointer = 0;
        string tempstringinter;
        getline(inputfile, tempstringinter); // consume end of count line
        getline(inputfile, tempstringinter);
        num_interstitial_elements = 0;
        vector<int>().swap(num_interstitial_atoms);
        vector<int>().swap(interstitial_type);       //clear
        while (pointer != string::npos) {
                pointer = tempstringinter.find(" ");
                element_name = tempstringinter.substr(0, pointer);
                if (element_name.size() == 0){
                    tempstringinter.erase(0, pointer + 1);
                    continue;
                }
                cout << "INTERSTITIAL ELEMENT: " << element_name << endl;
                atomicnum = elementnum(element_name);
                if (atomicnum == 0) {
                        cerr << "Interstitial element name illegal!";
                        inputfile.close();
                        return 0;
                }
                else {
                        interstitial_type.push_back(atomicnum);
                        num_interstitial_elements++;
                }
                tempstringinter.erase(0, pointer + 1);
        }

	// Read number of atoms of each interstitial species
	int num_interstitial_temp;
	for (i = 0; i < num_interstitial_elements; i++){
		inputfile >> num_interstitial_temp;
		cout << "Number of " << PERIODICTABLE[interstitial_type[i]-1] << " atoms: " << num_interstitial_temp << endl;
		num_interstitial_atoms.push_back(num_interstitial_temp);
	}

	// Read number of interstitial sites
        inputfile >> num_interstitial;
        cout << "Number of Interstitial site: "<<num_interstitial<<endl;
	string shuffleflag;
	getline(inputfile, shuffleflag);
	getline(inputfile, shuffleflag);

	// Read coordinates
	vector<Real> atom_temp(3);
	vector<Real> atom_pos_temp(3);
        vector<vector<Real>>().swap(pos);      // clear
        vector<int>().swap(atomtype);          // clear
	string coordtype;
	num_metallic_atoms = accumulate(num_atoms.begin(),num_atoms.end(),0);
	cout <<"Number of metallic atoms: " << num_metallic_atoms <<endl;
	getline(inputfile, coordtype);
	cout << coordtype;
	bool cartesian = ((coordtype.front() == 'C')||(coordtype.front() == 'c')||(coordtype.front() == 'K')||(coordtype.front() == 'k'));
	if (!cartesian){
		cout << "Fractional coords used"<<endl;
	}

        bool metal_has_occ = false;
        bool metal_format_checked = false;
        vector<int> metal_count_from_occ(num_elements, 0);
        string line;

        for (i = 0; i < num_metallic_atoms; i++){
            if (!read_nonempty_line(inputfile, line)){
                cerr << "Unexpected end of file while reading metallic coordinates" << endl;
                inputfile.close();
                return 0;
            }
            int occ = 0;
            bool has_occ = false;
            if (!parse_coord_occupation_line(line, atom_temp[0], atom_temp[1], atom_temp[2], occ, has_occ)){
                cerr << "Failed to parse metallic coordinate line: " << line << endl;
                inputfile.close();
                return 0;
            }
            if (!metal_format_checked){
                metal_has_occ = has_occ;
                metal_format_checked = true;
            }

            if (cartesian){
                atom_pos_temp[0] = atom_temp[0] * scalingfactor;
                atom_pos_temp[1] = atom_temp[1] * scalingfactor;
                atom_pos_temp[2] = atom_temp[2] * scalingfactor;
            }
            else{
                atom_pos_temp[0] = atom_temp[0] * cell_x1 + atom_temp[1] * cell_x2 + atom_temp[2] * cell_x3;
                atom_pos_temp[1] = atom_temp[0] * cell_y1 + atom_temp[1] * cell_y2 + atom_temp[2] * cell_y3;
                atom_pos_temp[2] = atom_temp[0] * cell_z1 + atom_temp[1] * cell_z2 + atom_temp[2] * cell_z3;
            }
            pos.push_back(atom_pos_temp);

            if (metal_has_occ){
                if (!has_occ){
                    cerr << "Mixed metallic coordinate format: missing occupation id on line " << i << endl;
                    inputfile.close();
                    return 0;
                }
                if (occ < 0 || occ >= num_elements){
                    cerr << "Metal occupation id out of range on line " << i << ": " << occ << endl;
                    inputfile.close();
                    return 0;
                }
                atomtype.push_back(occ);
                metal_count_from_occ[occ]++;
            }
            else{
                // Backward-compatible old format: metallic coordinates are grouped by element counts.
                int kind = 0;
                int cumulative = 0;
                for (int k = 0; k < num_elements; k++){
                    cumulative += num_atoms[k];
                    if (i < cumulative){
                        kind = k;
                        break;
                    }
                }
                atomtype.push_back(kind);
                metal_count_from_occ[kind]++;
            }
        }

        if (metal_has_occ){
            for (i = 0; i < num_elements; i++){
                if (metal_count_from_occ[i] != num_atoms[i]){
                    cerr << "Warning: metallic count mismatch for " << PERIODICTABLE[type[i]-1]
                         << ": header = " << num_atoms[i]
                         << ", coordinate occupations = " << metal_count_from_occ[i] << endl;
                    num_atoms[i] = metal_count_from_occ[i];
                }
            }
        }

        cout <<"Metallic structure read successfully" << endl;
	cout << coordtype << endl;

	// Read interstitial site positions and occupations
	vector<vector<Real>>().swap(interstitial_pos);      // clear
        vector<int>().swap(interstitial_postype);          // clear
        bool int_has_occ = false;
        bool int_format_checked = false;
        vector<int> int_count_from_occ(num_interstitial_elements, 0);

        for (i = 0; i < num_interstitial; i++){
            if (!read_nonempty_line(inputfile, line)){
                cerr << "Unexpected end of file while reading interstitial coordinates" << endl;
                inputfile.close();
                return 0;
            }
            int occ = -1;
            bool has_occ = false;
            if (!parse_coord_occupation_line(line, atom_temp[0], atom_temp[1], atom_temp[2], occ, has_occ)){
                cerr << "Failed to parse interstitial coordinate line: " << line << endl;
                inputfile.close();
                return 0;
            }
            if (!int_format_checked){
                int_has_occ = has_occ;
                int_format_checked = true;
            }

            if (cartesian){
                atom_pos_temp[0] = atom_temp[0] * scalingfactor;
                atom_pos_temp[1] = atom_temp[1] * scalingfactor;
                atom_pos_temp[2] = atom_temp[2] * scalingfactor;
            }
            else{
                atom_pos_temp[0] = atom_temp[0] * cell_x1 + atom_temp[1] * cell_x2 + atom_temp[2] * cell_x3;
                atom_pos_temp[1] = atom_temp[0] * cell_y1 + atom_temp[1] * cell_y2 + atom_temp[2] * cell_y3;
                atom_pos_temp[2] = atom_temp[0] * cell_z1 + atom_temp[1] * cell_z2 + atom_temp[2] * cell_z3;
            }
            interstitial_pos.push_back(atom_pos_temp);

            if (int_has_occ){
                if (!has_occ){
                    cerr << "Mixed interstitial coordinate format: missing occupation id on line " << i << endl;
                    inputfile.close();
                    return 0;
                }
                if (occ < -1 || occ >= num_interstitial_elements){
                    cerr << "Interstitial occupation id out of range on line " << i << ": " << occ << endl;
                    inputfile.close();
                    return 0;
                }
                interstitial_postype.push_back(occ);
                if (occ >= 0) int_count_from_occ[occ]++;
            }
            else{
                interstitial_postype.push_back(-1);
            }
        }

        if (!int_has_occ){
            // Backward-compatible old format: occupied interstitials are listed first by species counts.
            pointer = 0;
            for (i = 0; i < num_interstitial_elements; i++){
                for (j = 0; j < num_interstitial_atoms[i]; j++){
                    if (pointer >= num_interstitial){
                        cerr << "Interstitial occupation count exceeds number of sites" << endl;
                        inputfile.close();
                        return 0;
                    }
                    interstitial_postype[pointer] = i;
                    pointer++;
                }
            }
        }
        else{
            for (i = 0; i < num_interstitial_elements; i++){
                if (int_count_from_occ[i] != num_interstitial_atoms[i]){
                    cerr << "Warning: interstitial count mismatch for " << PERIODICTABLE[interstitial_type[i]-1]
                         << ": header = " << num_interstitial_atoms[i]
                         << ", coordinate occupations = " << int_count_from_occ[i] << endl;
                    num_interstitial_atoms[i] = int_count_from_occ[i];
                }
            }
        }

        // shuffle
	if (shuffleflag == "Shuffle"){
		cout << "Atomic coordinates shuffled" << endl;
		this->shuffle();
	}
	else{
		cout << "Atomic coordinates retained (No shuffle)" << endl;
	}

	inputfile.close();
	return(1);
}


void Structure::outputvasp(const char* filename){
	ofstream outputfile(filename);
        cout << "OUTPUT TO POSCAR\n";
	int i;
	for (i = 0; i < num_elements ; i++){
                outputfile << PERIODICTABLE[type[i]-1];
        }
	outputfile << " + ";
	for (i = 0; i < num_interstitial_elements; i++){
		outputfile << PERIODICTABLE[interstitial_type[i]-1];
	}
	outputfile << "\n";
	outputfile << "1.0\n";
	outputfile << cell_x1 << " "  << cell_y1 << " "  << cell_z1 << "\n";
	outputfile << cell_x2 << " "  << cell_y2 << " "  << cell_z2 << "\n";
	outputfile << cell_x3 << " "  << cell_y3 << " "  << cell_z3 << "\n";

        for (i = 0; i < num_elements; i++){
                outputfile << PERIODICTABLE[type[i]-1] <<" ";
        }
	for (i = 0; i < num_interstitial_elements - 1; i++){
                outputfile << PERIODICTABLE[interstitial_type[i]-1] << " ";
        }
        outputfile << PERIODICTABLE[interstitial_type[num_interstitial_elements-1]-1] << "\n";

        for (i = 0; i < num_elements; i++){
                outputfile << num_atoms[i] <<" ";
        }
	for (i = 0; i < num_interstitial_elements - 1; i++){
		outputfile << num_interstitial_atoms[i] << " ";	
	}
        outputfile << num_interstitial_atoms[num_interstitial_elements-1] << "\n";

	outputfile << "Cartesian\n";
	//output metallic atoms
	int kind;
	for (kind = 0; kind < num_elements; kind++){

		for (i = 0; i < num_metallic_atoms; i++){
               		if (atomtype[i] == kind){
				outputfile << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << "\n";
			}
		}
	}
	//output interstitial atoms
	for (kind = 0; kind < num_interstitial_elements; kind++){

                for (i = 0; i < num_interstitial; i++){
                        if (interstitial_postype[i] == kind){
                                outputfile << interstitial_pos[i][0] << " " << interstitial_pos[i][1] << " " << interstitial_pos[i][2] << "\n";
                        }
                }
        }
        outputfile.close();

}


void Structure::buildIntsiteMetalNeighborMap(Real cutoff){
        intsite_metal_neighbors.clear();
        intsite_metal_neighbors.resize(num_interstitial);

        vector<Real> a = make_vec3(cell_x1, cell_y1, cell_z1);
        vector<Real> b = make_vec3(cell_x2, cell_y2, cell_z2);
        vector<Real> c = make_vec3(cell_x3, cell_y3, cell_z3);
        Real inv[3][3];
        bool has_inv = inverse3x3_columns(a, b, c, inv);
        if (!has_inv){
                cerr << "Warning: cell matrix is singular; neighbor map will use direct Cartesian distances without PBC." << endl;
        }

        Real cutoff2 = cutoff * cutoff;
        for (int s = 0; s < num_interstitial; s++){
                for (int m = 0; m < num_metallic_atoms; m++){
                        vector<Real> d = make_vec3(interstitial_pos[s][0] - pos[m][0],
                                                   interstitial_pos[s][1] - pos[m][1],
                                                   interstitial_pos[s][2] - pos[m][2]);
                        if (has_inv) d = minimum_image_delta(d, a, b, c, inv);
                        Real d2 = dot3(d,d);
                        if (d2 <= cutoff2){
                                intsite_metal_neighbors[s].push_back(m);
                        }
                }
        }

        int min_n = num_metallic_atoms + 1;
        int max_n = 0;
        int empty_n = 0;
        Real avg_n = 0.0;
        for (int s = 0; s < num_interstitial; s++){
                int n = (int)intsite_metal_neighbors[s].size();
                min_n = min(min_n, n);
                max_n = max(max_n, n);
                if (n == 0) empty_n++;
                avg_n += n;
        }
        if (num_interstitial > 0) avg_n /= (Real)num_interstitial;
        cout << "Built intsite-metal neighbor map with cutoff = " << cutoff << endl;
        cout << "Neighbor count: min = " << min_n
             << ", max = " << max_n
             << ", avg = " << avg_n
             << ", empty sites = " << empty_n << endl;
}

void Structure::outputIntsiteMetalNeighborMap(const char* filename){
        ofstream outputfile(filename);
        if (!outputfile){
                cerr << "Cannot write intsite-metal neighbor map: " << filename << endl;
                return;
        }
        outputfile << "# intsite_id n_neighbors metal_site_ids...\n";
        for (int s = 0; s < (int)intsite_metal_neighbors.size(); s++){
                outputfile << s << " " << intsite_metal_neighbors[s].size();
                for (int j = 0; j < (int)intsite_metal_neighbors[s].size(); j++){
                        outputfile << " " << intsite_metal_neighbors[s][j];
                }
                outputfile << "\n";
        }
        outputfile.close();
}

int Structure::readIntsiteMetalNeighborMap(const char* filename){
        ifstream inputfile(filename);
        if (!inputfile){
                cerr << "Cannot read intsite-metal neighbor map: " << filename << endl;
                return 0;
        }
        intsite_metal_neighbors.clear();
        intsite_metal_neighbors.resize(num_interstitial);

        string line;
        while (getline(inputfile, line)){
                if (line.size() == 0) continue;
                if (line[0] == '#') continue;
                stringstream ss(line);
                int site_id, n_neighbors;
                if (!(ss >> site_id >> n_neighbors)) continue;
                if (site_id < 0 || site_id >= num_interstitial){
                        cerr << "Invalid intsite id in neighbor map: " << site_id << endl;
                        inputfile.close();
                        return 0;
                }
                intsite_metal_neighbors[site_id].clear();
                for (int k = 0; k < n_neighbors; k++){
                        int metal_id;
                        if (!(ss >> metal_id)){
                                cerr << "Missing metal id in neighbor map for intsite " << site_id << endl;
                                inputfile.close();
                                return 0;
                        }
                        if (metal_id < 0 || metal_id >= num_metallic_atoms){
                                cerr << "Invalid metal id in neighbor map: " << metal_id << endl;
                                inputfile.close();
                                return 0;
                        }
                        intsite_metal_neighbors[site_id].push_back(metal_id);
                }
        }
        inputfile.close();
        return 1;
}

int Structure::updateCoordinatesFromContcar(const char* contcar_filename, const char* intsite_neighbor_filename){
        vector<vector<Real>> contcar_coords;
        vector<Real> contcar_a, contcar_b, contcar_c;
        if (!read_poscar_cart_coords(contcar_filename, contcar_coords,
                                     &contcar_a, &contcar_b, &contcar_c)){
                return 0;
        }

        if (!readIntsiteMetalNeighborMap(intsite_neighbor_filename)){
                return 0;
        }

        vector<vector<Real>> old_pos = pos;
        vector<vector<Real>> old_int = interstitial_pos;
        vector<vector<Real>> new_pos = pos;
        vector<vector<Real>> new_int = interstitial_pos;
        vector<int> metal_seen(num_metallic_atoms, 0);

        int atom_index = 0;

        // Read ONLY metallic atom coordinates from CONTCAR.
        // Interstitial site coordinates, whether occupied or empty, are NOT read back
        // from CONTCAR. They are updated from the neighboring metal cage displacement.
        for (int kind = 0; kind < num_elements; kind++){
                for (int i = 0; i < num_metallic_atoms; i++){
                        if (atomtype[i] == kind){
                                if (atom_index >= (int)contcar_coords.size()){
                                        cerr << "CONTCAR has fewer atoms than expected while reading metal coordinates." << endl;
                                        return 0;
                                }
                                new_pos[i] = contcar_coords[atom_index];
                                metal_seen[i]++;
                                atom_index++;
                        }
                }
        }

        for (int m = 0; m < num_metallic_atoms; m++){
                if (metal_seen[m] != 1){
                        cerr << "Metal site " << m << " was mapped " << metal_seen[m]
                             << " times from CONTCAR. Expected exactly once." << endl;
                        return 0;
                }
        }

        // Sanity check: CONTCAR should contain metals plus occupied interstitial atoms.
        // We intentionally do not read interstitial coordinates back from CONTCAR,
        // but this warning helps catch broken or mismatched CONTCAR files.
        int expected_total_atoms = num_metallic_atoms;
        for (int k = 0; k < num_interstitial_elements; k++){
                expected_total_atoms += num_interstitial_atoms[k];
        }

        if ((int)contcar_coords.size() != expected_total_atoms){
                cerr << "Warning: CONTCAR atom count = " << contcar_coords.size()
                     << ", expected = " << expected_total_atoms << endl;
        }

        vector<Real> old_a = make_vec3(cell_x1, cell_y1, cell_z1);
        vector<Real> old_b = make_vec3(cell_x2, cell_y2, cell_z2);
        vector<Real> old_c = make_vec3(cell_x3, cell_y3, cell_z3);

        vector<Real> new_a = contcar_a;
        vector<Real> new_b = contcar_b;
        vector<Real> new_c = contcar_c;

        Real old_inv[3][3];
        Real new_inv[3][3];

        bool has_old_inv = inverse3x3_columns(old_a, old_b, old_c, old_inv);
        bool has_new_inv = inverse3x3_columns(new_a, new_b, new_c, new_inv);

        if (!has_old_inv || !has_new_inv){
                cerr << "Warning: cell matrix is singular; drift correction and displacement update will not use fractional PBC." << endl;
        }

        // ============================================================
        // Remove global fractional drift of the relaxed metal skeleton.
        //
        // Cell is allowed to relax, so Cartesian displacement contains
        // both cell deformation and actual rigid translation. Therefore
        // the global drift is estimated in fractional coordinates:
        //
        //     drift = average( frac_new_metal - frac_old_metal )
        //
        // with minimum-image treatment.
        //
        // This removes only the overall translation inside the cell,
        // while keeping cell relaxation and local atomic relaxation.
        // ============================================================
        vector<Real> frac_drift(3, 0.0);

        if (has_old_inv && has_new_inv){
                for (int m = 0; m < num_metallic_atoms; m++){
                        vector<Real> f_old = cart_to_frac(old_pos[m], old_inv);
                        vector<Real> f_new = cart_to_frac(new_pos[m], new_inv);

                        vector<Real> df = make_vec3(f_new[0] - f_old[0],
                                                    f_new[1] - f_old[1],
                                                    f_new[2] - f_old[2]);

                        df[0] -= round(df[0]);
                        df[1] -= round(df[1]);
                        df[2] -= round(df[2]);

                        frac_drift[0] += df[0];
                        frac_drift[1] += df[1];
                        frac_drift[2] += df[2];
                }

                frac_drift[0] /= (Real)num_metallic_atoms;
                frac_drift[1] /= (Real)num_metallic_atoms;
                frac_drift[2] /= (Real)num_metallic_atoms;

                cout << "Removing global fractional drift: "
                     << frac_drift[0] << " "
                     << frac_drift[1] << " "
                     << frac_drift[2] << endl;

                // Apply the same fractional drift correction to all relaxed metal atoms.
                for (int m = 0; m < num_metallic_atoms; m++){
                        vector<Real> f = cart_to_frac(new_pos[m], new_inv);

                        f[0] -= frac_drift[0];
                        f[1] -= frac_drift[1];
                        f[2] -= frac_drift[2];

                        f[0] -= floor(f[0]);
                        f[1] -= floor(f[1]);
                        f[2] -= floor(f[2]);

                        new_pos[m] = matvec_cell(f[0], f[1], f[2], new_a, new_b, new_c);
                }
        }

        // Calculate metal displacement after global drift removal.
        // For variable-cell relaxation, use fractional displacement:
        //
        //     df = frac_new_corrected - frac_old
        //     d_cart = new_cell * df
        //
        // This avoids mixing true cell deformation with rigid translation.
        vector<vector<Real>> metal_disp(num_metallic_atoms, vector<Real>(3, 0.0));

        for (int m = 0; m < num_metallic_atoms; m++){
                if (has_old_inv && has_new_inv){
                        vector<Real> f_old = cart_to_frac(old_pos[m], old_inv);
                        vector<Real> f_new = cart_to_frac(new_pos[m], new_inv);

                        vector<Real> df = make_vec3(f_new[0] - f_old[0],
                                                    f_new[1] - f_old[1],
                                                    f_new[2] - f_old[2]);

                        df[0] -= round(df[0]);
                        df[1] -= round(df[1]);
                        df[2] -= round(df[2]);

                        metal_disp[m] = matvec_cell(df[0], df[1], df[2],
                                                    new_a, new_b, new_c);
                }
                else{
                        metal_disp[m] = make_vec3(new_pos[m][0] - old_pos[m][0],
                                                  new_pos[m][1] - old_pos[m][1],
                                                  new_pos[m][2] - old_pos[m][2]);
                }
        }

        // Update ALL interstitial sites from neighboring metal cage displacement.
        // This applies to both occupied and empty sites.
        for (int s = 0; s < num_interstitial; s++){

                if (s >= (int)intsite_metal_neighbors.size() || intsite_metal_neighbors[s].size() == 0){
                        cerr << "Warning: no metal neighbors for interstitial site " << s
                             << "; its coordinate is unchanged." << endl;
                        continue;
                }

                vector<Real> avg(3, 0.0);
                for (int idx = 0; idx < (int)intsite_metal_neighbors[s].size(); idx++){
                        int m = intsite_metal_neighbors[s][idx];
                        avg[0] += metal_disp[m][0];
                        avg[1] += metal_disp[m][1];
                        avg[2] += metal_disp[m][2];
                }

                Real invn = 1.0 / (Real)intsite_metal_neighbors[s].size();

                new_int[s][0] = old_int[s][0] + avg[0] * invn;
                new_int[s][1] = old_int[s][1] + avg[1] * invn;
                new_int[s][2] = old_int[s][2] + avg[2] * invn;
        }

        // Wrap interstitial site positions into the new cell.
        // Metal atoms were already wrapped after drift correction.
        if (has_new_inv){
                for (int s = 0; s < num_interstitial; s++){
                        vector<Real> f = cart_to_frac(new_int[s], new_inv);

                        f[0] -= floor(f[0]);
                        f[1] -= floor(f[1]);
                        f[2] -= floor(f[2]);

                        new_int[s] = matvec_cell(f[0], f[1], f[2],
                                                 new_a, new_b, new_c);
                }
        }

        pos = new_pos;
        interstitial_pos = new_int;

        // SAVE carries the relaxed lattice vectors from CONTCAR.
        cell_x1 = contcar_a[0]; cell_y1 = contcar_a[1]; cell_z1 = contcar_a[2];
        cell_x2 = contcar_b[0]; cell_y2 = contcar_b[1]; cell_z2 = contcar_b[2];
        cell_x3 = contcar_c[0]; cell_y3 = contcar_c[1]; cell_z3 = contcar_c[2];

        return 1;
}

void Structure::outputsave(const char* filename){
        ofstream outputfile(filename);
        cout << "OUTPUT TO SAVE\n";
        int i;
        for (i = 0; i < num_elements ; i++){
                outputfile << PERIODICTABLE[type[i]-1];
        }
        outputfile << " + ";
        for (i = 0; i < num_interstitial_elements; i++){
                outputfile << PERIODICTABLE[interstitial_type[i]-1];
        }
        outputfile << "\n";
        outputfile << "1.0\n";
        outputfile << cell_x1 << " "  << cell_y1 << " "  << cell_z1 << "\n";
        outputfile << cell_x2 << " "  << cell_y2 << " "  << cell_z2 << "\n";
        outputfile << cell_x3 << " "  << cell_y3 << " "  << cell_z3 << "\n";
        for (i = 0; i < num_elements - 1; i++){
                outputfile << PERIODICTABLE[type[i]-1] <<" ";
        }
        outputfile << PERIODICTABLE[type[num_elements -1]-1] << "\n";
        for (i = 0; i < num_elements - 1; i++){
                outputfile << num_atoms[i] <<" ";
        }
        outputfile << num_atoms[num_elements - 1] << "\n";

        for (i = 0; i < num_interstitial_elements - 1; i++){
                outputfile << PERIODICTABLE[interstitial_type[i]-1] << " ";
        }
        outputfile << PERIODICTABLE[interstitial_type[num_interstitial_elements-1]-1] << "\n";

        for (i = 0; i < num_interstitial_elements - 1; i++){
                outputfile << num_interstitial_atoms[i] << " ";
        }
        outputfile << num_interstitial_atoms[num_interstitial_elements-1] << "\n";
        
        outputfile << num_interstitial << "\n";
        outputfile << "No Shuffle\n";
        outputfile << "Cartesian\n";

        // New PAIPAI SAVE format:
        // line order is the fixed site order; the 4th column is the occupation id.
        // metal site i:        x y z atomtype[i]
        // interstitial site s: x y z interstitial_postype[s]  (-1 means empty)
        for (i = 0; i < num_metallic_atoms; i++){
                outputfile << pos[i][0] << " " << pos[i][1] << " " << pos[i][2]
                           << " " << atomtype[i] << "\n";
        }

        for (i = 0; i < num_interstitial; i++){
                outputfile << interstitial_pos[i][0] << " " << interstitial_pos[i][1] << " " << interstitial_pos[i][2]
                           << " " << interstitial_postype[i] << "\n";
        }
        outputfile.close();
}




void Structure::shuffle(){
	int i,j,a,b;
	int shuffle_times=200;
	srand((unsigned)time(NULL));
	for (i=0;i<shuffle_times;i++){
		a=rand()%num_metallic_atoms;
		b=rand()%num_metallic_atoms;
		cout << "we swap  " << a << " & " <<b<<endl;
		swapMetal(a,b);
	}

	for (i = 0; i < num_interstitial; i++){
		interstitial_postype[i] = -1;
	}

	for (i = 0; i < num_interstitial_elements; i++){
		for (j = 0; j < num_interstitial_atoms[i]; j++){
			a = rand()%num_interstitial;
			while (interstitial_postype[a] > -1){
	                        a = rand()%num_interstitial;
			}
			cout << "we add interstitial at position " << a <<endl;
			interstitial_postype[a] = i;
		}
	}
}


int Structure::swapMetal(int a, int b){
	int temptype;
	if (a<0 || b<0 || a>=num_metallic_atoms || b>=num_metallic_atoms){
		cout << "Atom number exceeds when swaping positions of two metallic atoms\n";
        	return(2);	
	}
	if (atomtype[a] == atomtype[b]){
		cout << "Same atom type when swaping positions of two metallic atoms\n";
		return(0);
	}	
	temptype = atomtype[a];
	atomtype[a] = atomtype[b];
	atomtype[b] = temptype;
	return(1);
}

int Structure::exchangeMetal(int a, int type){
	if (a<0 || a>=num_metallic_atoms){
		cout << "Atom number exceeds when exchanging the type of metallic atom\n";
		return(2);	
	}
	if (type < 0 || type >= num_elements){
		cout << "Type number exceeds when exchanging the type of metallic atom\n";
		return(3);
	}
	if (atomtype[a] == type){
		cout << "Same atom type when exchanging when exchanging the type of metallic atom\n";
		return(0);
	}
	num_atoms[atomtype[a]]--;
	atomtype[a] = type;
	num_atoms[type]++;
	return(1);
}

int Structure::swapInterstitial(int a, int b){
        int temptype;
        if (a<0 || b<0 || a>=num_interstitial || b>=num_interstitial){
                cout << "Atom number exceeds when swaping positions of two interstitials\n";
                return(2);
        }       
        if (interstitial_postype[a] == interstitial_postype[b]){
                cout << "Same interstitial type when swaping type of two interstitials\n";
                return(0);
        }       
        temptype = interstitial_postype[a];
        interstitial_postype[a] = interstitial_postype[b];
        interstitial_postype[b] = temptype;
        return(1);
}

int Structure::exchangeInterstitial(int a, int type){
        if (a<0 || a>=num_interstitial){
                cout << "Atom number exceeds when exchanging the type of interstitial\n";
                return(2);
        }
        if (type < -1 || type >= num_interstitial_elements){
                cout << "Type number exceeds when exchanging the type of interstitial\n";
                return(3);
        }
        if (interstitial_postype[a] == type){
                cout << "Same interstitial type when exchanging the type of interstitial\n";
                return(0);
        }
	if (interstitial_postype[a] > -1){
        // if occupied by an atom, not an empty site
		num_interstitial_atoms[interstitial_postype[a]]--;
	}
	
	if (type > -1){
		num_interstitial_atoms[type]++;
	}
	interstitial_postype[a] = type;
	return(1);
}

Real Structure::calculateEnergyORB(){
	outputvasp("POSCAR");
	system("python calc.py");
	std::ifstream inputfile("energy");
	Real energy;
	inputfile >> energy;
	return (energy);
}

Real Structure::relaxedEnergyORB(){
        outputvasp("POSCAR");
	system("python relax.py");
        std::ifstream inputfile("energy");
        Real energy;
        inputfile >> energy;
        return (energy);
}


