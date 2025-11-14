#include "structure.h"
#include <string>
#include <fstream>
#include <vector>
#include <numeric>
#include <iostream>
#include <time.h>
#include <math.h>
#include <cstdlib>
using namespace std;
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
        vector<int>().swap(atomtype);       // element type of each atom (as the order in the Species Names, from 0 to num_elements-1)
        vector<int>().swap(num_atoms);          //number of atoms of each species
        vector<int>().swap(type);               //the atomic number of each species
        vector<vector<Real>>().swap(interstitial_pos); //interstitial positions, in Cartesian with scalingfactor = 1 (real position)
        vector<int>().swap(interstitial_type);  // the atomic number of each interstitial species
        vector<int>().swap(interstitial_postype);       // the type of each interstitial site occupation. -1 for not occupied, 0 to num_interstitial_elements-1 for the occupation of the corresponding type
        vector<int>().swap(num_interstitial_atoms);     //number of atoms of each interstitial species



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
        getline(inputfile, tempstring);
        getline(inputfile, tempstring);
        string element_name;
        num_elements = 0;
        int atomicnum;
	vector<int>().swap(type);       //clear
        while (pointer != string::npos) {
                pointer = tempstring.find(" ");
                
		element_name = tempstring.substr(0, pointer);
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
	vector<int>().swap(num_atoms);		//clear
	int i,j;
	num_metallic_atoms = 0;
	for (i = 0; i < num_elements; i++){
		inputfile >> tempnum;
		num_atoms.push_back(tempnum);
		num_metallic_atoms += tempnum;
	}

	//Read interstitial species
        pointer = 0;
        string tempstringinter;
        getline(inputfile, tempstringinter);
        getline(inputfile, tempstringinter);
        num_interstitial_elements = 0;
        vector<int>().swap(num_interstitial_atoms);
        vector<int>().swap(interstitial_type);       //clear
        while (pointer != string::npos) {
                pointer = tempstringinter.find(" ");
                element_name = tempstringinter.substr(0, pointer);
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
	// Read metallic atoms positions
	vector<Real> atom_temp(3);
	vector<Real> atom_pos_temp(3);
        vector<vector<Real>>().swap(pos);      // clear
        vector<int>().swap(atomtype);       // clear
	string coordtype;
	num_metallic_atoms = accumulate(num_atoms.begin(),num_atoms.end(),0);
	cout <<"Number of metallic atoms: " << num_metallic_atoms <<endl;
	getline(inputfile, coordtype);
	cout << coordtype;
	if ((coordtype.front() == 'C')||(coordtype.front() == 'c')||(coordtype.front() == 'K')||(coordtype.front() == 'k')){
		//Cartesian coords used
		
		for (i = 0; i < num_elements; i++){
			for (j = 0; j < num_atoms[i]; j++){
				inputfile >> atom_temp[0] >> atom_temp[1] >> atom_temp[2];
				atom_pos_temp[0] = atom_temp[0] * scalingfactor;
				atom_pos_temp[1] = atom_temp[1] * scalingfactor;
				atom_pos_temp[2] = atom_temp[2] * scalingfactor;
				pos.push_back(atom_pos_temp);
				atomtype.push_back(i);
			}
		}
	}
	else{
		cout << "Fractional coords used"<<endl;
		//Direct coords used
		for (i = 0; i < num_elements; i++){
                        for (j = 0; j < num_atoms[i]; j++){
                                inputfile >> atom_temp[0] >> atom_temp[1] >> atom_temp[2];
                                atom_pos_temp[0] = atom_temp[0] * cell_x1 + atom_temp[1] * cell_x2 + atom_temp[2] * cell_x3;
                                atom_pos_temp[1] = atom_temp[0] * cell_y1 + atom_temp[1] * cell_y2 + atom_temp[2] * cell_y3;
                                atom_pos_temp[2] = atom_temp[0] * cell_z1 + atom_temp[1] * cell_z2 + atom_temp[2] * cell_z3;
				pos.push_back(atom_pos_temp);
                                atomtype.push_back(i);
                        }
                }
	}
        cout <<"Metallic structure read successfully" << endl;
	cout << coordtype << endl;
	//Read interstitial site positions
	vector<vector<Real>>().swap(interstitial_pos);		//clear
        if ((coordtype.front() == 'C')||(coordtype.front() == 'c')||(coordtype.front() == 'K')||(coordtype.front() == 'k')){
                //Cartesian coords used
                for (i = 0; i < num_interstitial; i++){
                        inputfile >> atom_temp[0] >> atom_temp[1] >> atom_temp[2];
                        atom_pos_temp[0] = atom_temp[0] * scalingfactor;
                        atom_pos_temp[1] = atom_temp[1] * scalingfactor;
                        atom_pos_temp[2] = atom_temp[2] * scalingfactor;
                        interstitial_pos.push_back(atom_pos_temp);
                	interstitial_postype.push_back(-1);
		}
        }
        else{
                //Direct coords used
                for (i = 0; i < num_interstitial; i++){
                        inputfile >> atom_temp[0] >> atom_temp[1] >> atom_temp[2];
                        atom_pos_temp[0] = atom_temp[0] * cell_x1 + atom_temp[1] * cell_x2 + atom_temp[2] * cell_x3;
                        atom_pos_temp[1] = atom_temp[0] * cell_y1 + atom_temp[1] * cell_y2 + atom_temp[2] * cell_y3;
                        atom_pos_temp[2] = atom_temp[0] * cell_z1 + atom_temp[1] * cell_z2 + atom_temp[2] * cell_z3;
                        interstitial_pos.push_back(atom_pos_temp);
                        interstitial_postype.push_back(-1);
                }
        }
	
	// Set occupied interstitials
	pointer = 0;
	for (i = 0; i < num_interstitial_elements; i++){
		for (j = 0; j < num_interstitial_atoms[i]; j++){
			interstitial_postype[pointer] = i;
			pointer++;
		}
	}
        //shuffle
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
	//output other interstitial sites
	for (i = 0; i < num_interstitial; i++){
		if (interstitial_postype [i] == -1 ){
			outputfile << interstitial_pos[i][0] << " " << interstitial_pos[i][1] << " " << interstitial_pos[i][2] << "\n";
		}
	}
        outputfile.close();
}



void Structure::shuffle(){
	int i,j,a,b;
	int shuffle_times=10;
	srand((unsigned)time(NULL));
	for (i=0;i<shuffle_times;i++){
		a=rand()%num_metallic_atoms;
		b=rand()%num_metallic_atoms;
		cout << "we swap  " << a << " & " <<b<<endl;
		swapMetal(a,b);
	}

	for (i = 0; i < num_interstitial_elements; i++){
		for (j = 0; j < num_interstitial_atoms[i]; j++){
			a = rand()%num_interstitial;
			while (interstitial_postype[a] > -1){
	                        a = rand()%num_interstitial;
			}
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
        if (atomtype[a] == type){
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


