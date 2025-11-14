#ifndef STRUCTURE_H
#define STRUCTURE_H
#include <string>
#include <vector>
#include "element.h"
using namespace std;
class Structure{
	public:
		int num_metallic_atoms;
		int num_interstitial;
		int num_elements;
		int num_interstitial_elements;
		vector<int> type;
		vector<int> num_atoms;
		vector<int> atomtype;
		vector<int> interstitial_type;
		vector<int> interstitial_postype;
		vector<int> num_interstitial_atoms;
		vector<vector<Real>> interstitial_pos;
		vector<vector<Real>> pos;
		Real cell_x1, cell_x2, cell_x3, cell_y1, cell_y2, cell_y3, cell_z1, cell_z2, cell_z3;
		Structure();
		int readstruc(const char* filename);
		void outputvasp(const char* filename);
		void outputsave(const char* filename);
		void shuffle();
		Real calculateEnergyORB();
		Real relaxedEnergyORB();
		int swapMetal(int a, int b);
		int exchangeMetal(int a, int type);
		int swapInterstitial(int a, int b);
		int exchangeInterstitial(int a, int type);
};
#endif

