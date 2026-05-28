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
		vector<vector<int>> intsite_metal_neighbors;
		vector<vector<int>> intsite_hop_neighbors;
		Real cell_x1, cell_x2, cell_x3, cell_y1, cell_y2, cell_y3, cell_z1, cell_z2, cell_z3;

		Structure();
		int readstruc(const char* filename);
		void outputvasp(const char* filename);
		void outputsave(const char* filename);
		void shuffle();

		void buildIntsiteMetalNeighborMap(Real cutoff);
		void outputIntsiteMetalNeighborMap(const char* filename);
		int readIntsiteMetalNeighborMap(const char* filename);
		void buildIntsiteHopNeighborMap(Real cutoff);
		void outputIntsiteHopNeighborMap(const char* filename);
		int readIntsiteHopNeighborMap(const char* filename);
		int countInterstitialHopEdges() const;
		int localHopInterstitialRandom(int& a, int& b, int& forward_choices, int& reverse_choices);
		int updateCoordinatesFromContcar(const char* contcar_filename,
							 const char* intsite_neighbor_filename);
		int seedTrialCoordinatesFromContcar(const Structure& current_ref,
							 const char* contcar_filename,
							 const char* intsite_neighbor_filename);
		int reconcileInterstitialSitesFromContcar(const char* contcar_filename,
							 const char* intsite_neighbor_filename,
							 int allow_reassign,
							 Real max_site_distance,
							 const char* reordered_contcar_filename,
							 int& n_reassigned);

		Real calculateEnergyORB();
		Real relaxedEnergyORB();
		int swapMetal(int a, int b);
		int exchangeMetal(int a, int type);
		int swapInterstitial(int a, int b);
		int clusterSwapInterstitial(int a, int b);
		int exchangeInterstitial(int a, int type);
};
#endif
