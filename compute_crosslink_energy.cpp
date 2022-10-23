#include "header.h"



double compute_crosslink_energy();
double compute_crosslink_energy_diff(int ranchain,int ranmonomer,vector<double> & dx);
double crosslink_bond_energy(int c1,int c2,int m);


extern double r[M][N][dim];
extern vector <vector<pair<int,int> > > links;

double bond[dim]={0},ueff[dim]={0};
vector<double> tmp;

double compute_crosslink_energy_diff(int ranchain,int ranmonomer,vector<double> & dx){
	
	/*
	double dE=0;
	dE-=compute_crosslink_energy();
	for(int k=0;k<dim;k++){
		r[ranchain][ranmonomer][k]+=dx[k];
	}
	dE+=compute_crosslink_energy();
	for(int k=0;k<dim;k++){
		r[ranchain][ranmonomer][k]-=dx[k];
	}
	
	return dE;
	*/
	
	
	
	
	double dE=0;
	

	int imin=max(0,ranmonomer-1);  //bond perturbation can change the current cross section (ranmonomer),
	int imax=min(N-1,ranmonomer+1);  // or change the net filament direction of the crosssections above or below.
	for(int i=imin;i<=imax;i++){
		for(int b=0;b<links[i].size();b++){
			int c1=links[i][b].first;
			int c2=links[i][b].second;  //for every filament on this cross section,
			dE-=crosslink_bond_energy(c1,c2,i);
		}
	}
	for(int k=0;k<dim;k++){
		r[ranchain][ranmonomer][k]+=dx[k];  //here we're updating the positions of the new monomer.  We have to undo this later
	}
	//update the bond's position
	for(int i=imin;i<=imax;i++){
		for(int b=0;b<links[i].size();b++){
			int c1=links[i][b].first;
			int c2=links[i][b].second;
			dE+=crosslink_bond_energy(c1,c2,i);
		}
	}
	for(int k=0;k<dim;k++){
		r[ranchain][ranmonomer][k]-=dx[k];  	//undo the update

	}


	
	
	
	return dE;
	
}

double compute_crosslink_energy(){
	double E=0;
	for(int i=0;i<links.size();i++){
		for(int b=0;b<links[i].size();b++){
			int c1=links[i][b].first;
			int c2=links[i][b].second;
			E+=crosslink_bond_energy(c1,c2,i);
			//this computes the energy of each crosslink bond in the bundle, one at a time.
		}
	}
	return E;
}




double crosslink_bond_energy(int c1,int c2,int m){
	//c2 is the filament with the perturbation, NOT c1!
	double Ebond=0;
	
	
	double mag=0;
	for(int k=0;k<dim;k++){
		bond[k]=r[c1][m][k]-r[c2][m][k];
		mag+=bond[k]*bond[k];
	}
	mag=sqrt(mag);
	for(int k=0;k<dim;k++){
		bond[k]/=mag;
	}
	Ebond+=kcs/2*(mag-across)*(mag-across);
	//we have computed the stretching energy.
	
	
	return Ebond;
}
