#include "header.h"



double compute_filament_energy();
double compute_filament_energy_diff(int ranchain,int ranmonomer,vector<double> & dx);


extern double r[M][N][dim];




double compute_filament_energy(){
	double E=0;
	double u1[dim]={0},u2[dim]={0};
	for(int i=0;i<M;i++){
		double sep=0;
		for(int k=0;k<dim;k++){
			u1[k]=r[i][1][k]-r[i][0][k];
			sep+=u1[k]*u1[k];
		}
		sep=sqrt(sep);
		for(int k=0;k<dim;k++){
			u1[k]/=sep;
		}
		//have first normalized bond vector.
		E+=kfs/2*(sep-aback)*(sep-aback);
		
		for(int j=1;j<N-1;j++){
			double sep=0;
			for(int k=0;k<dim;k++){
				u2[k]=r[i][j+1][k]-r[i][j][k];
				sep+=u2[k]*u2[k];
			}
			sep=sqrt(sep);
			for(int k=0;k<dim;k++){
				u2[k]/=sep;
			}
			//have second normalized bond vector.
			
			E+=kfs/2*(sep-aback)*(sep-aback);
			for(int k=0;k<dim;k++){
				E-=kfb*u1[k]*u2[k];
				//dot product of the bond vectors
			}
			E+=kfb;  //arbitrary offset, if the bonds are straight the energy is zero.
			
			for(int k=0;k<dim;k++){
				u1[k]=u2[k];
			}
			//swap the new bond vector into the old.

		}
		
	}
	
	return E;
	
	
}

double compute_filament_energy_diff(int ranchain,int ranmonomer,vector<double> & dx){
	double dE=0;
	double u1[dim]={0},u2[dim]={0};
	int imin=max(0,ranmonomer-2);
	int imax=min(N-1,ranmonomer+2);
	
	double sep=0;
	
	sep=0;
	for(int k=0;k<dim;k++){
		u1[k]=r[ranchain][imin+1][k]-r[ranchain][imin][k];
		sep+=u1[k]*u1[k];
	}
	sep=sqrt(sep);
	for(int k=0;k<dim;k++){
		u1[k]/=sep;
	}
	//have first normalized bond vector.
	dE-=kfs/2*(sep-aback)*(sep-aback);
	
	for(int i=imin+1;i<imax;i++){
		sep=0;
		for(int k=0;k<dim;k++){
			u2[k]=r[ranchain][i+1][k]-r[ranchain][i][k];
			sep+=u2[k]*u2[k];
		}
		sep=sqrt(sep);
		for(int k=0;k<dim;k++){
			u2[k]/=sep;
		}
		//have second normalized bond vector.
		
		dE-=kfs/2*(sep-aback)*(sep-aback);
		for(int k=0;k<dim;k++){
			dE+=kfb*u1[k]*u2[k];
			//dot product of the bond vectors
		}
		
		for(int k=0;k<dim;k++){
			u1[k]=u2[k];
		}
		//swap the new bond vector into the old.
	}
	
	for(int k=0;k<dim;k++){
		r[ranchain][ranmonomer][k]+=dx[k];  //here we're updating the positions of the new monomer.  We have to undo this later
	}
	//update the bond's position

	
	
	sep=0;
	for(int k=0;k<dim;k++){
		u1[k]=r[ranchain][imin+1][k]-r[ranchain][imin][k];
		sep+=u1[k]*u1[k];
	}
	sep=sqrt(sep);
	for(int k=0;k<dim;k++){
		u1[k]/=sep;
	}
	//have first normalized bond vector.
	dE+=kfs/2*(sep-aback)*(sep-aback);
	
	for(int i=imin+1;i<imax;i++){
		sep=0;
		for(int k=0;k<dim;k++){
			u2[k]=r[ranchain][i+1][k]-r[ranchain][i][k];
			sep+=u2[k]*u2[k];
		}
		sep=sqrt(sep);
		for(int k=0;k<dim;k++){
			u2[k]/=sep;
		}
		//have second normalized bond vector.
		
		dE+=kfs/2*(sep-aback)*(sep-aback);
		for(int k=0;k<dim;k++){
			dE-=kfb*u1[k]*u2[k];
			//dot product of the bond vectors
		}
		for(int k=0;k<dim;k++){
			u1[k]=u2[k];
		}
		//swap the new bond vector into the old.
	}
	
	
	for(int k=0;k<dim;k++){
		r[ranchain][ranmonomer][k]-=dx[k];  //undo the update
	}

	
	
	
	
	
	return dE;
	
	
}



























	
	
