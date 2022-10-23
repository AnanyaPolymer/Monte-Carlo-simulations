#include "header.h"


extern double r[M][N][dim];
extern vector <vector<pair<int,int> > > links;
extern double ftop,fbottom;

double compute_forces(double *fall){
	
	//there are four terms that affect the forces at the top and bottom:
	// (1) the backbone link between 0 and 1 or N-1 and N-2 on each filament
	// (2) the bending term between 0,1,2 and N-1,N-2,N-3
	// (3) the LJ interactions between 0 or N-1 and the rest of the chain
	// (4) the force due to 0 or N-1 changing the bond direction for crosslinks.
	
	//double ft=0,fb=0;
	for(int k=0;k<2;k++){
		fall[k] = 0.0; //fall[0] = fb, fall[1] = ft
	}
	double mag,mag1,mag2,dp;
	double v[dim]={0};

	for(int i=0;i<M;i++){
		
		
		mag=0;
		for(int k=0;k<dim;k++){
			mag+=(r[i][1][k]-r[i][0][k])*(r[i][1][k]-r[i][0][k]);
		}
		mag=sqrt(mag);  //this is |r0-r1|

		//fb+=kfs*(across/mag-1)*(r[i][0][0]-r[i][1][0]);
		fall[0]+=kfs*(across/mag-1)*(r[i][0][0]-r[i][1][0]);
		//for streching of the filament,
		// -du/dz0=k(across-|r0-r1|)(z0-z1)/|r0-r1|
		
		mag=0;
		for(int k=0;k<dim;k++){
			mag+=(r[i][N-1][k]-r[i][N-2][k])*(r[i][N-1][k]-r[i][N-2][k]);
		}
		mag=sqrt(mag);  //this is |r0-r1|
		//ft+=kfs*(across/mag-1)*(r[i][N-1][0]-r[i][N-2][0]);
		fall[1]+=kfs*(across/mag-1)*(r[i][N-1][0]-r[i][N-2][0]);
		
		//we have dealt with the filament stretching contribution.
		

		
		mag1=0;
		mag2=0;
		dp=0;
		for(int k=0;k<dim;k++){
			mag1+=(r[i][1][k]-r[i][0][k])*(r[i][1][k]-r[i][0][k]);
			mag2+=(r[i][2][k]-r[i][1][k])*(r[i][2][k]-r[i][1][k]);
			dp+=(r[i][2][k]-r[i][1][k])*(r[i][1][k]-r[i][0][k]);
		}
		mag1=sqrt(mag1);
		mag2=sqrt(mag2);
		//computes magnitude and dot product of displacements
		//fb+=kfb/mag1/mag2*(r[i][1][0]-r[i][2][0]+dp*(r[i][1][0]-r[i][0][0])/mag1/mag1);
		fall[0]+=kfb/mag1/mag2*(r[i][1][0]-r[i][2][0]+dp*(r[i][1][0]-r[i][0][0])/mag1/mag1);
		//for bending of the filament at the endpoint,
		// -du/dz0 = k[ (z1-z2) + u0.u1 (z1-z0)/|u0|^2]/|u0||u1|
		
		mag1=0;
		mag2=0;
		dp=0;
		for(int k=0;k<dim;k++){
			mag1+=(r[i][N-2][k]-r[i][N-1][k])*(r[i][N-2][k]-r[i][N-1][k]);
			mag2+=(r[i][N-3][k]-r[i][N-2][k])*(r[i][N-3][k]-r[i][N-2][k]);
			dp+=(r[i][N-3][k]-r[i][N-2][k])*(r[i][N-2][k]-r[i][N-1][k]);
		}
		mag1=sqrt(mag1);
		mag2=sqrt(mag2);
		//ft+=kfb/mag1/mag2*(r[i][N-2][0]-r[i][N-3][0]+dp*(r[i][N-2][0]-r[i][N-1][0])/mag1/mag1);
		fall[1]+=kfb/mag1/mag2*(r[i][N-2][0]-r[i][N-3][0]+dp*(r[i][N-2][0]-r[i][N-1][0])/mag1/mag1);

		
		
		/*double rm7;
		for(int i2=0;i2<M;i2++){//for every other filament in the bundle
			for(int j=1;j<LJrange;j++){//and for every other monomer on those filaments in range
				mag=0;
				for(int k=0;k<dim;k++){
					mag+=(r[i][0][k]-r[i2][j][k])*(r[i][0][k]-r[i2][j][k]);
				}
				mag=sqrt(mag);
				//this is the magnitude of the vector between each filament
				rm7=1/mag;
				rm7=rm7*rm7*rm7;//1/mag^3
				rm7=rm7*rm7;//1/mag^6
				rm7=rm7/mag;
				fb+=eLJ*12*(r[i][0][0]-r[i2][j][0])*rm7;
				// -du/dz0 = 12 eps (z0-z')/|r-r'|^7
			}

			
			
			for(int j=N-2;j>N-1-LJrange;j--){
				mag=0;
				for(int k=0;k<dim;k++){
					mag+=(r[i][N-1][k]-r[i2][j][k])*(r[i][N-1][k]-r[i2][j][k]);
				}
				mag=sqrt(mag);
				//this is the magnitude of the vector between each filament
				rm7=1/mag;
				rm7=rm7*rm7*rm7;//1/mag^3
				rm7=rm7*rm7;//1/mag^6
				rm7=rm7/mag;
				ft+=eLJ*12*(r[i][N-1][0]-r[i2][j][0])*rm7;
				// -du/dz0 = 12 eps (z0-z')/|r-r'|^7
			}
		}
		*/
		//bottom
		for(int b=0;b<links[1].size();b++){
			//for every bond at the next cross=section, we need to compute a force.
			int i1=links[1][b].first;
			int i2=links[1][b].second;
			if((i1==i)|(i2==i)){
				//if the current filament is involved in the bond, we must compute a force
				int iother=i1;
				if(i==iother){
					iother=i2;
				}
				//pick out the other filament.
				
				for(int k=0;k<dim;k++){
					v[k]=0;
				}
				mag1=0;
				mag2=0;
				dp=0;
				for(int k=0;k<dim;k++){
					v[k]=+r[i1][2][k]-r[i1][0][k];
					v[k]=+r[i2][2][k]-r[i2][0][k];
					dp+=v[k]*(r[i][1][k]-r[iother][1][k]);
					mag1+=v[k]*v[k];
					mag2+=(r[i1][1][k]-r[i2][1][k])*(r[i1][1][k]-r[i2][1][k]);
				}//this is the filament direction.
				mag1=sqrt(mag1);  //this is the normalization for the filament
				mag2=sqrt(mag2);
				//fb+=kcs*((r[i][1][0]-r[iother][1][0])-dp*v[0]/mag1/mag1)/mag1/mag2;
				fall[0]+=kcs*((r[i][1][0]-r[iother][1][0])-dp*v[0]/mag1/mag1)/mag1/mag2;
			}
			
		}
		
		
		//top
		for(int b=0;b<links[N-2].size();b++){
			//for every bond at the next cross=section, we need to compute a force.
			int i1=links[N-2][b].first;
			int i2=links[N-2][b].second;
			if((i1==i)|(i2==i)){
				//if the current filament is involved in the bond, we must compute a force
				int iother=i1;
				if(i==iother){
					iother=i2;
				}
				//pick out the other filament.
				
				for(int k=0;k<dim;k++){
					v[k]=0;
				}
				mag1=0;
				mag2=0;
				dp=0;
				for(int k=0;k<dim;k++){
					v[k]=+r[i1][N-1][k]-r[i1][N-3][k];
					v[k]=+r[i2][N-1][k]-r[i2][N-3][k];
					dp+=v[k]*(r[i][N-2][k]-r[iother][N-2][k]);
					mag1+=v[k]*v[k];
					mag2+=(r[i1][N-2][k]-r[i2][N-2][k])*(r[i1][N-2][k]-r[i2][N-2][k]);
				}//this is the filament direction.
				mag1=sqrt(mag1);  //this is the normalization for the filament
				mag2=sqrt(mag2);
				//ft+=kcs*((r[i][N-2][0]-r[iother][N-2][0])-dp*v[0]/mag1/mag1)/mag1/mag2;
				fall[1]+=kcs*((r[i][N-2][0]-r[iother][N-2][0])-dp*v[0]/mag1/mag1)/mag1/mag2;
			}
			
		}
		 
		
	}
	
	
	ftop=fall[1];
	fbottom=fall[0];

	
	
	
	return 0;
}

