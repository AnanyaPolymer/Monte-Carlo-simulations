#include "header.h"


template<typename T>
std::string ctos(T q){
	std::stringstream ss;
	ss<<q;
	return ss.str();
}

extern double r[M][N][dim];


bool generate_initial_configuration(double H){
	
	
	
#ifdef RING_GEOMETRY
		
		
		if(M==1){
			for(int i=0;i<N;i++){
				r[0][i][0]=((double)i)*H*aback;
				for(int j=1;j<dim;j++){
					r[0][i][j]=0;
				}
			}
			return true;
		}
		else if(M==2){
			for(int i=0;i<N;i++){
				r[0][i][0]=((double)i)*H*aback;
				r[1][i][0]=((double)i)*H*aback;
				r[0][i][1]=-across/2;
				r[1][i][1]=across/2;
				for(int j=2;j<dim;j++){
					r[0][i][j]=0;
					r[1][i][j]=0;

				}
			}
			return true;
		}
		else{
			std::cout<<"M>=3 disabled!";
			return false;
		}
		/*
		else{		//if M=2, we just have a pair of filaments so put them next to each other
			for(int i=0;i<M;i++){  //i is the ith filament
				double c=cos(((double)i)/((double)M)*2*pi);  //draw a circle for x
				double s=sin(((double)i)/((double)M)*2*pi);  //draw a circle for y
				double Rcircle=across/(2*sin(pi/M));
				for(int j=0;j<N;j++){  // j is the jth monomer in filament i
					r[i][j][2]=c*Rcircle;  //pick the radius across/(2sin(pi/M))
					r[i][j][1]=s*Rcircle;
					r[i][j][0]=((double)j)*H*aback;
				}
			}
		}
		return true;
		 */
	
#endif
	
	
	
	
#ifdef HEXAGON_GEOMETRY
	/*
	int filCnt=0;
	//this block builds a hexagonal lattice.  
	for(int j=-NUM_LAYERS;j<=NUM_LAYERS;j++){
		int kmin=-NUM_LAYERS; //
		if(j<0){
			kmin-=j; //kmin = kmin - j
		}
		int kmax=NUM_LAYERS;
		if(j>0){
			kmax-=j;
		}
		for(int k=kmin;k<=kmax;k++){
			r[filCnt][0][0]=((double)(j+2*k))*across/2;
			r[filCnt][0][1]=sqrt(3)*((double)j)*across/2;
			r[filCnt][0][2]=0;
			filCnt+=1;
		}
	}
	//this generates a hexagonal lattice.  The particular filament index is not meanginful (e.g. the 0th filament is not in the center).  links will still be generated between all filaments that are sufficiently close.
	
	for(int i=0;i<M;i++){  //i is the filament index
		for(int j=1;j<N;j++){  //j is the monomer index
			r[i][j][0]=r[i][0][0];  //x component is held fixed
			r[i][j][1]=r[i][0][1];	//y component is held fixed
			r[i][j][2]=r[i][j-1][2]+H*aback;  //z component is the previous, + H*aback
		}
	}
	
	return true;
	 */
	std::cout<<"HEXAGON DISABLED!\n";
	return false;
#endif
	
	cout<<"you have not specified a geometry.  Initial conditions are not set up, and the program will halt now!\n";
	return false;
	
	
}

