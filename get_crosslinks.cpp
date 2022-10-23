#include "header.h"
extern mt19937 rgen; //if rgen is changing this line is changing too
extern uniform_real_distribution<double> uran; //if uran is changing this                                                  line must change too



extern double r[M][N][dim];
extern vector< vector <pair<int,int> > > links;


void generate_crosslinks(double eta,double H){
    
    
    
    for(int i1=0;i1<M-1;i1++){
        for(int i2=i1+1;i2<M;i2++){
            double sep=0;
            for(int k=0;k<dim;k++){
                sep+=(r[i1][0][k]-r[i2][0][k])*(r[i1][0][k]-r[i2][0][k]);
            }
            sep=sqrt(sep);
            if(sep<across*1.00001){  //the small factor avoids rounding errors
                //what we now know is that the filaments i1 and i2 are separated by less than a crosslink bond length.  Because of that, we should try to link them.
    
                pair<int,int> p=make_pair(i1,i2);
                
                for(int j=0;j<N;j++){
                    if(uran(rgen)<eta){
                        links[j].push_back(p);
                        //with probability eta, we link the two filaments in this cross-section.
                    }
                    
                }
            }
        }
    }
}

