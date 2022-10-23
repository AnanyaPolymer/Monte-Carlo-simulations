
#include "header.h"
mt19937 rgen;
uniform_real_distribution<double> uran (0.0, 1.0);
normal_distribution<double> gran (0.0,1.0);
template<typename T>
std::string ctos(T q){
	std::stringstream ss;
	ss<<q;
	return ss.str();
}


extern bool generate_initial_configuration(double H);  //initial_conditions.cpp
extern void generate_crosslinks(double eta,double H);  //get_crosslinks.cpp
extern bool metropolis(double dE);  //metropolis.cpp

extern string energy_file_name(int run,int nstart,int nstop,double H);								//print_data.cpp
extern void print_energy(int nout,ofstream & energyFile,double H);		//print_data.cpp
extern void print_crosslinks(int run,int nstart,int nstop,double H);									//print_data.cpp
extern void print_configurations(int nout,int run,int nstart,int nstop,double H);	
extern string rsq_file_name(int run,int nstart,int nstop,double H);
extern void print_rsq(int nout,ofstream & rsqFile,double H);
extern string dist_file_name(int run,int nstart,int nstop,double H);
extern void print_dist(int nout,ofstream & distFile,double H);
extern string len_file_name(int run,int nstart,int nstop,double H);
extern void print_len(int nout,double *bl,ofstream & lenFile,double H);
extern string force_file_name(int run,int nstart,int nstop,double H);
extern void print_force(int nout,ofstream & forceFile,double H);
extern string y_file_name(int run,int nstart,int nstop,double H);
extern void print_y(int nout,ofstream & yFile,double H);


extern int read_in_data(int run,int nstart,int nstop,double Hread,double H);	//read_data.cpp
extern double compute_filament_energy_diff(int ranchain,int ranmonomer,vector<double> & dx);  //compute_filament_energy.cpp
extern double compute_filament_energy(); //compute_filament_energy.cpp
//extern double compute_LJ_energy_diff(int ranchain,int ranmonomer,vector<double> & dx);  //compute_LJ_energy.cpp
//extern double compute_LJ_energy(); //compute_LJ_energy.cpp
extern double compute_crosslink_energy_diff(int ranchain,int ranmonomer,vector<double> & dx);  //compute_crosslink_energy.cpp
extern double compute_crosslink_energy();  //compute_crosslink_energy.cpp
extern double compute_forces(double *fall);  //forces.cpp



double r[M][N][dim]={0};
//these are the positions of the monomers.
// r[i] refers to chain, r[i][j] refers to jth monomer of chain, r[i][j][k] refers to coordinates.
// r[i][j][0] = x coordinate of the jth monomer in the ith filament
// r[i][j][1] = y coordinate of the jth monomer in the ith filament
// r[i][j][2] = z coordinate of the jth monomer in the ith filament
vector<pair<int,int> > ptemp;  //an empty list of pairs
vector <vector<pair<int,int> > > links(N,ptemp);
		//these are the crosslinks.  This vector has length N.
		//links[i].first = chain1
		//links[i].second =  chain2
		//the length of links is variable, so it's ok to have 0 crosslinks between filaments for the ith monomer
vector <double> dx(dim,0);  //kick we are giving each monomer


double E1=0,E2=0,E3=0;
//energies for filaments, LJ, and crosslinks
double ftop=0,fbottom=0;
//forces on the top and bottom, to compute stress
double nacc=0;
//number of accepted trials
int seed,seedin;

int main(int argc, char* argv[]){

	if (argc<2){
		//std::cerr << "Usage: you must enter the variables in the following order"<<endl;
		//std::cerr << "dim, H, N, M, kfb, nPrint, niter, nSample, nEquil, nrun"<<endl;
		return 0;
	}
	
	int niter = atoi(argv[1]);
	int nstart = 100*(niter-1);
	int nstop = niter*100 - 1;
	
	double H;
	double Hread;

	double eta=1.0;  //this is the fraction of bonds that exist
    double bl[N-1];
    double fall[2];
    bool readIn=true;  //true means "read in data if it exists."  false means "overwrite (delete) data if it exists."
    
	
	
auto t0=chrono::steady_clock::now(); 


for(int comp = 1;comp<2;comp++){
		H  = 1 - 0.01*comp;
		Hread = H + 0.01;

    std::string fldr = "/project/morrison/amondal/pin_data_new_dim"+ctos(dim)+"_kfs"+ctos(kfs)+"_nPrint"+ctos(nPrint)+"_nrun"+ctos(nrun)+"_H"+ctos(H);
    fldr = fldr+"/start"+ctos(nstart)+"_stop"+ctos(nstop);
	std::string mkfldr = "mkdir -p "+fldr;
	system(mkfldr.c_str());


for(int run=nstart;run<=nstop;run++){

		readIn = true;

		seedin = 76;
    	seed = -seedin - ((int)(899037199*niter)) -((int)(89778917*comp)) - ((int)(15678911*run));
    	rgen = mt19937(seed);
		
		std::string dir1 = fldr+ "/configs_nrun"+ctos(run);
		std::string mkdir1 = "mkdir -p "+dir1;
		std::cout<<dir1<<endl;
		system(mkdir1.c_str());	

		std::string dir2 = fldr + "/data_nrun"+ctos(run);
		std::string mkdir2 = "mkdir -p "+dir2;
		//std::cout<<dir2<<endl;
		system(mkdir2.c_str());	



		ofstream energyFile;
		ofstream distFile;
		ofstream lenFile;
		ofstream yFile;
		
		string energyFileName=energy_file_name(run,nstart,nstop,H);  //get the file name for the energy.
		string distFileName = dist_file_name(run,nstart,nstop,H);
		string lenFileName = len_file_name(run,nstart,nstop,H);
		string yFileName = y_file_name(run,nstart,nstop,H);
		

	
	


	int nout=0;  //this is the number of configurations we have printed so far.
	int nstep=-nEquil;  //this is the number of steps we have run so far.
						//a negative value means we will run steps before printing the "first" configuration.

	if(readIn){  //if we want to read in the data, we don't want to create an initial state.  Instead, we want to use the existing data.
		nout=read_in_data(run,nstart,nstop,Hread,H);  // this reads in the configuration and link file for this system.
		//cout<<nout<<"\n";
		if(nout>-1){
			nstep=0;  //we don't equilibrate if we read in data.
			print_configurations(0,run,nstart,nstop,H);
			energyFile.open(energyFileName.c_str(),ios_base::app);  
			distFile.open(distFileName.c_str(),ios_base::app); //prints distance between any two monomers
			lenFile.open(lenFileName.c_str(),ios_base::app); 
			yFile.open(yFileName.c_str(),ios_base::app); 
			print_energy(0,energyFile,H); //open energy file to append to previous file.
			print_y(0,yFile,H); // prints the transverse amplitude at the midpoint
			print_len(0,bl,lenFile,H); // prints the contour length of the chain
			std::cout<<"reading data = "<<nout<<","<<Hread<<endl;
		}
		else{
			nout = 0;
			readIn=false;  //if nout=0, then read_in_data didn't find any data!  We have to generate it.
			std::cout<<"did not find any data!"<<endl;
		}
	}
	if(!readIn){//if we're not reading in data, we need to generate it.
		bool didInitialize=generate_initial_configuration(H);//this constructs the filaments
		if(!didInitialize){
			return 0;  //if we failed to initialize, give up!
		}
		generate_crosslinks(eta,H);  //this connects the different fillaments
        //print_crosslinks(run,nstart,nstop,H);  //this prints the list of connections.
		print_configurations(0,run,nstart,nstop,H);  //this prints the initial configuration.
		energyFile.open(energyFileName.c_str());  //this opens a new energy file.
		distFile.open(distFileName.c_str());
		lenFile.open(lenFileName.c_str());
		yFile.open(yFileName.c_str());
		print_energy(0,energyFile,H);
		print_y(0,yFile,H);
		print_len(0,bl,lenFile,H);

		std::cout<<"generating initial conditions"<<endl;

	}
	
	


	E1=compute_filament_energy();
	//E2=compute_LJ_energy();
	E3=compute_crosslink_energy();
	//compute the different energy contributions from the initial configuration
	
	nout =0;
	while(nout<nPrint){  //we want to print nPrint configurations.  nout tracks how many we printed so far.

		
		
		int ranchain=((int)((uran(rgen)*(M))));  //this is the chain we're changing
		int ranmonomer=((int)((uran(rgen)*(N-pinned_at_bottom-pinned_at_top))))+pinned_at_bottom;  //bead to change.  The minimum selected is pinned_at_bottom, the maximum is N-pinned_at_top.
	
		
		
		
	
		for(int k=0;k<dim;k++){
			dx[k]=dxfactor*aback*gran(rgen);  //this is a normally distributed number with mean 0 and variance aback^2*dxfactor^2
		}
		if((ranmonomer==0)|(ranmonomer==N-1)){
		 	dx[0]=0;  //can't move up or down at the endpoints, but can move laterally unless frozen (using pin_at_top or pin_at_bottom)
		 }
		
		
		
		double dE1=compute_filament_energy_diff(ranchain,ranmonomer,dx);
		//double dE2=compute_LJ_energy_diff(ranchain,ranmonomer,dx);
		double dE3=compute_crosslink_energy_diff(ranchain,ranmonomer,dx);
		//compute the energy differences.
		//double dE=dE1+dE2+dE3;
		double dE=dE1+dE3;
		//total change in energy
		
		if(metropolis(dE)){  //if metropolis criterion says the move should be acepted
			for(int k=0;k<dim;k++){
				r[ranchain][ranmonomer][k]+=dx[k];	//update position
			}
			E1+=dE1;
			//E2+=dE2;
			E3+=dE3;
			//update energy
			nacc+=1;
			//update number of accepted
		}
		
		
		
		if((nstep>=0)&(nstep==nSample)){  //nSample is the number of trials between printings.  If nstep==nSample, it's time to print a configuration.
			nstep=0;//reset nstep to 0
			nout+=1;  //update nout, since we're printing something new.
			
			//print_configurations(nout,run,nstart,nstop,segiter,H);
			//print the bundle positions
			//std::cout<<"printing "<<nout<<endl;
			compute_forces(fall);
			//compute forces at the top and bottom
			print_energy(nout,energyFile,H);
			//print energies and forces
			//print_rsq(nout,rsqFile);
			print_dist(nout,distFile,H);
			print_len(nout,bl,lenFile,H);
			//print_force(nout,forceFile);
			print_y(nout,yFile,H);
			nacc=0;
			
			
			

		}
		nstep+=1;
		//cout<<nout<<"\t"<<nstep<<"\n";
		

		
	

	}
	print_configurations(nout,run,nstart,nstop,H);
	//print_rsq(nout,rsqFile);
	std::cout<<"only printing "<<nout<<endl;
	energyFile.close();  //close the energy output file
	//rsqFile.close();
	distFile.close();
	lenFile.close();
	//forceFile.close();
	yFile.close();
	}
}
  

	auto tf=chrono::steady_clock::now();//compute the final time.
	auto timediff=chrono::duration_cast<chrono::seconds>(tf-t0);
	std::cout<<"program took "<<((double)(timediff.count()))<<" s\n"; 
	
	return 0;
}
		

