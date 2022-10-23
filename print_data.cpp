#include "header.h"

template<typename T>
std::string ctos(T q){
	std::stringstream ss;
	ss<<q;
	return ss.str();
}

extern double r[M][N][dim];
extern vector <vector<pair<int,int> > > links;
extern double E1,E2,E3;
extern double nacc;
extern double ftop,fbottom;

string energy_file_name(int run,int nstart,int nstop,double H);
string rsq_file_name(int run,int nstart,int nstop,double H);
string dist_file_name(int run,int nstart,int nstop,double H);
string len_file_name(int run,int nstart,int nstop,double H);
string y_file_name(int run,int nstart,int nstop,double H);
string force_file_name(int run,int nstart,int nstop,double H);
void print_crosslinks(int run,int nstart,int nstop,double H);
void print_energy(int nout,ofstream & energyFile,double H);
void print_configurations(int nout,int run,int nstart,int nstop,double H);
void print_rsq(int nout,ofstream & rsqFile,double H);
void print_dist(int nout,ofstream & distFile,double H);
void print_len(int nout,double *bl,ofstream & lenFile,double H);
void print_y(int nout,ofstream & yFile,double H);
void print_force(int nout,ofstream & forceFile,double H);


string energy_file_name(int run,int nstart,int nstop,double H){
	string maindir = "/project/morrison/amondal/pin_data_new_dim"+ctos(dim)+"_kfs"+ctos(kfs)+"_nPrint"+ctos(nPrint)+"_nrun"+ctos(nrun)+"_H"+ctos(H);
    maindir = maindir+"/start"+ctos(nstart)+"_stop"+ctos(nstop);
	string ename= maindir + "/energy_acc_force_M"+ctos(M)+"_N"+ctos(N)+"_nrun"+ctos(run)+".txt";  //this is the energy file name.
	return ename;
}

void print_energy(int nout,ofstream & energyFile,double H){
	energyFile<<nout<<","<<E1<<","<<E3<<","<<100.0*nacc/((double)nSample)<<","<<ftop<<","<<fbottom<<endl;  //this prints the energies, the acceptance rate, and the forces.
}


void print_crosslinks(int run,int nstart,int nstop,double H){
	
	//let's print the crosslinks!
	string maindir = "/project/morrison/amondal/pin_data_new_dim"+ctos(dim)+"_kfs"+ctos(kfs)+"_nPrint"+ctos(nPrint)+"_nrun"+ctos(nrun)+"_H"+ctos(H);
    maindir = maindir+"/start"+ctos(nstart)+"_stop"+ctos(nstop);
	ofstream linkfile;
	string linkname= maindir + "/links_M"+ctos(M)+"_N"+ctos(N)+"_nrun"+ctos(run)+".txt";  //filename for the links
	linkfile.open(linkname.c_str());  //open the file
	for(int i=0;i<links.size();i++){
		for(int j=0;j<links[i].size();j++){
			linkfile<<links[i][j].first<<","<<links[i][j].second<<","<<i<<endl;  //this prints all links.  The output format is Filament1,Filament2,monomer, meaning the two filaments are connected at that monomer.
		}
	}
	linkfile.close();//close the file
	
}


void print_configurations(int nout,int run,int nstart,int nstop,double H){
	string maindir = "/project/morrison/amondal/pin_data_new_dim"+ctos(dim)+"_kfs"+ctos(kfs)+"_nPrint"+ctos(nPrint)+"_nrun"+ctos(nrun)+"_H"+ctos(H);
    maindir = maindir+"/start"+ctos(nstart)+"_stop"+ctos(nstop);
	string dir = maindir + "/configs_nrun"+ctos(run)+"/";

	ofstream outfile;

	string outname= dir+"configuration_M"+ctos(M)+"_N"+ctos(N)+"_out"+ctos(nout)+".txt";  //this is the file name
	outfile.open(outname.c_str());  //open the file
	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			outfile<<r[i][j][0];
			for(int k=1;k<dim;k++){
				outfile<<","<<r[i][j][k];  //print the coordinates
			}
			outfile<<"\n";
		}			
	}
	outfile.close();  //close the file
	//std::cout<<outname<<endl;
}

string rsq_file_name(int run,int nstart,int nstop,double H){
	string maindir = "/project/morrison/amondal/pin_data_new_dim"+ctos(dim)+"_kfs"+ctos(kfs)+"_nPrint"+ctos(nPrint)+"_nrun"+ctos(nrun)+"_H"+ctos(H);
    maindir = maindir+"/start"+ctos(nstart)+"_stop"+ctos(nstop);
	string dir = maindir + "/data_nrun"+ctos(run)+"/";
	string rsqname= dir + "rsq_M"+ctos(M)+"_N"+ctos(N)+"_nrun"+ctos(run)+".txt";  //this is the energy file name.
	return rsqname;
}

void print_rsq(int nout,ofstream & rsqFile,double H){
	double rsq = 0.0;
	for(int k=0;k<dim;k++){
			rsq += (r[0][N-1][k] - r[0][0][k])*(r[0][N-1][k] - r[0][0][k]);
			}
	rsqFile<<nout<<","<<rsq<<endl;
}

string dist_file_name(int run,int nstart,int nstop,double H){
	string maindir = "/project/morrison/amondal/pin_data_new_dim"+ctos(dim)+"_kfs"+ctos(kfs)+"_nPrint"+ctos(nPrint)+"_nrun"+ctos(nrun)+"_H"+ctos(H);
    maindir = maindir+"/start"+ctos(nstart)+"_stop"+ctos(nstop);
	string dir = maindir + "/data_nrun"+ctos(run)+"/";
	string distname= dir + "dist_M"+ctos(M)+"_N"+ctos(N)+"_nrun"+ctos(run)+".txt";  //this is the energy file name.
	return distname;
}

void print_dist(int nout,ofstream & distFile,double H){
	double diff = 0.0;
	for(int k =0;k<dim;k++){
		diff += (r[0][74][k] - r[0][24][k])*(r[0][74][k] - r[0][24][k]);
	}
	distFile<<nout<<","<<diff<<endl;
}

string len_file_name(int run,int nstart,int nstop,double H){
	string maindir = "/project/morrison/amondal/pin_data_new_dim"+ctos(dim)+"_kfs"+ctos(kfs)+"_nPrint"+ctos(nPrint)+"_nrun"+ctos(nrun)+"_H"+ctos(H);
    maindir = maindir+"/start"+ctos(nstart)+"_stop"+ctos(nstop);
	string dir = maindir + "/data_nrun"+ctos(run)+"/";
	string lenname= dir + "len_M"+ctos(M)+"_N"+ctos(N)+"_nrun"+ctos(run)+".txt";  //this is the energy file name.
	return lenname;
}

void print_len(int nout,double *bl,ofstream & lenFile,double H){
	for(int i =0;i<N-1;i++){
		bl[i]=0.0;
	}
	for(int k=0;k<dim;k++){
		bl[0] += (r[0][1][k] - r[0][0][k])*(r[0][1][k] - r[0][0][k]);
		bl[N-2] += (r[0][N-1][k] - r[0][N-2][k])*(r[0][N-1][k] - r[0][N-2][k]);
	}
	bl[0] = sqrt(bl[0]);
	bl[N-2] = sqrt(bl[N-2]);

	
	for(int i=1;i<N-2;i++){
		for(int k=0;k<dim;k++){
			bl[i] += (r[0][i+1][k] - r[0][i][k])*(r[0][i+1][k] - r[0][i][k]);
			}
			bl[i] = sqrt(bl[i]);
	}
	double len =0.0;
	for(int i =0;i<N-1;i++){
		len += bl[i]; 
		}
		lenFile<<nout<<","<<len<<","<<r[0][49][1]/len<<endl;
}

string y_file_name(int run,int nstart,int nstop,double H){
	string maindir = "/project/morrison/amondal/pin_data_new_dim"+ctos(dim)+"_kfs"+ctos(kfs)+"_nPrint"+ctos(nPrint)+"_nrun"+ctos(nrun)+"_H"+ctos(H);
    maindir = maindir+"/start"+ctos(nstart)+"_stop"+ctos(nstop);
	string dir = maindir + "/data_nrun"+ctos(run)+"/";
	string yname= dir + "y_M"+ctos(M)+"_N"+ctos(N)+"_nrun"+ctos(run)+".txt";  //this is the energy file name.
	return yname;
}

void print_y(int nout,ofstream & yFile,double H){
	yFile<<nout<<","<<r[0][0][0]<<","<<r[0][24][0]<<","<<r[0][49][0]<<","<<r[0][74][0]<<","<<r[0][99][0]<<","<<r[0][0][1]<<","<<r[0][24][1]<<","<<r[0][49][1]<<","<<r[0][74][1]<<","<<r[0][99][1]<<","<<r[0][0][2]<<","<<r[0][24][2]<<","<<r[0][49][2]<<","<<r[0][74][2]<<","<<r[0][99][2]<<endl;
}

//Forces on the chain
string force_file_name(int run,int nstart,int nstop,double H){
	string maindir = "/project/morrison/amondal/pin_data_new_dim"+ctos(dim)+"_kfs"+ctos(kfs)+"_nPrint"+ctos(nPrint)+"_nrun"+ctos(nrun)+"_H"+ctos(H);
    maindir = maindir+"/start"+ctos(nstart)+"_stop"+ctos(nstop);
	string dir = maindir + "/data_nrun"+ctos(run)+"/";
	string forcename= dir + "force_M"+ctos(M)+"_N"+ctos(N)+"_nrun"+ctos(run)+".txt";  //this is the energy file name.
	return forcename;
}

void print_force(int nout,ofstream & forceFile,double H){
	forceFile<<nout<<","<<fbottom<<","<<ftop<<endl;
}

