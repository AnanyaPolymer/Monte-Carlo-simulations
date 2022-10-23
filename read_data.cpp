#include "header.h"


template<typename T>
std::string ctos(T q){
	std::stringstream ss;
	ss<<q;
	return ss.str();
}

extern double r[M][N][dim];
extern vector< vector <pair<int,int> > > links;


int read_in_data(int run,int nstart,int nstop,double Hread,double H){

	string dir = "/project/morrison/amondal/pin_data_new_dim"+ctos(dim)+"_kfs"+ctos(kfs)+"_nPrint"+ctos(nPrint)+"_nrun"+ctos(nrun)+"_H"+ctos(Hread);
    dir = dir+"/start"+ctos(nstart)+"_stop"+ctos(nstop);
	string filelisttxt="fns_start"+ctos(nstart)+"_stop"+ctos(nstop)+".txt";
	string ls="ls "+dir+"/configs_nrun"+ctos(run)+" > "+filelisttxt;
	system(ls.c_str());
	ifstream fnlist;
	fnlist.open(filelisttxt.c_str());
	string fn;
	int outmax=-1;
	while(getline(fnlist,fn)){
		stringstream ss(fn);
		//the file names have the form
		//  configuration_M[]_N[]_out[].txt.
		//we want to choose the largest out value correponding
		//to our M,N values
		string data;
		getline(ss,data,'_');  //here we're splitting on "_"
								//data= "configuration".  we don't care.
		getline(ss,data,'_');  //data = "M[]".  We do care.
		data.erase(data.begin());
		int m=atoi(data.c_str());  //m is value of M for this file.
		getline(ss,data,'_');  //data = "N[]".  We do care.
		data.erase(data.begin());
		int n=atoi(data.c_str());  //n is value of N for this file.
		
		getline(ss,data,'_');  //data = "out[].txt".  We do care.
		data=data.substr(3,data.length()-4);  //drop the .txt;
		int out=atoi(data.c_str());  //n is value of N for this file.

		if((m==M)&(n==N)){
			if(out>outmax){
				outmax=out;
			}
		}
	}
	fnlist.close();
	
	string delstr="rm "+filelisttxt;
	system(delstr.c_str());
	ifstream restartfile;
	string restartname=dir+"/configs_nrun"+ctos(run)+"/configuration_M"+ctos(M)+"_N"+ctos(N)+"_out"+ctos(outmax)+".txt";
	restartfile.open(restartname.c_str());
	string restartdata;
	int i=0;
	int j=0;
	double sfac = H/Hread;
	while(getline(restartfile,restartdata)){
		stringstream ss(restartdata);
		string x,y,z;
		for(int k=0;k<dim;k++){
			getline(ss,x,',');  //get xyz components;
			r[i][j][k]=atof(x.c_str());
		}
		r[i][j][0] = r[i][j][0]*(sfac);
		j+=1;
		if(j==N){
			j=0;
			i+=1;
		}
	}
	restartfile.close();

	

	

	string dir1 = dir + "/data_nrun"+ctos(run);
	restartname=dir1+"/links_M"+ctos(M)+"_N"+ctos(N)+"_nrun"+ctos(run)+".txt";
	restartfile.open(restartname.c_str());
	string lk;
	while(getline(restartfile,lk)){
		stringstream ss(lk);
		string f1,f2,mon;
		getline(ss,f1,',');  //get first filament;
		getline(ss,f2,',');  //get second filament;
		getline(ss,mon,',');  //get monomer;
		pair<int,int> p=make_pair(atoi(f1.c_str()),atoi(f2.c_str()));
		links[atoi(mon.c_str())].push_back(p);
	}
	restartfile.close();
	
	
		
	return outmax;

}
