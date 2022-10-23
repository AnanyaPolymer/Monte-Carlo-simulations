#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <map>
#include <utility>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <random>
#include <chrono>
using namespace std;


 #define dim 3  //three dimensional system
 #define N 100 // number of monomers per filament.  All filaments equal length
 //#define H 0.99

#define RING_GEOMETRY	//this specifies a hollow ring
#define FIL_IN_RING 1	// use 1 or 2 for the current project

//#define HEXAGON_GEOMETRY//this specifies a hexagonal cross section.
//#define NUM_LAYERS 2



#ifdef RING_GEOMETRY
#define M (FIL_IN_RING)    //number of filaments in the system
#endif
#ifdef HEXAGON_GEOMETRY
#define M (1+6*(NUM_LAYERS*(NUM_LAYERS+1)/2))    //number of filaments in the system.  NUM_LAYERS=0 means just one filament, at the center.
#endif







#define pinned_at_top 1
#define pinned_at_bottom 1
//this is the number of COMPLETELY immobile beads at the top or bottom, beyond the constrained bead.  Use 1 if you wish to constrain uz=+-1 at the walls.

 #define nPrint 5000	//number of configurations to actually generate
 #define nSample 500000	//number of reorganization steps between configurations
 #define nEquil 10000000		//number of samples go generate before first print
 #define nrun 50000 //100 runs in per batch. Total 500 batches



#define dxfactor 0.25  //variance of a trial


#define kfs 10		//strength of the backbone bonds
#define kcs 0		//strength of the crosslink bond
#define aback 7			//length of a backbone bond
#define across 12		//length of a crosslink bond

#define kfb 99.5  //resistance to bending
#define kcb 0  //could also use 234 as the weaker link per Wang paper

//#define eLJ 4			//strength of LJ repulsion
//#define LJrange 10		//number of bonds to search for the interaction.


//#define pi 3.141592653589
#define pi acos(-1)