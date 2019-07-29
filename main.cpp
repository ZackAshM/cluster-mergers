/* -------------------------------
   N-body 2D Star Cluster Merger Simulation
   Zachary Martin
   Physics 305
   Begin Date: 10 April 2019
   To compile: g++ cluster_generator.cpp VecFRK4.cpp main.cpp -o ~/bin/main
   requires the files "Vector3D.h" and "QuadTree.h"
   ------------------------------- */

using namespace std;
#include <iostream>
#include <iomanip>
#include <fstream>
#define USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <vector>

// for star cluster generation
#define seed0 4102019213
// #define Nmax 10000		// default = 10000; use small Nmax for testing
// Star cluster properties:
#define Mcl 1.0e3		// solar masses; typical open cluster mass
#define Rcl 15.0		// ly; typical open cluster radius
//#define Mcl2 2.0e2		// solar masses; typical open cluster mass
//#define Rcl2 5.0		// ly; typical open cluster radius
#define Mcl2 Mcl		// solar masses; typical open cluster mass
#define Rcl2 Rcl		// ly; typical open cluster radius
//#define Rcl 1.5813e-5
#define clust_vel 0.0003	// ly/(100 yr) total collision velocity of cluster
// #define shape1 1		// for cluster_gen(); box=0, sphere=1, disk=2
// #define shape2 1

#define notes false		// turn notes on/off

// this includes spatial vector definitions and operations
#include "Vector3D.h"
// this includes Barnes-Hut QuadTree methods for N-body calculation
// this also includes the Body struct
#define scale 4.5*Rcl		// ly, length of tree's root node
#define G 1.56700e-9		// ly^3 Msun^-1 (100 yr)^-2, gravitational constant
#include "QuadTreeClass.h"


vector<Vector3D> forces;
int seed;

/* This function is found in an external file: cluster_generator.cpp
   It should generate a list of N<Nmax 7D arrays (listed using a vector).
   The amount of stars, N, is generated via cluster mass Mcl.
   To extrapolate the postion and masses of each star generated:
   star position Vector3D's: 0th-2nd elements of the index'th 7D array, cluster_gen(index,...)[0,1,2]
   star masses: 3rd element of the index'th 7D array, cluster_gen(index,...)[3]
   star velocity Vector3D's: 4th-7th elements of the index'th 7D array, cluster_gen(index,...)[4,5,6]
   extract number of stars using cluster_gen(Nmax,...)[0] = N
*/
vector<double> cluster_gen(int, int, int, double, double, int);

/* This function is found in an external file: VecFRK4.cpp
   It should calculate the 4th-order Runge-Kutta (RK4) correction
   terms in solving the differential equations of motion.
*/
Vector3D VecFRK4( int , Vector3D (*f_x)(int, double, Vector3D, Vector3D),
		  Vector3D (*f_v)(int, double, Vector3D, Vector3D),
		  int , double , Vector3D, Vector3D , double );
		  
Vector3D f_r(int, double, Vector3D, Vector3D);
Vector3D f_v(int, double, Vector3D, Vector3D);


/* --------------------------- BEGIN MAIN ------------------------------- */
int main(int argc, char *argv[])
{

/* ------------------ ARGUMENTS ------------------- */

    // argument check
    if(argc<8){
	cerr << "usage: main [seed (0 for default)][Star Nmax][shape 1][shape 2][cluster outfile][Tmax][main outfile]" << endl;
	cerr << "For shape: box=0, sphere=1, disk=2" << endl;
	exit(0);
    }
	
    cerr << "Initializing Arguments..." << endl;
	
    if(atoi(argv[1]) == 0){			// seed
	seed = seed0;
    }else seed = atoi(argv[1]);
	
    int Nmax = atoi(argv[2]);		// max no. of stars
    int shape1 = atoi(argv[3]);		// shape of first cluster
    int shape2 = atoi(argv[4]);		// shape of second cluster
	
    ofstream clustfile;
    clustfile.open(argv[5]);		// file with initial clusters
	
    int Tmax = atof(argv[6]);		// max time of collision
    double dt = 1.0;
	
    ofstream outfile;
    outfile.open(argv[7]);			// main outfile for time evolution
	
/* ------------------------------------------------ */
	
/* --------------- CLUSTER GENERATION ------------- */

    cerr << "Beginning Cluster Generation..." << endl;

    // First cluster
    vector<Vector3D> pos_star1(Nmax);
    vector<Vector3D> vel_star1(Nmax);
    vector<double> mass_star1(Nmax);
	
    cerr << "Generating first cluster..." << endl;
    // extract number of stars generated
    int N1 = min((int)cluster_gen(Nmax, Nmax, seed, Mcl, Rcl, shape1)[0], Nmax);
    cerr << "Number of Stars (1) = " << N1 << endl;

    clustfile << "# (1) N  STAR POSITIONS		STAR VELOCITIES 		MASSES" << endl;
    for(int i=0;i<N1;i++){
	pos_star1[i] = Vector3D( cluster_gen(i, Nmax, seed, Mcl, Rcl, shape1)[0] - Rcl, 
				 cluster_gen(i, Nmax, seed, Mcl, Rcl, shape1)[1], 
				 0.0 /*cluster_gen(i, Nmax, seed, Mcl, Rcl, shape1)[2]*/);
	vel_star1[i] = Vector3D( cluster_gen(i, Nmax, seed, Mcl, Rcl, shape1)[4], 
				 cluster_gen(i, Nmax, seed, Mcl, Rcl, shape1)[5], 
				 0.0 /*cluster_gen(i, Nmax, seed, Mcl, Rcl, shape1)[6]*/);
	mass_star1[i] = cluster_gen(i, Nmax, seed, Mcl, Rcl, shape1)[3];
	
	if((i+1)%(100) == 0){
	    cerr << "Creating star # " << i+1 << "..." << endl;
	}
	clustfile << i << "  " << pos_star1[i].Getx() << "  " << pos_star1[i].Gety() << "  " 
		  << pos_star1[i].Getz() << "  " << vel_star1[i].Getx() << "  " 
		  << vel_star1[i].Gety() << "  " << vel_star1[i].Getz() << "  " 
		  << mass_star1[i] << endl;
    }
	
    clustfile << endl;
	
    // Second cluster
    vector<Vector3D> pos_star2(Nmax);
    vector<Vector3D> vel_star2(Nmax);
    vector<double> mass_star2(Nmax);
	
    cerr << "Generating second cluster..." << endl;
    // extract number of stars generated
    int N2 = min((int)cluster_gen(Nmax, Nmax, seed*7, Mcl2, Rcl2, shape2)[0], Nmax);	
    cerr << "Number of Stars (2) = " << N2 << endl;

    clustfile << "# (2) N  STAR POSITIONS		STAR VELOCITIES 		MASSES" << endl;
    for(int i=0;i<N2;i++){
	pos_star2[i] = Vector3D( cluster_gen(i, Nmax, seed*7, Mcl2, Rcl2, shape2)[0] + Rcl, 
				 cluster_gen(i, Nmax, seed*7, Mcl2, Rcl2, shape2)[1], 
				 0.0 /*cluster_gen(i, Nmax, seed*7, Mcl2, Rcl2, shape2)[2]*/ );
	vel_star2[i] = Vector3D( cluster_gen(i, Nmax, seed*7, Mcl2, Rcl2, shape2)[4], 
				 cluster_gen(i, Nmax, seed*7, Mcl2, Rcl2, shape2)[5], 
				 0.0 /*cluster_gen(i, Nmax, seed*7, Mcl2, Rcl2, shape2)[6]*/);
	mass_star2[i] = cluster_gen(i, Nmax, seed*7, Mcl2, Rcl2, shape2)[3];
	
	if((i+1)%(100) == 0){
	    cerr << "Creating star # " << i+1 << "..." << endl;
	}
	clustfile << i << "  " << pos_star2[i].Getx() << "  " << pos_star2[i].Gety() << "  " 
		  << pos_star2[i].Getz() << "  " << vel_star2[i].Getx() << "  " 
		  << vel_star2[i].Gety() << "  " << vel_star2[i].Getz() << "  " 
		  << mass_star2[i] << endl;
    }
	
    int Ntot = N1 + N2;		// total number of stars and also body count
	
    clustfile.close();
	
    cerr << "Star Cluster Generation Finished." << endl;
    cerr << "Assigning Bodies for QuadTree..." << endl;

    // assign bodies
    
//    double bigstarmass = 150.;
//    Vector3D bigstarpos = (0.5*Rcl, 0.0, 0.0);
//    Vector3D bigstarvel = (-clust_vel, 0.0, 0.0);
    
    vector<Body> bodies;
    
//    Body bigstar(bigstarmass, bigstarpos, bigstarvel);
//    bodies.push_back(bigstar);
    
    for(int i=0; i<N1; i++){
	Body temp(mass_star1[i], pos_star1[i], vel_star1[i] + Vector3D(clust_vel, 0.0, 0.0));
	bodies.push_back(temp);
    }
    for(int i=0; i<N2; i++){
	Body temp(mass_star2[i], pos_star2[i], vel_star2[i] - Vector3D(clust_vel, 0.0, 0.0));
	bodies.push_back(temp);
    }
	
    cerr << "Bodies assigned." << endl;
	
/* ---------- END OF CLUSTER GENERATION ------------ */
	
    cerr << "Beginning time evolution..." << endl;
	
    // this is the primary loop
    int ii = 0;						// index for printing every x yr
    for (double t = 0; t < Tmax; t += dt) {
	static int div = (int)(100/dt);		// print only every x yr
		
	// create a new QuadTree
	if(t==0) cerr << "Creating 1st QuadTree..." << endl;	// checking if this works the first time
	QuadTreeNode* tree = new QuadTreeNode();
	tree->root_node = true;

    // reserve enough memory
    forces.reserve(bodies.size());

	// insert all of the masses into the tree
	if(t==0) cerr << "Inserting Bodies into 1st QuadTree..." << endl;
	for (int i = 0; i < bodies.size(); i++) {
	    if(t==0){
		if(i%50 == 0){
		    cerr << "Inserting " << i << "th body..." << endl;
		}
	    }
	    tree->insert(bodies[i]);
	}
		
	if(t==0) cerr << "1st Body insertion finished." << endl;
	if(t==0) cerr << "Beginning RK4 Calculation..." << endl;
	
	if(ii%div == 0){
		outfile << t+dt << " ";
	}
			
	for(int j = 0; j < bodies.size(); j++){
	    Vector3D rtold, vtold, rt, vt;
		
		
		
	    // compute force on each body
//			vector<Vector3D> forces;	// declared globally for f_v
	    for (int i = 0; i < bodies.size(); i++) {
	    	Vector3D dvdt = (1./bodies[i].mass)*tree->netforce_on(bodies[i]);
			forces.push_back(dvdt);
			if(notes) cerr << "Output accel = " << forces[i].printt() << endl;
	    }
		
	    /* These lines update the position and velocity of jth body */
	    rtold = bodies[j].position;
	    vtold = bodies[j].velocity;
	    rt = rtold + VecFRK4(0,f_r,f_v,j,t,rtold,vtold,dt);
	    vt = vtold + VecFRK4(1,f_r,f_v,j,t,rtold,vtold,dt);

	    /* set the old values to the new updated values for next iteration */
	    bodies[j].position = rt;
	    bodies[j].velocity = vt;
			
		if(ii%div == 0){
	    	outfile << rt.Getx() << " " << rt.Gety() << " " << vt.Getx() << " " << vt.Gety() << " " << bodies[j].mass << " ";
		}
			
    	}
    	
   		if(ii%div == 0){
			outfile << endl;
		}
    	
    	if(ii%(10000) == 0){		// note every 50 years
	    cerr << "Calculating... t = " << t+dt << " ..." << endl;
	}
	
//		if(ii<6){
//			tree->print();
//		}
	
    	ii++;
    	
    	// these lines delete and clear values no longer needed
    	// in order to make room in memory
    	forces.clear();
    	tree->free();
    	delete tree;
    
    }


}

/* ---------------------------- END OF MAIN ------------------------------- */

Vector3D f_r(int i, double t, Vector3D r, Vector3D v)
{
    return(v);
}

Vector3D f_v(int i, double t, Vector3D r, Vector3D v)
{
    return(forces[i]);
}


