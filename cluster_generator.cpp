/* -------------------------------
N-body Star Cluster Merger Simulation:
	Cluster Generating Function
Zachary Martin
Physics 305
------------------------------- */

using namespace std;
#include <iostream>
#define USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <vector>
#include "Vector3D.h"

/* typical defines to use for this function
#define Nmax 10000		// avoid very large N for calculations
#define Mcl 5.0e2		// solar masses
#define Rcl 30.0		// ly, radius of cluster
#define total_vel 1.0e-11	// ly/s total velocity of cluster
*/
// #define G 6.67e-11		// gravitational constant
#define G 1.56700e-9		// ly^3 Msun^-1 (100 yr)^-2, gravitational constant


// Here is the primary function

vector<double> cluster_gen(int index, int Nmax, int seed, double Mcl, double Rcl, int shape)

/*	index calls the index'th 7D array
	Nmax is the upper limit for amount of stars
	seed is the random generator seed for generating random masses/positions
	Mcl and Rcl indicate the total cluster mass and size
	shape=0,1,2 indicates the cluster shape: 0-box, 1-sphere, and 2-disk ~10ly thick
	
	This function should generate a list of N<Nmax 7D arrays (listed using a vector).
	  The amount of stars, N, is generated via cluster mass Mcl.
	  To extrapolate the postion, velocity and masses of each star generated:
	    star position Vector3D's: 0th-2nd elements of the index'th 7D array
	    star masses: 3rd element of the index'th 7D array
	    star velocity Vector3D's: 4th-7th elements of the index'th 7D array
	    extract number of stars using cluster_gen(Nmax,...)[0] = N
	    
	EXAMPLE OF USAGE:
	vector<Vector3D> pos_star(Nmax);
	vector<Vector3D> vel_star(Nmax);
	vector<double> mass_star(Nmax);
	for(int i=0;i<N;i++){
		pos_star[i] = Vector3D(cluster_gen(i,...)[0], 
			      cluster_gen(i,...)[1], cluster_gen(i,...)[2]);
		pos_star[i] = Vector3D(cluster_gen(i,...)[4], 
			      cluster_gen(i,...)[5], cluster_gen(i,...)[6]);
		mass_star[i] = cluster_gen(i,...)[3];
	outfile << pos_star[i].Getx() << "  " << pos_star[i].Gety() << "  " 
		<< pos_star[i].Getz() << "  " << vel_star[i].Getx() << "  " 
		<< vel_star[i].Gety() << "  " << vel_star[i].Getz() << "  " 
		<< mass_star[i] << endl;
	}
*/


{
	srand48(seed);			// seed for random number generator
	
	double *x, *y, *z, *vx, *vy, *vz, *mass;
	
	
	// create masses based on Mcl; total mass ~Mcl
	mass = new double [Nmax]();			// array of masses
	int j = 0;					// (star) counter
	double Mtot = 0.0;				// total mass; to compare to Mcl
	while(Mtot<Mcl){
		if(j==Nmax){break;}			// stop if max number reached
		mass[j] = 9.5*drand48() + 0.5;		// typical cluster star masses 
							// range from 0.5-10Msun
		Mtot += mass[j];			// accumulate total mass of cluster
		j++;					// records number of stars made
	}
	
	int N = j;		// setting number of stars based on mass
	
	if(index==Nmax){
		return vector<double> {N*1.0};		// used to extract number of stars generated
	}
	
	
	x = new double [N]();		// list of x values
	y = new double [N]();		// list of y values
	z = new double [N]();		// list of z values
	
	// list of velocity values
	double orb_vel = pow(0.5*G*Mtot/Rcl,0.5);
	vx = new double [N]();	
	vy = new double [N]();
	vz = new double [N]();
	
	// create position elements 0-2 for ith array (x, y, z)
	//  this will depend on set shape of the distribution:
	//  0: box
	//  1: sphere
	//  2: disk ~10ly thick
	//  note: cluster is generated with respect to origin,
	//  need to displace the 2nd cluster when extrapolating position vector
	if(shape == 0){	// box side length 2*Rcl
		for(int i=0;i<N;i++){
			x[i] = 2.0*Rcl*drand48() - Rcl;
			y[i] = 2.0*Rcl*drand48() - Rcl;
			z[i] = 2.0*Rcl*drand48() - Rcl;
	}} else if(shape == 1){		// sphere
		for(int i=0;i<N;i++){
			double theta = 2.0*M_PI*drand48();		// Correct distribution of angles (cred: wolfram)
        	double phi = acos(2.0*drand48() - 1.0);
        	double scale = drand48()*Rcl;
        	
			x[i] = scale*sin(phi)*cos(theta);
			y[i] = scale*sin(phi)*sin(theta);
			z[i] = scale*cos(phi);
			
			vx[i] = -(orb_vel)*sin(theta);
			vy[i] = (orb_vel)*cos(theta);
			vz[i] = 0.0;
			
	}} else if(shape == 2){		// disk ~10ly thick
		for(int i=0;i<N;i++){
			double theta = 2.0*M_PI*drand48();
			double thickness = Rcl/3.0;
			double scale = drand48()*Rcl;
			x[i] = scale*cos(theta);
			y[i] = scale*sin(theta);
			z[i] = 2.0*thickness*drand48() - thickness;
	}} else {
		cerr << "Error: no cluster shape specified. Terminating program..." << endl;
		exit(0);
	}
	
	// now we create the arrays based on the generated values
	vector<vector<double>> clust_vect(N);		// list of arrays (each called clust_list)
	vector<double> clust_list(7);			// 7D arrays contained in clust_vect
	
	// assign elements into each array
	for(int k=0;k<N;k++){
		clust_list[0] = x[k];
		clust_list[1] = y[k];
		clust_list[2] = z[k];
		clust_list[3] = mass[k];
		clust_list[4] = vx[k];
		clust_list[5] = vy[k];
		clust_list[6] = vz[k];
	
	// assign arrays into vector list
		clust_vect[k] = clust_list;
	}
	
	delete x, y, z, vx, vy, vz, mass;	// delete pointer arrays
	
	return clust_vect[index];		// return the index'th 7D array

}


