/* VecFRK4: a function to calculate fourth-order Runge-Kutta
   for 1-dimensional F = dp/dt = f(x(t),v(t)) (Newton force law) 
   
   Both x and v increments must be calculated at each step in 
   the integration of the equation, and xold and vold must be
   reset after each step 
 --P. Gorham: updated 3/20/2007 for Physics 305, add explicit time dependence,
correct error in previous version  */
// this version for c++ Vector3D class usage 4/9/14

using namespace std;       
#include <iostream>        
#define _USE_MATH_DEFINES  
#include <cmath>  
#include "Vector3D.h"

// Here is the primary function

Vector3D VecFRK4(int ytype, Vector3D (*f_x)(int, double, Vector3D, Vector3D),
			  Vector3D (*f_v)(int, double, Vector3D, Vector3D),
			  int i, double t, Vector3D xold, Vector3D vold, double dt)

/*	ytype = 0 for x-equation, 1 for v-equation
  	f_x = user supplied function for dx/dt=v
	f_v: force equation as noted above
	i: some generic index used by parent function
	t=time, xold, vold, prior position and velocity
   	dt = step size

   Vector3D (*f_xv)(double, Vector3D, Vector3D) ==> pointer to a function of
   two Vector3Ds that returns a Vector3D object. Either dv/dt or dx/dt.

	4th order general form:
    y(t+dt) = y(t) + 1/6 * (k1 + 2k2 + 2k3 + k4)
			k1 = dt * f(t,y)
            k2 = dt * f(t+dt/2, y(t) + k1/2)
        	k3 = dt * f(t+dt/2, y(t) + k2/2)
            k4 = dt * f(t+dt, y(t) + k3)
*/

{
	Vector3D k1x,k2x,k3x,k4x,k1v,k2v,k3v,k4v;

		k1x = dt * (*f_x)(i, t, xold, vold);
		k1v = dt * (*f_v)(i, t, xold, vold);
		
		k2x = dt * (*f_x)(i, t+dt/2.0, xold+0.5*k1x, vold+0.5*k1v);
		k2v = dt * (*f_v)(i, t+dt/2.0, xold+0.5*k1x, vold+0.5*k1v);
		
		k3x = dt * (*f_x)(i, t+dt/2.0, xold+0.5*k2x, vold+0.5*k2v);
		k3v = dt * (*f_v)(i, t+dt/2.0, xold+0.5*k2x, vold+0.5*k2v);
		
		k4x = dt * (*f_x)(i, t+dt, xold+k3x, vold+k3v);
		k4v = dt * (*f_v)(i, t+dt, xold+k3x, vold+k3v);


	if(ytype==0){
		return( (1.0/6.0) * (k1x+2.0*k2x+2.0*k3x+k4x) ); 
		} else { 
		return( (1.0/6.0) * (k1v+2.0*k2v+2.0*k3v+k4v) ); 
	}
}


