// Three dimensional Vector c++ class, based on original from
//  http://www.gamedev.net/topic/487050-c-vector-class/
// modified for use in Physics 305
// 4/3/2013  P. Gorham
// some updates 3/10/15 PG
// 4/2/15: changed way unit vectors are created
// 4/3/19 fixed bug in update method, needed a return value

#ifndef _VECTOR
#define _VECTOR

using namespace std;
#include <iostream>
#include <cmath>
#include <string>

class Vector3D
{
private:
	double x, y, z;  // these are only visible by members of class
public:
// public members can be used by external calls
//default constructor, note all components nitialized to 0.0 unless
// other values are present in parentheses
	Vector3D(double X = 0.0, double Y = 0.0, double Z = 0.0)
	{
		x = X;
		y = Y;
		z = Z;
	}
	~Vector3D(){};  // destructor, removes the object from memory
	
// get x component
	double Getx() 
	{
		return(x);
	}
	
// get y component
	double Gety() 
	{
		return(y);
	}

// get z component
	double Getz()
	{
		return(z);
	}

//calculate and return the magnitude of this vector
	double GetMagnitude()
	{
		return sqrt(x*x + y*y + z*z);
	}

//multiply this vector by a scalar
// the "const" at the end of the function declaration
// is necessary to ensure that the function does not
// change the class object it operates with
	Vector3D operator*(double num) const    // overloads the *
	{
		return Vector3D(x*num, y*num, z*num); 
	}
	
//multiply a vector by a scalar, friend functions allow for less
// restrictive conditions, eg., numerical values * vectors etc.
// friend functions have access to private variables but are
// not directly part of the class. In this case, we need the
// friend operator to ensure that reverse ordering of the "*" works
// with the double class
    friend Vector3D operator*(double num, Vector3D const &vec)
    {
         return Vector3D(vec.x * num, vec.y * num, vec.z * num);
    }
	

//add two vectors
	Vector3D operator+(const Vector3D &vec) const   // overloads the +
	{
		return Vector3D(x+vec.x, y+vec.y, z+vec.z);
	}

//subtract two vectors
	Vector3D operator-(const Vector3D &vec) const   // overloads the -
	{
		return Vector3D(x-vec.x, y-vec.y, z-vec.z);
	}

// get a new unit vector based on the original vector
	Vector3D unitVector3D()
	{
		double magnitude = sqrt(x*x + y*y + z*z);
		double tmpx = x/magnitude;
		double tmpy = y/magnitude;
		double tmpz = z/magnitude;
		return Vector3D(tmpx,tmpy,tmpz);
	}
// change a vector's components
	Vector3D update(double X , double Y, double Z )
	{
		x = X;
		y = Y;
		z = Z;
        return *this; // return a reference to the current object
	}
	
//calculate and return dot product
	double dot(const Vector3D &vec) const
	{
		return x*vec.x + y*vec.y + z*vec.z;
	}
	
//calculate and return cross product
	Vector3D cross(const Vector3D &vec) const
	{
		return Vector3D(y*vec.z - z*vec.y,
				z*vec.x - x*vec.z,
				x*vec.y - y*vec.x);
	}

//print a 3vector
	void print()
	{
		cout << x <<", "<< y <<", "<< z <<endl;
	}	
	
//print using cout
	string printt()
	{
		return "(" + to_string(x) + ", " + to_string(y) + ", " + to_string(z) + ")";
	}
	
};

#endif

