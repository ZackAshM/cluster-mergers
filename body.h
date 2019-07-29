#pragma once

#include "Vector3D.h"

class Body
{
 public:

  // Body attributes
  double MASS;
  Vector3D POSITION;
  Vector3D VELOCITY;

  // class constructor
  Body( double mass = 0.0, Vector3D position = Vector3D(), Vector3D velocity = Vector3D() ) {
    MASS = mass;
    POSITION = position;
    VELOCITY = velocity;
  }
  ~Body(){};
  
};
