#include "body.h"
#include "Vector3D.h"

// calculates gravitational force on body temp2 by temp1
auto Fg ( Body temp1, Body temp2 ) -> Vector3D {
  // in quadtree, temp1 represents this body
  // and temp2 represents the next body

  // get gravitation constant in units of ly^3 Msun^-1 (100 yr)^-2
  static constexpr double G{ 1.56700e-9 };
  
  // get the masses
  auto M1{ temp1.MASS };
  auto M2{ temp2.MASS };

  // this takes the position vector difference
  // from temp1 to temp2
  auto dr{ temp2.POSITION - temp1.POSITION };
  // and this is the magnitude
  auto dr_mag{ pow(dr.GetMagnitude(),2.) + 0.0001 };
  // 0.0001 for softening
		
  // here is the gravitational force on temp2
  auto F{ ( - G * M1 * M2 / ( pow( dr_mag, 1.5 ) ) ) * dr };
			
  return(F);
		
}
