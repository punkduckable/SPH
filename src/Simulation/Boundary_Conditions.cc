#include "Simulation.h"
#include "Body/Body.h"
#include "Vector/Vector.h"
#include "Particle/Particle.h"
#include "Array.h"
#include <assert.h>

void Simulation::Apply_General_BCs(Body & Body_In, Array<Boundary_Condition> & BCs_In) {
  /* Cycle through the components of the BC array */
  unsigned Num_BCs = BCs_In.Get_Length();
  Vector Position;
  unsigned Num_Particles = Body_In.Get_Num_Particles();

  for(unsigned i = 0; i < Num_BCs; i++) {
    /* Cycle through the particles in Body_In. If a particle satisifies the
    condition then apply the effect. */

    double Mag_Normal_Vector = Magnitude(BCs_In[i].Condition_Plane_Normal_Vector);

    for(unsigned j = 0; j < Num_Particles; j++) {
      // Find distance from ith particle to the BC plane.
      if(BCs_In[i].Condition_Position_Type == Position_Type::Reference) {
        Position = Body_In[j].Get_X();
      } // if(BCs_In[i].Condition_Position_Type == Position_Type::Reference) {
      else {
        Position = Body_In[j].Get_x();
      } // else {
      double Distance_To_Plane = Dot_Product(Position, BCs_In[i].Condition_Plane_Normal_Vector)/Mag_Normal_Vector;


      // Now check if the distance satisifies the condition
      if( ((BCs_In[i].Condition_Inequality == Inequality::GE) &&
           (Distance_To_Plane >= BCs_In[i].Condition_Plane_Distance))

          ||

          ((BCs_In[i].Condition_Inequality == Inequality::LE) &&
           (Distance_To_Plane <= BCs_In[i].Condition_Plane_Distance)) ) {

        // If so, apply the effect to the ith particle.
        if(BCs_In[i].Effect_x == true) { Body_In[j].V[0] = BCs_In[i].Effect_Vector[0]; }
        if(BCs_In[i].Effect_y == true) { Body_In[j].V[1] = BCs_In[i].Effect_Vector[1];  }
        if(BCs_In[i].Effect_z == true) { Body_In[j].V[2] = BCs_In[i].Effect_Vector[2];  }
      } // if (Condition == True)
    } // for(unsigned j = 0; j < Num_Particles; j++) {
  } // for(unsigned i = 0; i < Num_BCs; i++) {
} // void Simulation::Apply_General_BCs(Body & Body_In, Array<Boundary_Condition> & BCs_In[i]) {



void Simulation::Apply_Box_BCs(Body & Box, Box_Properties & Box_Parameters) {
  /* This function sets the BC's for the six sides of the cube. The faces are
  x_plus, x_minus, y_plus, y_minus, z_plus, and z_minus. These faces are named
  based on the direction that the normal vector to that face points.


                            Y     X
                            |    /
                            |   /
                            |  /
                            | /
                        _ _ |/_ _ _ _ _ _ _ Z
                            /
                           /|
                          / |

  For example: the normal vector to the 'x_plus' and 'x_minus' faces point in
  the +x and -x directions, respectivly. */

  // This body must be a box.
  assert(Box.Get_Is_Box() == true);

  // Determine side lengths
  unsigned X_SIDE_LENGTH = Box.Get_X_SIDE_LENGTH();
  unsigned Y_SIDE_LENGTH = Box.Get_Y_SIDE_LENGTH();
  unsigned Z_SIDE_LENGTH = Box.Get_Z_SIDE_LENGTH();
  unsigned i,j,k;


  // +x face (i = X_SIDE_LENGTH-1)
  i = X_SIDE_LENGTH-1;
  for(j = 0; j < Y_SIDE_LENGTH; j++) {
    for(k = 0; k < Z_SIDE_LENGTH; k++) {
      Apply_Box_Particle_BCs(Box[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j], Box_Parameters.x_plus_BC);
    } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(j = 0; j < Y_SIDE_LENGTH; j++) {

  // -x face (i = 0)
  i = 0;
  for(j = 0; j < Y_SIDE_LENGTH; j++) {
    for(k = 0; k < Z_SIDE_LENGTH; k++) {
      Apply_Box_Particle_BCs(Box[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j], Box_Parameters.x_minus_BC);
    } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(j = 0; j < Y_SIDE_LENGTH; j++) {



  // +y face (j = y_Side_len-1)
  j = Y_SIDE_LENGTH-1;
  for(i = 0; i < X_SIDE_LENGTH; i++) {
    for(k = 0; k < Z_SIDE_LENGTH; k++) {
      Apply_Box_Particle_BCs(Box[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j], Box_Parameters.y_plus_BC);
    } //for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < X_SIDE_LENGTH; i++) {

  // -y face (j = 0)
  j = 0;
  for(i = 0; i < X_SIDE_LENGTH; i++) {
    for(k = 0; k < Z_SIDE_LENGTH; k++) {
      Apply_Box_Particle_BCs(Box[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j], Box_Parameters.y_minus_BC);
    } // for(k = 0; k < Z_SIDE_LENGTH; k++) {
  } // for(i = 0; i < X_SIDE_LENGTH; i++) {



  // +z face (k = Z_SIDE_LENGTH-1)
  k = Z_SIDE_LENGTH-1;
  for(i = 0; i < X_SIDE_LENGTH; i++) {
    for(j = 0; j < Y_SIDE_LENGTH; j++) {
      Apply_Box_Particle_BCs(Box[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j], Box_Parameters.z_plus_BC);
    } // for(j = 0; j < Y_SIDE_LENGTH; j++) {
  } // for(i = 0; i < X_SIDE_LENGTH; i++) {

  // -z face (k = 0)
  k = 0;
  for(i = 0; i < X_SIDE_LENGTH; i++) {
    for(j = 0; j < Y_SIDE_LENGTH; j++) {
      Apply_Box_Particle_BCs(Box[i*Y_SIDE_LENGTH*Z_SIDE_LENGTH + k*Y_SIDE_LENGTH + j], Box_Parameters.z_minus_BC);
    } // for(j = 0; j < Y_SIDE_LENGTH; j++) {
  } // for(i = 0; i < X_SIDE_LENGTH; i++) {
} // void Simulation::Apply_Box_BCs(Body & Box, Box_Properties & Box_Parameters) {



void Simulation::Apply_Box_Particle_BCs(Particle & P_In, Vector BC) {
  for(unsigned i = 0; i < 3; i++) {
    if(BC[i] == Free_BC_Box) { continue; }
    else { P_In.V[i] = BC[i]; }
  } // for(unsigned i = 0; i < 3; i++) {
} // void Simulation::Apply_Box_Particle_BCs(Particle & P_In, Vector BC) {
