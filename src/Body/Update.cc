#include "Body.h"
#include "Simulation/Simulation.h"
#include "Particle/Particle.h"
#include "Vector/Vector.h"
#include "List.h"
#include "Diagnostics/Operation_Count.h"
#include "Errors.h"
#include <math.h>
#include <assert.h>

// Static prototypes.
static void Calculate_Force(Vector & F,
                            const double V_j,
                            const Tensor & T1,
                            const Tensor & T2,
                            const Vector & Grad_Wj);

static double Calculate_Delta(const Tensor & F,
                              const Vector & R_j,
                              const Vector & r_j,
                              const double Mag_rj);



////////////////////////////////////////////////////////////////////////////////
// Update methods

void Body::Update_P(const double dt) {
  /* The purpose of this function is to calculate the First Piola-Kirchhoff
  stress tensor for each particle in a Body.

  This function assumes that the position each particle in Particle's has been
  updated to the previous time step. It also assumes that each particle has
  a neighbor list, has calculate A^(-1), and has populated its dynamic array
  members. This function should not be called until these assumptions are valid.

  what are the arguments? This function accepts a the desired time step. dt
  is used to calculate the viscosity. We must have dt > 0 */
  assert(dt > 0);

  // First, let's declare some local variables.
  Tensor F;                                      // Deformation gradient                 : unitless Tensor
  Tensor C;                                      // Richt-Cauchy stress tensor           : unitless Tensor
  Tensor S;                                      // Second Poila-Kirchhoff stress tensor : Mpa Tensor
  Tensor I = {1,0,0,                             // Identity tensor
              0,1,0,
              0,0,1};


  List<unsigned> Damaged_Particle_List;          // Keeps track of which particles are newly damaged

  Tensor F_Prime;                                // F time derivative                    : 1/s Tensor
  Tensor L;                                      // symmetric part of velocity gradient  : 1/s Tensor
  Tensor Visc;                                   // Viscosity correction term for P      : Mpa*s Tensor
  Vector rj;                                     // Displacement vector of jth neighbor   : mm Vector


  // Material parameters
  const double Lame = Body_Material.Lame;        // Lame parameter                        : Mpa
  const double mu0 = Body_Material.mu0;          // Shear modulus                        : Mpa

  // Let's update each particle's stress tensors, keeping track of which
  // particles are damaged (in the Damaged_Particle_List)
  #pragma omp for
  for(unsigned i = 0; i < Num_Particles; i++) {
    // First, Check if the current particle is damaged (if so, we skip this particle)
    if(Particles[i].D >= 1) { continue; }


    ////////////////////////////////////////////////////////////////////////////
    /* Now, we can calculate F for the current particle by cycling through
    that particle's neighbors. First, however, we need to reset F (from the
    previous iteration) */
    F = {0,0,0,
         0,0,0,
         0,0,0};

    const unsigned Num_Neighbors = Particles[i].Num_Neighbors;
    for(unsigned j = 0; j < Num_Neighbors; j++) {
      // Get neighbor ID and volume of jth particle
      unsigned Neighbor_ID = Particles[i].Neighbor_IDs[j]; // Index of jth neighbor.
      double V_j = Particles[Neighbor_ID].Volume;          // Volume of jth neighbor     : mm^3

      // Now find spatial distnace to jth particle
      rj = Particles[Neighbor_ID].x - Particles[i].x;                          //        : mm Vector

      // Calculate deformation gradient from jth particle.
      F += Dyadic_Product(rj, V_j*Particles[i].Grad_W[j]);                     //        : unitless Tensor
    } // for(unsigned j = 0; j < Num_Neighbors; j++) {

    // Deformation gradient with correction
    F *= Particles[i].A_Inv;                                                   //        : unitless Tensor


    ////////////////////////////////////////////////////////////////////////////
    /* Calculate Damage:
    To do this, we first need to find the principle stretch. To do this, we
    need to find the square root of the biggest eigenvalue of the (right)
    Cauchy Green strain tensor C. Note, C will be used for later calculations */

    C = (F^(T))*F;                               // Right Cauchy-Green strain tensor     : unitless Tensor
    double J = Determinant(F);                   // J is det of F                        : unitless

    // Calculate current principle stretch
    double Max_EigenValue = Max_Component(Eigenvalues(C));                     //        : unitless
    double Stretch_Max_Principle = sqrt(Max_EigenValue);                       //        : unitless

    // If this stretch is greater than max stretch, update particle's Max stretch.
    Particles[i].Stretch_M = Stretch_Max_Principle;
    if(Stretch_Max_Principle > Particles[i].Stretch_H) {
      Particles[i].Stretch_H = Stretch_Max_Principle;
    } // if(Stretch_Max_Principle > Particles[i].Stretch_H) {

    // if damage is enabled and Max is greater than crticial then start adding damage
    if((*this).Is_Damageable == true) {
      if(Particles[i].Stretch_H > Particles[i].Stretch_Critical) {
        /* Note: Set_Tau requires that Tau != 0. Therefore, dividing by Tau
        should be well defined. */
        Particles[i].D = exp(((Particles[i].Stretch_H - Particles[i].Stretch_Critical)*(Particles[i].Stretch_H - Particles[i].Stretch_Critical))/(Tau*Tau)) - 1.;

        #ifdef OPERATION_COUNT
          // 3 subtractions, 2 multiplications, 1 division, and one exp in the calculation above.
          OP_Count::Subtraction += 3;
          OP_Count::Multiplication += 2;
          OP_Count::Division += 1;
          OP_Count::Exp += 1;
        #endif
      } // if(Particles[i].Stretch_H > Particles[i].Stretch_Critical) {

      // If particle is fully damaged, add this particle to the damaged list
      // and move on (we'll remove it later)
      if(Particles[i].D >= 1) {
        Particles[i].D = 1;  // Set D to 1 (paraview gets angry if D is too big)
        Damaged_Particle_List.Push_Back(Particles[i].ID);
        continue;
      } // if(Particles[i].D >= 1) {
    } // if((*this).Is_Damageable == true) {


    ////////////////////////////////////////////////////////////////////////////
    /* Now that we have calculated the deformation gradient, we need to calculate
    the first Piola-Kirchhoff stess tensor. To do this, however, we need to
    find the Second Piola-Kirchhoff stress tensor and the Viscosity tensor. */

    /* Calculate Second Piola-Kirchoff stress tensor:
    It should be noted that calculating this tensor requires taking the log of J.
    In theory, J will always be positive.... However, it is theoretically possible
    for this to not be the case. Thus, before calculating S, we check if J is
    non-positive. If it is, then we treat this particle as damaged. */
    if(J <= 0) {
      Particles[i].D = 1;
      printf("Particle %d in %s has a non-positive Jacobian, J =  %lf.\n",Particles[i].ID, Name.c_str(), J);

      // Let's also print this Particle's neighbor ID's (this helps debug issues)
      printf("Neighbor ID's: ");
      for(unsigned j = 0; j < Particles[i].Num_Neighbors; j++) { printf("%u ",Particles[i].Neighbor_IDs[j]); }
      printf("\n");

      // Now we add the bad particles to the Damaged_Particle list (we'll remove
      // it once we have cycled through all particles)
      Damaged_Particle_List.Push_Back(Particles[i].ID);
      continue;
    } // if(J <= 0) {

    /* Note: C should be invertible if J != 0.
    C = (F^T)F. Thus, det(C) = det(F^T)det(F) = det(F)^2 = J^2. Thus, if J != 0
    then det(C) != 0, and C^(-1) is well defined. */
    S = (1 - Particles[i].D)*(mu0*I + (-mu0 + 2.*Lame*log(J))*(C^(-1)));       //        : Mpa Tensor

    #ifdef OPERATION_COUNT
      /* 1 addition, 1 subtraction, 1 log, 2 multiplications in the calculation above.
      There are many other operations in there, but they're all done with operator
      overloading and are, therefore, counted elsewhere. */
      OP_Count::Subtraction += 1;
      OP_Count::Addition += 1;
      OP_Count::Log += 1;
      OP_Count::Multiplication += 2;
    #endif



    ////////////////////////////////////////////////////////////////////////////
    /* Calculate viscosity tensor:
    To do this, we need to calculate the time derivative of the deformation
    gradient. To get O(h^2) accuracy, we use a three point approximation for
    calculating this derivative. To use this three point derivative, we need
    to know the last two deformation gradients. These are stored in each particle
    in the 'F' array. The 'F_Index' in the body keeps track
    of which member of this array was last updated (therefore telling us which
    element of the array is F from the last time step, F(t-dt), and which is from two
    time steps ago, F(t-2*dt)). The current deformation gradient, F(t), is stored
    in the local F variable (We just calculated it). We then use the following
    three point apporixmation to calculate (dF/dt),
          dF/dt  = (1/(2*dt))*(3*F(t) - 4*F(t-dt) + F(t - 2dt)) + O(h^2)
    Using this, we can then calculate the spatial velocity gradient l with
          l = F'*F^-1.
    Which can then be used to calculate d (the rate of strain tensor),
          d = (1/2)*(l + L^T)
    Which can then be used to calculate the Viscosity tensor,
           Viscosity = 2*J*Mu*d*F^(-T)
    */
    unsigned char i_dt, i_2dt;
    i_dt = F_Index;
    if(i_dt == 0) { i_2dt = 1; }
    else {          i_2dt = 0; }

    /* Note: Update_P requires that dt > 0 (see assertion near the top of this
    fucntion). Thus, dividing by dt should be ewell defined.

    Note: F^(-T) should be well defined as long as J != 0 (which we check for
    above). J = det(F), so if J != 0 then F should be invertible. */
    F_Prime = (1./(2.*dt))*(3*F - 4*Particles[i].F[i_dt] + Particles[i].F[i_2dt]);      // 1/s Tensor
    L = F_Prime*(F^(-1));                                                      //        : 1/s Tensor
    Visc = (J*mu)*(L + (L^(T))*(F^(-T)));                                      //        : Mpa Tensor

    #ifdef OPERATION_COUNT
      /* F_Prime : 1 division, 1 multiplication (all other operations use operator overloading)
      Visc       : 1 multiplication (other operations use operator overloading) */
      OP_Count::Multiplication += 2;
      OP_Count::Division += 1;
    #endif

    ////////////////////////////////////////////////////////////////////////////
    /* Calculate P (First Piola-Kirchhoff stress tensor). Here we also update
    the current particle's P and F members. Note, however, that we need to be
    careful about which F we update (since each particle stores the last two
    F's). We really want to overwrite the older deformation gradient. Therefore,
    we overwrite the current particle's F tensor for two time steps ago with the
    new F. This way, the current particle stores F(t) (in Particles[i].F[i_2dt])
    and F(t-dt) (in Particles[i].F[i_dt]). */
    Particles[i].P = (F*S + Visc)*Particles[i].A_Inv;                          //         : Mpa Tensor
    Particles[i].F[i_2dt] = F;                                                 //         : unitless Tensor
    Particles[i].Visc = Visc*Particles[i].A_Inv;
  } // for(unsigned i = 0; i < Num_Particles; i++) {

  // Now we need to remove the damaged particles. To do this, we can one by one
  // have each thread remove its damaged particles
  #pragma omp critical
  Remove_Damaged_Particles(Damaged_Particle_List);

  // Explicit barrier to ensure that all broken particles have been removed
  // before any thread is allowed to move on.
  #pragma omp barrier
} // void Body::Update_P(const double dt) {



void Body::Update_x(const double dt) {
  /* This function is used to update the position (spatial position) of every
  partile in a Body.

  To do this, we cycle through the particles in the Body. For each
  particle, we calculate the internal, hourglassing, and external forces that
  are applied to the particle. Once the forces have been found, the particle's
  acceleration is calculated. From the acceleration, we use leapfrog integration
  to get a new velocity and position of the particle.

  This function assumes that every particle in the Body has an updated P tensor.
  This function should not be run until this assumption is valid. */

  // Current (ith) particle properties
  Vector Force_Internal;                         // Internal Force vector                : N Vector
  Vector Force_Hourglass;                        // Hour-glass force                     : N Vector
  Vector Force_Viscosity;                        // Viscosity force.                     : N Vector
  Vector a;                                      // acceleration                         : mm/s^2 Vector

  // Neighboring (jth) particle properties
  Tensor P_j;                                    // First Piola-Kirchhoff stress tensor  : Mpa Tensor
  Tensor F_j;                                    // Deformation gradient                 : unitless Tensor
  Vector rj;                                     // Displacement vector                  : mm Vector

  // Material parameters
  const double E = (*this).Body_Material.E;      // Hourglass stiffness/Youngs Modulus   : Mpa

  // Damage variables
  List<unsigned> Damaged_Particle_List;      // Keeps track of which particles have broken

  // Now loop through each particle in the particles' array
  #pragma omp for
  for(unsigned i = 0; i < Num_Particles; i++) {
    // First, check if particle is damaged (if so, we skip this particle)
    if( Particles[i].Get_D() >= 1) { continue; }

    // Now reset the force vectors
    Force_Internal = {0,0,0};
    Force_Hourglass = {0,0,0};
    Force_Viscosity = {0,0,0};

    // Set up current particle properties
    double V_i = Particles[i].Get_Volume();      // volume of current particle           : mm^3
    const Tensor & F_i = Particles[i].Get_F(F_Index);
    const Tensor & P_i = Particles[i].Get_P();
    const Tensor & Visc = Particles[i].Visc;

    Vector * R = Particles[i].R;                 // Reference displacement array         : mm Vector
    double * Mag_R = Particles[i].Mag_R;         // Mag of reference displacment array   : mm
    double * W = Particles[i].W;                 // Kernel function array                : 1/mm^3
    Vector * Grad_W = Particles[i].Grad_W;       // Grad_W array                         : 1/mm^4 Vector

    const unsigned Num_Neighbors = Particles[i].Num_Neighbors;

    for(unsigned j = 0; j < Num_Neighbors; j++) {
      // Update Neighbor
      unsigned Neighbor_ID = Particles[i].Neighbor_IDs[j];       // ID of current neighbor particle

      //////////////////////////////////////////////////////////////////////////
      /* Calculate Internal force */

      /* Note, each term in the internal force sum is multiplied by Vi. If  we
      we were to multiply through by Vi in this loop, we'd peerform this operation
      Num_Neighbors times. By moving it out of the summation (mutiplying Force_Internal
      by Vi after the loop) we reduce the number of multiplications to 1, thereby
      reducing the number of FLOPs required to calculate the internal force and
      speeding up the program. */


      double V_j = Particles[Neighbor_ID].Volume;// Volume of jth particle               : mm^3
      P_j = Particles[Neighbor_ID].Get_P();                                    //        : Mpa Tensor

      // Force_Internal += V_j*((P_i + P_j)*Grad_W[j]);
      Calculate_Force(Force_Internal , V_j, P_i , P_j                        , Grad_W[j]);         // N Vector

      // Force_Viscosity += V_j*((Visc + Particles[Neighbor_ID].Visc)*Grad_W[j]);
      Calculate_Force(Force_Viscosity, V_j, Visc, Particles[Neighbor_ID].Visc, Grad_W[j]);         // N Vector

      //////////////////////////////////////////////////////////////////////////
      /* Calculate Hour Glass force */

      /* Here we calculate delta_ij.
      Before discussing this, let us establish the following definitions
            r_ij = rj - ri
            R_ij = Rj - Ri
            rk = spacial position of kth particle
            Rk = ref position of kth particle
            Fk = deformation gradient of kth particle)
      delta_ij is given by,
            delta_ij = (Error_ij dot r_ij)/|r_ij|
      where
            Error_ij = Fi*R_ij - r_ij.
      since the dot product is distributive, notice that
            delta_ij = ((Fi*R_ij - r_ij) dot (r_ij))/|r_ij|
                    = (Fi*R_ij dot r_ij)/|r_ij| - (r_ij dot r_ij)/|r_ij|
                    = (Fi*R_ij dot r_ij)/|r_ij| - |r_ij|^2/|r_ij|
                    = (Fi*R_ij dot r_ij)/|r_ij| - |r_ij|
      Calcualating delta this way actually uses fewer floating point operations
      and should therefore perform better. */
      rj = Particles[Neighbor_ID].x - Particles[i].x;                          //        : mm Vector
      double Mag_rj = Magnitude(rj);                                           //        : mm
      if(Mag_rj == 0) {
        char Buf[512];
        sprintf(Buf,
                "Divide By Zero Exception: thrown by %s while calling Body::Update_x\n"
                "In %s, Particle %u and %u have the same x coordinate. As such, rj = 0,\n"
                "which means that we can calculate delta_ij (which is needed to update x).\n",
                (*this).Name.c_str(), (*this).Name.c_str(), Neighbor_ID, i);
        throw Divide_By_Zero(Buf);
      } // if(Mag_rj == 0) {

      // double delta_ij = Dot_Product(F_i*R[j], rj)/(Mag_rj) - Mag_rj;
      double delta_ij = Calculate_Delta(F_i, R[j], rj, Mag_rj);                //        : mm

      /* Here we calculate delta_ji.
            delta_ji = ( Error_ji dot r_ji )/|r_ji|
      With
            Error_ji = F_j*R_ji - r_ji
      Notice that we need the jth particles deformation gradient to calculate
      Error_ji. We also need the jth particles R_ji and r_ij. We could calculate all
      of these, but we can save some time by making a few clever observations.
      From the definion of R,
            R_ji = X_i - X_j = -(X_j - X_i) = -R_ij
      Likewise, from the defintion of r, we can decduce
            r_ij = -r_ji
      Thus,
          Error_ji = F_j*R_ji - r_ji
                   = F_j*(-R_ij) + r_ij
                   = -F_j*(R_ij) + r_ij
      Then, using the fact that the dot product is distributive, we get
            delta_ji = ( Error_ji dot r_ji ) / |r_ji|
                     = ( Error_ji dot -r_ij ) / |r_ij|
                     = -( Error_ji dot r_ij ) / |r_ij|
                     = -( (-F_j*R_ij + r_ij ) dot r_ij )/|r_ij|
                     = -( -F_j*R_ij dot r_ij) / |r_ij| - (r_ij dot r_ij)/|r_ij|
                     = (F_j*R_ij dot r_ij) / |r_ij| - |r_ij|^2/|r_ij|
                     = (F_j*R_ij dot r_ij) / |r_ij| - |r_ij|
      Computing delta_ji this way uses fewer arithmetic operations and should
      therefore improve performance.

      Note: this calculation requires that Mag_rj != 0. We check for this
      condition above, however. */
      F_j = Particles[Neighbor_ID].F[F_Index];                                 //        : unitless Tensor

      // double delta_ji = Dot_Product(F_j*R[j], rj)/(Mag_rj) - Mag_rj;
      double delta_ji = Calculate_Delta(F_j, R[j], rj, Mag_rj);                //        : mm

      /* Finally, we calculate the hour glass force. However, it should be
      noted that each term of Force_Hourglass is multiplied by -(1/2), E, alpha,
      and Vi. However, these four quantities are constants. We can therefore
      pull these multiplications out of the summations (thereby saving
      several thousand floating point operations per particle!)

      Note: this assumes that Mag_R and Mag_rj are non-zero. We already checked
      for the latter. The former should be true so long as the bodies were setup
      properly. */
      Force_Hourglass += (((V_j*W[j])/(Mag_R[j]*Mag_R[j]*Mag_rj))*             //        : (1/mm) Vector
                         (delta_ij + delta_ji))*rj;


      #ifdef OPERATION_COUNT
        /* 3 multiplications, 1 division, 1 addition in the calculation above
        (the other multiplication uses operator overloading) */
        OP_Count::Multiplication += 3;
        OP_Count::Division += 1;
        OP_Count::Addition += 1;
      #endif
    } // for(unsigned j = 0; j < Num_Neighbors; j++) {
    Force_Hourglass *= -.5*E*V_i*alpha;  // Each term in F_Hg is multiplied by this. Pulling out of sum improved runtime : N Vector
    Force_Internal *= V_i;               // Each term in F_Int is multiplied by Vi, pulling out of sum improved runtime  : N Vector
    Force_Viscosity *= V_i;              // Viscious force

    /* Compute acceleration of particle at new position a(t_i+1).
    Note that all the forces we have calculated have been in units of Newtons.
    Our mass is in units of grams and we want the acceleration in units of
    mm/s^2. To get that, we note that 1N = 10^6(g*mm/s^2). Therefore, if we
    multiply our force, in Newtons, by 10^6 and then divide by the mass, in grams,
    then we get acceleration in mm/s^2.

    Note: This assumes that Particles[i].Mass != 0. However, the Set_Mass
    function of the Particle class requires that Mass != 0, and Get_Mass will
    only work if a mass has been set. Thus, this assumption should hold. */
    a = ((1e+6)*(1./Particles[i].Get_Mass()))*(Force_Internal                  //        : mm/s^2 Vector
                                             + Particles[i].Force_Contact
                                             + Particles[i].Force_Friction
                                             + Force_Hourglass);
    /* If gravity is enabled, add that in. */
    if((*this).Gravity_Enabled == true) { a += Body::g; }

    #ifdef OPERATION_COUNT
      /* Force_Hourglass : 3 multiplications
      a                  : 1 division, 1 multiplication (eveything else is operator oveloading)*/
      OP_Count::Multiplication += 4;
      OP_Count::Division += 1;
    #endif

    /* Now update the velocity, position vectors. This is done using the
    'leap-frog' integration scheme. However, during the first step of this
    scheme, we need to use forward euler to get the initial velocity.*/
    if(First_Time_Step == true) {
      First_Time_Step = false;
      Particles[i].V += (dt/2.)*a;                                             // velocity starts at t_i+1/2           : mm/s Vector

      #ifdef OPERATION_COUNT
        /* 1 multiplication to calculate V above (multiplication by a uses
        operator oveloading) */
        OP_Count::Division += 1;
      #endif
    } // if(First_Time_Step == true) {

    /* Before updating the velocity/position, let's check if the particle has
    diverged. This happens whenever any of the components of the acceleration
    vector are 'nan'. If they are, then we damage and remove this Particle.
    It should be noted that this is sort of a last resort mechanism. */
    if(std::isnan(a[0]) || std::isnan(a[1]) || std::isnan(a[2])) {
      Particles[i].Print();
      printf("Particle %d in %s has a nan acceleration :(\n",Particles[i].ID, Name.c_str());
      Particles[i].D = 1;
      Damaged_Particle_List.Push_Back(Particles[i].ID);
      continue;
    } //  if(std::isnan(a[0]) || std::isnan(a[1]) || std::isnan(a[2])) {

    Particles[i].x += dt*Particles[i].V;         // x_i+1 = x_i + dt*v_(i+1/2)           : mm Vector
    Particles[i].V += (dt)*a;                    // V_i+3/2 = V_i+1/2 + dt*a(t_i+1)      : mm/s Vector
    Particles[i].a = a;                          // update acceleration vector           : mm/s^2 Vector

    if(Simulation::Print_Particle_Forces == true || Simulation::Print_Body_Forces == true || Simulation::Print_Body_Torques == true) {
      Particles[i].Force_Internal  = Force_Internal;    // update Internal force                : N Vector
      Particles[i].Force_Hourglass = Force_Hourglass;   // update Hourglassing force            : N Vector
      Particles[i].Force_Viscosity = Force_Viscosity;   // update Viscosity force               : N Vector
    } // if(Simulation::Print_Particle_Forces == true || Simulation::Print_Body_Forces == true || Simulation::Print_Body_Torques == true) {
  } // for(int i = 0; i < Num_Particles; i++) {

  // Now we need to remove the damaged particles. To do this, we can one by one
  // have each thread remove its damaged particles
  #pragma omp critical
  Remove_Damaged_Particles(Damaged_Particle_List);

  /* Note, there is no explicit barrier here because the next kernel, which
  updates each particle's timestep counter, does not use any data being
  written in the critial second above and because that kernel contains
  its own implied barrier. Therefore, once the data above is called upon, we
  can be sure that all particles have been damaged. Therefore, we don't
  need a barrier. */
} // void Body::Update_x(const double dt) {



////////////////////////////////////////////////////////////////////////////////
// Helper functions.
/* These are functions that compute quantities that would otherwise be computed
using operator overloading and would otherwise incur temporary variables. */

static void Calculate_Force(Vector & F,
                            const double V_j,
                            const Tensor & T1,
                            const Tensor & T2,
                            const Vector & Grad_Wj) {
  /* This function is used to calculate Force_Internal and Force_Viscosity due
  to one of a particle's neighbors. These quantities are calculated by the
  following expression:
      Force_Internal  += V_j*((P_i + P_j)*Grad_W[j])
      Force_Viscosity += V_j*((Visc + Particles[Neighbor_ID].Visc)*Grad_W[j])
  Thus, this function computes the following:
      F += V_j*((T1 + T2)*GradW_j) without any
  The goal is to eliminate any use of temporary objects and (hopefully) improve
  runtime.

  update_x is the only thing that should call this function. */

  /* Note: Tensors are stored in ROW MAJOR ordering. Thus, we want to change
  rows as infrequently. */
  const double * T1_Ar = T1.Get_Ar();
  const double * T2_Ar = T2.Get_Ar();
  const double * Grad_Wj_Ar = Grad_Wj.Get_Ar();


  F[0] += + V_j*( (T1_Ar[0*3 + 0] + T2_Ar[0*3 + 0])*Grad_Wj_Ar[0] +
                  (T1_Ar[0*3 + 1] + T2_Ar[0*3 + 1])*Grad_Wj_Ar[1] +
                  (T1_Ar[0*3 + 2] + T2_Ar[0*3 + 2])*Grad_Wj_Ar[2] );

  F[1] += + V_j*( (T1_Ar[1*3 + 0] + T2_Ar[1*3 + 0])*Grad_Wj_Ar[0] +
                  (T1_Ar[1*3 + 1] + T2_Ar[1*3 + 1])*Grad_Wj_Ar[1] +
                  (T1_Ar[1*3 + 2] + T2_Ar[1*3 + 2])*Grad_Wj_Ar[2] );

  F[2] += + V_j*( (T1_Ar[2*3 + 0] + T2_Ar[2*3 + 0])*Grad_Wj_Ar[0] +
                  (T1_Ar[2*3 + 1] + T2_Ar[2*3 + 1])*Grad_Wj_Ar[1] +
                  (T1_Ar[2*3 + 2] + T2_Ar[2*3 + 2])*Grad_Wj_Ar[2] );

  #ifdef OPERATION_COUNT
    /* Each component above uses 6 additions and 4 multiplications */
    OP_Count::Multiplication += 12;
    OP_Count::Addition += 18;
  #endif
} // static void Calculate_Force(Tensor & F,...



static double Calculate_Delta(const Tensor & F,
                              const Vector & R_j,
                              const Vector & rj,
                              const double Mag_rj) {
  /* This function computes delta_ij and delta_ji. By definition,
        delta_ji = Dot_Product(F_j*R[j], rj)/(Mag_rj) - Mag_rj
        delta_ij = Dot_Product(F_i*R[j], rj)/(Mag_rj) - Mag_rj
  Thus, this function returns the following quantity:
        Dot_Product(F*R_j, rj)/(Mag_rj) - Mag_rj
  The goal is to eliminate any use of temporary objects.

  update_x is the ONLY thing that should call this function. */

  const double * F_Ar  = F.Get_Ar();
  const double * Rj_Ar = R_j.Get_Ar();
  const double * rj_Ar = rj.Get_Ar();


  double FRj_0 = F_Ar[0*3 + 0]*Rj_Ar[0] + F_Ar[0*3 + 1]*Rj_Ar[1] + F_Ar[0*3 + 2]*Rj_Ar[2];
  double FRj_1 = F_Ar[1*3 + 0]*Rj_Ar[0] + F_Ar[1*3 + 1]*Rj_Ar[1] + F_Ar[1*3 + 2]*Rj_Ar[2];
  double FRj_2 = F_Ar[2*3 + 0]*Rj_Ar[0] + F_Ar[2*3 + 1]*Rj_Ar[1] + F_Ar[2*3 + 2]*Rj_Ar[2];

  #ifdef OPERATION_COUNT
    // This includes the calculations above and below (return statement)
    OP_Count::Multiplication += 12;
    OP_Count::Addition += 8;
    OP_Count::Division += 1;
    OP_Count::Subtraction += 1;
  #endif

  return (FRj_0*rj_Ar[0] + FRj_1*rj_Ar[1] + FRj_2*rj_Ar[2])/(Mag_rj) - Mag_rj;
} // static double Calculate_Delta(const Tensor & F,
