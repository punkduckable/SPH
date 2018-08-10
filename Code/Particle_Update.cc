
#if !defined(PARTICLE_UPDATE)
#define PARTICLE_UPDATE

#include "Particle_Helpers.h"

////////////////////////////////////////////////////////////////////////////////
// Friend functions (Update P, Update particle position)

void Particle_Helpers::Update_P(Particle & P_In, Particle_Array & Particles, const double dt) {
  // Check if particle is damaged (if so, we skip this particle)
  if(P_In.D >= 1)
    return;

  /* The purpose of this function is to calculate the First Piola-Kirchhoff
  stress tensor for the particle P_In.

  This function assumes that the position of each of P_In's neighbors has
  been updated to the previous time step. It also assumes that each particle has
  a neighbor list, has calculate A^(-1), and has populated its dynamic array
  members. Finally, it assumes that the static member variables k1, k2,
  and mu0 have been set.This function should not be called until these
  assumptions are valid.

  what are the arguments? This function accepts a Particle (P_In), a list of all
  particles in the current body, and the desired time step. This function uses
  these arguments to calculate P (the first Piola-Kirchhoff stress tensor) */

  /* First, let's set up the local variables that will be used to update the
  particle's position */
  double V_j;                                    // Volume of jth neighbor               : mm^3
  unsigned int Neighbor_ID;                      // Index of jth neighbor.

  Tensor F = {0,0,0,                             // Deformation gradient                 : unitless Tensor
              0,0,0,
              0,0,0};

  Tensor C;                                      // Richt-Cauchy stress tensor           : unitless Tensor
  double J;                                      // Deformation gradient determinant     : unitless
  Tensor S;                                      // Second Poila-Kirchhoff stress tensor : Mpa Tensor
  Tensor I = {1,0,0,
              0,1,0,
              0,0,1};                            // Identity tensor

  double Stretch_Max_Principle;                                                //        : unitless
  const double Tau = Particles.Get_Tau();                                       //        : unitless

  const double Lame = Particles.Get_Lame();      // Lame paramater                       : Mpa
  const double mu0 = Particles.Get_mu0();        // Shear modulus                        : Mpa
  const double mu = Particles.Get_mu();          // Viscosity                            : Mpa*s

  Tensor F_Prime;                                // F time derivative                    : 1/s Tensor
  Tensor L;                                      // symmetric part of velocity gradient  : 1/s Tensor
  Tensor Visc;                                   // Viscosity correction term for P      : Mpa*s Tensor
  Vector *Grad_W = P_In.Grad_W;                  // Pointer to P_In's Grad_W array.      : 1/mm Vector
  Vector rj;                                     // Displacemtn vector of jth neighbor   : mm Vector

  //////////////////////////////////////////////////////////////////////////////
  /* Now, we can calculate F by cycling through the neighbors. The contribution
  to F by the jth neighbor is dj (Dyadic Product) V_j Grad_W(Rj, h) */
  for(unsigned int j = 0; j < P_In.Num_Neighbors; j++) {
    Neighbor_ID = P_In.Neighbor_IDs[j];
    V_j = Particles[Neighbor_ID].Vol;                                          //        : mm^3
    rj = Particles[Neighbor_ID].x - P_In.x;                                    //        : mm Vector
    F += Dyadic_Product(rj, V_j*Grad_W[j]);                                    //        : unitless Tensor
  } // for(unsigned int j = 0; j < Num_Neighbors; j++) {

  // Deformation gradient with correction
  F *= P_In.A_Inv;                                                             //        : unitless Tensor

  //////////////////////////////////////////////////////////////////////////////
  /* Calculate Damage:

  To do this, we first need to find the principle stretch. To do this, we
  need to find the square root of the biggest eigenvalue of the (right)
  Cauchy Green strain tensor C. Note, C will be used for later calculations */
  C = (F^(T))*F;                                 // Right Cauchy-Green strain tensor     : unitless Tensor
  J = Determinant(F);                            // J is det of F                        : unitless Tensor

  // Calculate current principle stretch
  double Max_EigenValue = Max_Component(Eigenvalues(C, 'F'));

  Stretch_Max_Principle = sqrt(Max_EigenValue);

  // If this stretch is greater than max stretch, update particle's Max stretch.
  P_In.Stretch_M = Stretch_Max_Principle;
  if(Stretch_Max_Principle > P_In.Stretch_H)
    P_In.Stretch_H = Stretch_Max_Principle;

    // if damage is enabled and Max is greater than crticial then start adding damage
  if(Particles.Get_Damagable() == true) {
    if(P_In.Stretch_H > P_In.Stretch_Critical)
      P_In.D = exp(((P_In.Stretch_H - P_In.Stretch_Critical)*(P_In.Stretch_H - P_In.Stretch_Critical))/(Tau*Tau)) - 1;

    // If particle is fully damaged, remove it from array.
    if(P_In.D >= 1) {
      Remove_Damaged_Particle(P_In, Particles);
      return;
    } // if(P_In.D >= 1) {
  } //   if(Particles.Get_Damagable() == true) {

  //////////////////////////////////////////////////////////////////////////////
  /* Now that we have calculated the deformation gradient, we need to calculate
  the first Piola-Kirchhoff stess tensor. To do this, however, we need to
  find the Second Piola-Kirchhoff stress tensor and the Viscosity term. */

  /* Calculate Second Piola-Kirchoff stress tensor:
  It should be noted that calculating this tensor requires taking the log of J.
  In theory, J will always be positive.... However, it is theoretically possible
  for this to not be the case. Thus, before calculating S, we check if J is
  non-positive. If it is, then we treat this particle as damaged and remove it */
  if(J <= 0) {
    P_In.D = 1;
    printf("Particle %d has a non-positive Jacobian, J =  %lf.\n",P_In.ID, J);

    // Let's also print this Particle's neighbor ID's (this helps debug issues)
    printf("Neighbor ID's: ");
    for(unsigned int i = 0; i < P_In.Num_Neighbors; i++)
      printf("%u ",P_In.Neighbor_IDs[i]);
    printf("\n");

    // Now let's remove this bad particle
    Remove_Damaged_Particle(P_In, Particles);
    return;
  } //   if(J < 0) {

  S = (1-P_In.D)*(mu0*I + (-mu0 + 2.*Lame*log(J))*(C^(-1)));                   //        : Mpa Tensor

  /* Calculate viscosity tensor:
  To do this, we need to calculate the time derivative of the deformation
  gradient. To get O(h^2) accuracy, we use a three point approximation for
  calculating this derivative. To use this three point derivative, we need
  to know the last two deformation gradients. These are stored in each particle
  in the 'F' array. The 'F_Counter' in the current particle array keeps track
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
  i_dt = Particles.Get_F_Counter();
  if(i_dt == 0)
    i_2dt = 1;
  else
    i_2dt = 0;

  F_Prime = (1./(2.*dt))*(3*F - 4*P_In.F[i_dt] + P_In.F[i_2dt]);               //        : 1/s Tensor
  L = F_Prime*(F^(-1));                                                        //        : 1/s Tensor
  Visc = (J*mu)*(L + (L^(T))*(F^(-T)));                                        //        : Mpa Tensor

  /* Calculate P (First Piola-Kirchhoff stress tensor), send it and F to P_In.
  Note, however, that we need to be careful about which F we update (since P_In
  has a two element array that stores the last two F's). We know that the most
  recent F is stored in P_In.F[i_dt]. We really want to overwrite the 'old'
  element of P_In's F array. therefore, we make P_In.F[i_2dt] equal to F (the
  local F that we just calculated). This way, P_In stores F(t) (in P_In.F[i_2dt])
  and F(t-dt) (in P_In.F[i_dt]).*/
  P_In.P = (F*S + Visc)*P_In.A_Inv;                                            //         : Mpa Tensor
  P_In.F[i_2dt] = F;                                                           //         : unitless Tensor
  P_In.Visc = Visc*P_In.A_Inv;                                                 // For debugging

} // void Particle_Helpers::Update_P(Particle & P_In, Particle_Array & Particles, const double dt) {



void Particle_Helpers::Update_x(Particle & P_In, Particle_Array & Particles, const double dt) {
  // Check if particle is damaged (if so, we skip this particle)
  if( P_In.D >= 1)
    return;

  /* This function assumes that every particle in the Particle's array has
  an updated P tensor. Likewise, it assumes that the E and alpha static
  member variables have been set. This function should not be run until
  these assumptions are valid. */

  P_In.Force_Int = {0,0,0};                      // Internal Force vector                : N Vector
  P_In.Force_HG = {0,0,0};                       // Hour-glass force                     : N Vector
  P_In.Force_Visc = {0,0,0};                     // For debugging

  const Vector g = {0,-9810,0};                  // Gravity                              : mm/s^2 Vector

  unsigned int Neighbor_ID;                      // ID of current neighbor particle (in paritlce's array)

  /* Jth particle variables */
  double V_j;                                    // Volume of jth particle               : mm^3
  Tensor P_j;                                    // First Piola-Kirchhoff stress tensor  : Mpa Tensor
  Tensor F_j;                                    // Deformation gradient                 : unitless Tensor
  Vector rj;                                     // Displacement vector                  : mm Vector

  /* P_In aliases (ith particle variables).
  notice that P_In/P_i does not chnage throughout this function. Therefore,
  all P_In variables are declared as consts to avoid accidential modification */
  const double alpha = Particles.Get_alpha();    // alpha member of the particles array  : unitless
  const double E = Particles.Get_E();            // Hourglass stiffness                  : Mpa

  const double V_i = P_In.Vol;                   // Volume of P_In                       : mm^3

  const Tensor P_i = P_In.P;                     // First Piola-Kirchhoff stress tensor  : Mpa Tensor
  const unsigned char Current_F = Particles.Get_F_Counter();         // Keeps track of which F was most recently updated
  const Tensor F_i = P_In.F[Current_F];          // Deformation gradient                 : unitless Tensor

  const Vector * R = P_In.R;                     // Reference displacement array         : mm Vector
  const Vector * Grad_W = P_In.Grad_W;           // Grad_W array                         : 1/mm Vector

  /* Hour glass variables */
  double Mag_rj;                                                               //        : mm
  double delta_ij;                                                             //        : mm
  double delta_ji;                                                             //        : mm

  for(unsigned int j = 0; j < P_In.Num_Neighbors; j++) {
    // Update Neighbor
    Neighbor_ID = P_In.Neighbor_IDs[j];

    ////////////////////////////////////////////////////////////////////////////
    /* Calculate Internal force */

    /* Note, each term in the internal force sum is multiplied by Vi. If  we
    we were to multiply through by Vi in this loop, we'd peerform this operation
    Num_Neighbors times. By moving it out of the summation (mutiplying Force_Int
    by Vi after the loop) we reduce the number of multiplications to 1, thereby
    reducing the number of FLOPs required to calculate the internal force and
    speeding up the program. */

    V_j = Particles[Neighbor_ID].Vol;                                          //        : mm^3
    P_j = Particles[Neighbor_ID].P;                                            //        : Mpa Tensor
    P_In.Force_Int += (V_j)*((P_i + P_j)*Grad_W[j]);                           //        : N Vector
    P_In.Force_Visc += (V_j)*((P_In.Visc + Particles[Neighbor_ID].Visc)*Grad_W[j]); // For debugging

    ////////////////////////////////////////////////////////////////////////////
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
    rj = Particles[Neighbor_ID].x - P_In.x;                                    //        : mm Vector
    Mag_rj = Magnitude(rj);                                                    //        : mm
    delta_ij = Vector_Dot_Product(F_i*R[j], rj)/(Mag_rj) - Mag_rj;             //        : mm

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
    therefore improve performnace.
    */

    F_j = Particles[Neighbor_ID].F[Current_F];                                 //        : unitless Tensor
    delta_ji = Vector_Dot_Product(F_j*R[j], rj)/(Mag_rj) - Mag_rj;//: mm

    /* Finally, we calculate the hour glass force. However, it should be
    noted that each term of Force_HG is multiplied by -(1/2), E, alpha,
    and Vi. However, these four quantities are constants. We can therefore
    pull these multiplications out of the summations (thereby saving
    several thousand floating point operations per particle!)*/
    P_In.Force_HG += (((V_j*P_In.W[j])/(P_In.Mag_R[j]*P_In.Mag_R[j]*Mag_rj))*  //        : (1/mm) Vector
                (delta_ij + delta_ji))*(rj);
  } // for(unsigned int j = 0; j < Num_Neighbors; j++) {
  P_In.Force_HG *= -.5*E*V_i*alpha;    // Each term in F_Hg is multiplied by this. Pulling out of sum improved runtime : N Vector
  P_In.Force_Int *= V_i;               // Each term in F_Int is multiplied by Vi, pulling out of sum improved runtime  : N Vector
  P_In.Force_Visc *= V_i;              // for debugging

  /* Compute acceleration of particle at new position a(t_i+1).
  Note that all the forces we have calculated have been in units of Newtons.
  Our mass is in units of grams and we want the acceleration in units of
  mm/s^2. To get that, we note that 1N = 10^6(g*mm/s^2). Therefore, if we
  multiply our force, in Newtons, by 10^6 and then divide by the mass, in grams,
  then we get acceleration in mm/s^2. */
  P_In.a = ((1e+6)*(1./P_In.Mass))*(P_In.Force_Int                             //        : mm/s^2 Vector
                                  + P_In.Force_Contact
                                  + P_In.Force_HG)
                                  + g;                               // gravity

  /* Now update the velocity, position vectors. This is done using the
  'leap-frog' integration scheme. However, during the first step of this
  scheme, we need to use forward euler to get the initial velocity.*/
  if(P_In.First_Time_Step == true) {
    P_In.First_Time_Step = false;
    P_In.V += (dt/2.)*P_In.a;              // velocity starts at t_i+1/2           : mm/s Vector
  } //   if(P_In.First_Time_Step == true) {

  /* Before updating the velocity/position, let's check if the particle has
  diverged. This happens whenever any of the components of the acceleration
  vector are 'nan'. If they are, then we damage and remove this Particle.
  It should be noted that this is sort of a last resort mechanism. */
  if(std::isnan(P_In.a[0]) || std::isnan(P_In.a[1]) || std::isnan(P_In.a[2])) {
    printf("Particle %d's acceleration is nan :(\n",P_In.ID);
    P_In.D = 1;
    Remove_Damaged_Particle(P_In, Particles);
    return;
  } //  if(std::isnan(acceleration[0]) || std::isnan(acceleration[1]) || std::isnan(acceleration[2])) {

  P_In.x += dt*P_In.V;                           // x_i+1 = x_i + dt*v_(i+1/2)           : mm Vector
  P_In.V += dt*P_In.a;                           // V_i+3/2 = V_i+1/2 + dt*a(t_i+1)      : mm/s Vector
} // void Particle_Helpers::Update_x(Particle & P_In, const Particle_Array & Particles, const double dt) {

#endif
