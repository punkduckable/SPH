
#if !defined(PARTICLE_UPDATE)
#define PARTICLE_UPDATE

#include "Particle_Helpers.h"

////////////////////////////////////////////////////////////////////////////////
// Friend functions (Update P, Update particle position)

void Particle_Helpers::Update_P(Particle_Array & Particles, const double dt) {
  /* The purpose of this function is to calculate the First Piola-Kirchhoff
  stress tensor for each particle in the Particles Particle_Array.

  This function assumes that the position each partilce in Particle's has been
  updated to the previous time step. It also assumes that each particle has
  a neighbor list, has calculate A^(-1), and has populated its dynamic array
  members. This function should not be called until these assumptions are valid.

  what are the arguments? This function accepts a Particle (P_In), a list of all
  particles in the current body, and the desired time step. This function uses
  these arguments to calculate P (the first Piola-Kirchhoff stress tensor). dt
  is used in calculating the viscosity. */

  // First, let's declare some local variables.
  Tensor F;                                      // Deformation gradient                 : unitless Tensor
  Tensor C;                                      // Richt-Cauchy stress tensor           : unitless Tensor
  Tensor S;                                      // Second Poila-Kirchhoff stress tensor : Mpa Tensor
  Tensor I = {1,0,0,                             // Identity tensor
              0,1,0,
              0,0,1};

  const double Tau = Particles.Get_Tau();                                      //        : unitless
  List<unsigned int> Damaged_Particle_List;      // Keeps track of which particles are newly damaged


  const double Lame = Particles.Get_Lame();      // Lame paramater                       : Mpa
  const double mu0 = Particles.Get_mu0();        // Shear modulus                        : Mpa
  const double mu = Particles.Get_mu();          // Viscosity                            : Mpa*s
  const unsigned int Num_Particles = Particles.Get_Num_Particles();

  Tensor F_Prime;                                // F time derivative                    : 1/s Tensor
  Tensor L;                                      // symmetric part of velocity gradient  : 1/s Tensor
  Tensor Visc;                                   // Viscosity correction term for P      : Mpa*s Tensor
  Vector rj;                                     // Displacemtn vector of jth neighbor   : mm Vector

  // Let's update each particle's stress tensors, keeping track of which
  // particles are damaged (in the Damaged_Particle_List)
  #pragma omp for
  for(unsigned int i = 0; i < Num_Particles; i++) {
    // First, Check if the current particle is damaged (if so, we skip this particle)
    if(Particles[i].D >= 1)
      continue;

    ////////////////////////////////////////////////////////////////////////////
    /* Now, we can calculate F for the current particle by cycling through
    that particle's neighbors. First, however, we need to reset F (from the
    previous iteration) */
    F = {0,0,0,
         0,0,0,
         0,0,0};

    const unsigned int Num_Neighbors = Particles[i].Num_Neighbors;
    for(unsigned int j = 0; j < Num_Neighbors; j++) {
      // Get neighbor ID and volume of jth particle
      unsigned Neighbor_ID = Particles[i].Neighbor_IDs[j]; // Index of jth neighbor.
      double V_j = Particles[Neighbor_ID].Vol;             // Volume of jth neighbor     : mm^3

      // Now find spatial distnace to jth particle
      rj = Particles[Neighbor_ID].x - Particles[i].x;                          //        : mm Vector

      // Calculate deformation gradient from jth particle.
      F += Dyadic_Product(rj, V_j*Particles[i].Grad_W[j]);                     //        : unitless Tensor
    } // for(unsigned int j = 0; j < Num_Neighbors; j++) {

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
    double Max_EigenValue = Max_Component(Eigenvalues(C, 'F'));                //        : unitless
    double Stretch_Max_Principle = sqrt(Max_EigenValue);                       //        : unitless

    // If this stretch is greater than max stretch, update particle's Max stretch.
    Particles[i].Stretch_M = Stretch_Max_Principle;
    if(Stretch_Max_Principle > Particles[i].Stretch_H)
      Particles[i].Stretch_H = Stretch_Max_Principle;

    // if damage is enabled and Max is greater than crticial then start adding damage
    if(Particles.Get_Damagable() == true) {
      if(Particles[i].Stretch_H > Particles[i].Stretch_Critical)
        Particles[i].D = exp(((Particles[i].Stretch_H - Particles[i].Stretch_Critical)*(Particles[i].Stretch_H - Particles[i].Stretch_Critical))/(Tau*Tau)) - 1;

      // If particle is fully damaged, add this particle to the damaged list
      // and move on (we'll remove it later)
      if(Particles[i].D >= 1) {
        Damaged_Particle_List.Add_Back(Particles[i].ID);
        continue;
      } // if(Particles[i].D >= 1) {
    } //   if(Particles.Get_Damagable() == true) {

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
      printf("Particle %d in %s has a non-positive Jacobian, J =  %lf.\n",Particles[i].ID, Particles.Get_Name().c_str(), J);

      // Let's also print this Particle's neighbor ID's (this helps debug issues)
      printf("Neighbor ID's: ");
      for(unsigned int j = 0; j < Particles[i].Num_Neighbors; j++)
        printf("%u ",Particles[i].Neighbor_IDs[j]);
      printf("\n");

      // Now we add the bad particles to the Damaged_Particle list (we'll remove
      // it once we have cycled through all partilces)
      Damaged_Particle_List.Add_Back(Particles[i].ID);
      continue;
    } //   if(J < 0) {

    S = (1- Particles[i].D)*(mu0*I + (-mu0 + 2.*Lame*log(J))*(C^(-1)));        //        : Mpa Tensor

    ////////////////////////////////////////////////////////////////////////////
    /* Calculate viscosity tensor:
    To do this, we need to calculate the time derivative of the deformation
    gradient. To get O(h^2) accuracy, we use a three point approximation for
    calculating this derivative. To use this three point derivative, we need
    to know the last two deformation gradients. These are stored in each particle
    in the 'F' array. The 'F_Index in the current particle array keeps track
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
    i_dt = Particles.Get_F_Index();
    if(i_dt == 0)
      i_2dt = 1;
    else
      i_2dt = 0;

    F_Prime = (1./(2.*dt))*(3*F - 4*Particles[i].F[i_dt] + Particles[i].F[i_2dt]);      // 1/s Tensor
    L = F_Prime*(F^(-1));                                                      //        : 1/s Tensor
    Visc = (J*mu)*(L + (L^(T))*(F^(-T)));                                      //        : Mpa Tensor

    ////////////////////////////////////////////////////////////////////////////
    /* Calculate P (First Piola-Kirchhoff stress tensor). Here we also update
    the current particle's P and F members. Note, however, that we need to be
    careful about which F we update (since each particle stores the last two
    F's). We really want to overwrite the older deformation gradient. Therefore,
    we overwrite the current particle's F tensor for two time steps ago with the
    new F. This way, the current particle stores F(t) (in Particles[i].F[i_2dt])
    and F(t-dt) (in Particles[i].F[i_dt]).*/
    Particles[i].P = (F*S + Visc)*Particles[i].A_Inv;                          //         : Mpa Tensor
    Particles[i].F[i_2dt] = F;                                                 //         : unitless Tensor
    //Particles[i].Visc = Visc*Particles[i].A_Inv;                               // For debugging
  } // for(unsigned int i = 0; i < Num_Particles; i++) {

  // Now we need to remove the damaged particles. To do this, we can one by one
  // have each thread remove its damaged particles
  #pragma omp critical
  while(Damaged_Particle_List.Node_Count() != 0)
    Particle_Helpers::Remove_Damaged_Particle(Particles[Damaged_Particle_List.Remove_Back()], Particles);

  // Explicit barrier to ensure that all broken particles have been removed
  // before any thread is allowed to move on.
  #pragma omp barrier
} // void Particle_Helpers::Update_P(Particle_Array & Particles, const double dt) {



void Particle_Helpers::Update_x(Particle_Array & Particles, const double dt) {
  /* This function is used to update the position (spatial position) of every
  partile in some particle array.

  To do this, we cycle through each particle in the particle array. For each
  particle, we calculate the internal, hourglassing, and external forces that
  are applied to the particle. Once the forces have been found, the particle's
  acceleration is calculated. From the acceleration, we use leapfrog integration
  to get a new velocity and position of the particle.

  This function assumes that every particle in the Particle's array has
  an updated P tensor. This function should not be run until this assumption
  is valid. */

  // First, define some local variables.
  const Vector g = {0,0,0};                  // Gravity                              : mm/s^2 Vector

  // Current (ith) particle properties
  Vector Force_Int;                              // Internal Force vector                : N Vector
  Vector Force_HG;                               // Hour-glass force                     : N Vector
  //Vector Force_Visc;                             // For debugging
  Tensor F_i;                                    // Deformation gradient                 : unitless Tensor
  Tensor P_i;                                    // First Piola-Kirchhoff stress tensor  : Mpa Tensor
  Vector a;                                      // acceleration                         : mm/s^2 Vector

  // Neighboring (jth) particle properties
  Tensor P_j;                                    // First Piola-Kirchhoff stress tensor  : Mpa Tensor
  Tensor F_j;                                    // Deformation gradient                 : unitless Tensor
  Vector rj;                                     // Displacement vector                  : mm Vector

  /* Particle array parameters:
  Note: these are all declared as constants to avoid accidential modification */
  const double alpha = Particles.Get_alpha();    // alpha member of the particles array  : unitless
  const double E = Particles.Get_E();            // Hourglass stiffness                  : Mpa
  const unsigned char F_Index = Particles.Get_F_Index();   // Keeps track of which F was most recently updated
  const unsigned int Num_Particles = Particles.Get_Num_Particles();

  // Damage variables
  List<unsigned int> Damaged_Particle_List;      // Keeps track of which particles have broken

  // Now loop through each partilce in the Partilces' array
  #pragma omp for
  for(unsigned int i = 0; i < Num_Particles; i++) {
    // First, check if particle is damaged (if so, we skip this particle)
    if( Particles[i].Get_D() >= 1)
      continue;

    // Now reset the force vectors
    Force_Int = {0,0,0};
    Force_HG = {0,0,0};
    //Force_Visc = {0,0,0};

    // Set up current particle properties
    double V_i = Particles[i].Get_Vol();         // volume of current particle           : mm^3
    F_i = Particles[i].Get_F(F_Index);
    P_i = Particles[i].Get_P();
    Vector * R = Particles[i].R;                 // Reference displacement array         : mm Vector
    double * Mag_R = Particles[i].Mag_R;         // Mag of reference displacment array   : mm
    double * W = Particles[i].W;                 // Kernel function array                : 1/mm^3
    Vector * Grad_W = Particles[i].Grad_W;       // Grad_W array                         : 1/mm^4 Vector

    const unsigned int Num_Neighbors = Particles[i].Num_Neighbors;

    for(unsigned int j = 0; j < Num_Neighbors; j++) {
      // Update Neighbor
      unsigned int Neighbor_ID = Particles[i].Neighbor_IDs[j];       // ID of current neighbor particle

      //////////////////////////////////////////////////////////////////////////
      /* Calculate Internal force */

      /* Note, each term in the internal force sum is multiplied by Vi. If  we
      we were to multiply through by Vi in this loop, we'd peerform this operation
      Num_Neighbors times. By moving it out of the summation (mutiplying Force_Int
      by Vi after the loop) we reduce the number of multiplications to 1, thereby
      reducing the number of FLOPs required to calculate the internal force and
      speeding up the program. */


      double V_j = Particles[Neighbor_ID].Vol;   // Volume of jth particle               : mm^3
      P_j = Particles[Neighbor_ID].P;                                          //        : Mpa Tensor
      Force_Int += (V_j)*((P_i + P_j)*Grad_W[j]);                              //        : N Vector
      //Force_Visc += (V_j)*((Particles[i].Visc + Particles[Neighbor_ID].Visc)*Grad_W[j]);// For debugging

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
      double delta_ij = Vector_Dot_Product(F_i*R[j], rj)/(Mag_rj) - Mag_rj;    //        : mm

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

      F_j = Particles[Neighbor_ID].F[F_Index];                                 //        : unitless Tensor
      double delta_ji = Vector_Dot_Product(F_j*R[j], rj)/(Mag_rj) - Mag_rj;    //        : mm

      /* Finally, we calculate the hour glass force. However, it should be
      noted that each term of Force_HG is multiplied by -(1/2), E, alpha,
      and Vi. However, these four quantities are constants. We can therefore
      pull these multiplications out of the summations (thereby saving
      several thousand floating point operations per particle!)*/
      Force_HG += (((V_j*W[j])/(Mag_R[j]*Mag_R[j]*Mag_rj))*                    //        : (1/mm) Vector
                  (delta_ij + delta_ji))*(rj);
    } // for(unsigned int j = 0; j < Num_Neighbors; j++) {
    Force_HG *= -.5*E*V_i*alpha;       // Each term in F_Hg is multiplied by this. Pulling out of sum improved runtime : N Vector
    Force_Int *= V_i;                  // Each term in F_Int is multiplied by Vi, pulling out of sum improved runtime  : N Vector
    //Force_Visc *= V_i;                 // for debugging

    /* Compute acceleration of particle at new position a(t_i+1).
    Note that all the forces we have calculated have been in units of Newtons.
    Our mass is in units of grams and we want the acceleration in units of
    mm/s^2. To get that, we note that 1N = 10^6(g*mm/s^2). Therefore, if we
    multiply our force, in Newtons, by 10^6 and then divide by the mass, in grams,
    then we get acceleration in mm/s^2. */
    a = ((1e+6)*(1./Particles[i].Mass))*(Force_Int                             //        : mm/s^2 Vector
                                        + Particles[i].Force_Contact
                                        + Particles[i].Force_Friction
                                        + Force_HG)
                                        + g;                                   // gravity

    /* Now update the velocity, position vectors. This is done using the
    'leap-frog' integration scheme. However, during the first step of this
    scheme, we need to use forward euler to get the initial velocity.*/
    if(Particles.Get_First_Time_Step() == true) {
      Particles.Set_First_Time_Step(false);
      Particles[i].V += (dt/2.)*a;                                             // velocity starts at t_i+1/2           : mm/s Vector
    } // if(Particles.Get_First_Time_Step() == true) {

    /* Before updating the velocity/position, let's check if the particle has
    diverged. This happens whenever any of the components of the acceleration
    vector are 'nan'. If they are, then we damage and remove this Particle.
    It should be noted that this is sort of a last resort mechanism. */
    if(std::isnan(a[0]) || std::isnan(a[1]) || std::isnan(a[2])) {
      printf("Particle %d in %s has a nan acceleration :(\n",Particles[i].ID, Particles.Get_Name().c_str());
      Particles[i].D = 1;
      Damaged_Particle_List.Add_Back(Particles[i].ID);
      continue;
    } //  if(std::isnan(a[0]) || std::isnan(a[1]) || std::isnan(a[2])) {

    Particles[i].x += dt*Particles[i].V;         // x_i+1 = x_i + dt*v_(i+1/2)           : mm Vector
    Particles[i].V += dt*a;                      // V_i+3/2 = V_i+1/2 + dt*a(t_i+1)      : mm/s Vector
    Particles[i].a = a;                          // update acceleration vector           : mm/s^2 Vector

    if(Simulation::Print_Forces == true) {
      Particles[i].Force_Int = Force_Int;          // update Internal force                : N Vector
      Particles[i].Force_HG = Force_HG;            // update Hourglassing force            : N Vector
      //Particles[i].Force_Visc = Force_Visc;        // update Viscosity force               : N Vector
    } // if(Simulation::Print_Forces == true) {
  } // for(int i = 0; i < Num_Particles; i++) {

    // Now we need to remove the damaged particles. To do this, we can one by one
    // have each thread remove its damaged particles
    #pragma omp critical
    while(Damaged_Particle_List.Node_Count() != 0)
      Particle_Helpers::Remove_Damaged_Particle(Particles[Damaged_Particle_List.Remove_Back()], Particles);

    /* Note, there is no explicit barrier here because the next kernel, which
    updates each particle's timestep counter, does not use any data being
    written in the critial second above and because that kernel contains
    its own implied barrier. Therefore, once the data above is called upon, we
    can be sure that all particles have been damaged. Therefore, we don't
    need a barrier. */
} // void Particle_Helpers::Update_x(const Particle_Array & Particles, const double dt) {

#endif
