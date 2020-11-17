#include "Simulation.h"
#include "Body/Body.h"
#include "Particle/Particle.h"
#include "Vector/Vector.h"
#include "IO/Save_Simulation.h"
#include "IO/Load_Simulation.h"
#include "Errors.h"
#if defined(_OPENMP)
  #include <omp.h>
#endif

// Prototypes for functions that are local to this file
namespace Simulation {
  void Export_Bodies_Data(Body * Bodies,
                          const unsigned Num_Bodies,
                          const unsigned time_steps);
} // namespace Simulation {



void Simulation::Run(void) {
  //////////////////////////////////////////////////////////////////////////////
  // Simulation variables

  // Loop indicies
  unsigned p,b,t;

  // Simulation time measurement variables
  TIME_TYPE time1,
            time2,
            simulation_time = 0,
            update_BC_time = 0,
            update_P_time = 0,
            contact_time = 0,
            update_x_time = 0,
            Print_time = 0;

  Body * Bodies;                                 // Will point to the Bodies's for this simulation


  //////////////////////////////////////////////////////////////////////////////
  // Simulation start up.
  Setup(&Bodies);



  //////////////////////////////////////////////////////////////////////////////
  // Run time steps
  printf(         "\nRunning %d time steps....\n", Num_Time_Steps);
  time1 = Get_Time();

  // time step loop.
  #pragma omp parallel default(shared) private(b, p, t) firstprivate(Num_Bodies, Num_Time_Steps, dt, TimeSteps_Between_Prints)
  {
    for(t = 0; t < Num_Time_Steps; t++) {
      //////////////////////////////////////////////////////////////////////////
      // Export Bodies data
      #pragma omp single nowait
      { time2 = Get_Time(); }

      if(t%Simulation::TimeSteps_Between_Prints == 0) {
        #pragma omp single nowait
        { printf("%d time steps complete\n",t); }

        Simulation::Export_Bodies_Data(Bodies, Num_Bodies, t);
      } // if(t%TimeSteps_Between_Prints == 0) {

      #pragma omp single nowait
      { Print_time += Time_Since(time2); }



      //////////////////////////////////////////////////////////////////////////
      // Apply Boundary conditions
      #pragma omp single nowait
      { time2 = Get_Time(); }

      // Note: apply BCs uses a parallel for loop
      for(b = 0; b < Num_Bodies; b++) { Bodies[b].Apply_BCs();  }

      #pragma omp single nowait
      { update_BC_time += Time_Since(time2); }



      //////////////////////////////////////////////////////////////////////////
      // Update Stress tensor (P)

      #pragma omp single nowait
      { time2 = Get_Time(); }

      for(b = 0; b < Num_Bodies; b++) {
        // Note: We don't update P for Bodys that are fixed in place
        if(Bodies[b].Get_Is_Fixed() == true) { continue; }

        else {
          /* Update each Particles's P tensor.

          Note: the Update_P method has an orphaned for loop (and takes care of
          removing damaged particles in parallel, damaged particles are not
          removed until every particle's P tensor has been updated. This makes
          the code parallelizable and determinstic) */
          Bodies[b].Update_P(dt);
        } // else
      } // for(b = 0; b < Num_Bodies; b++) {

      #pragma omp single nowait
      { update_P_time += Time_Since(time2); }



      //////////////////////////////////////////////////////////////////////////
      // Contact
      /* Here we enable particle-particle contact. To do this, we cycle through
      each Body. For the mth body, we check if any of its particles are
      in contact with any of the particles in the ith body for i > m. We only
      use i > m so that we only run the contact algorithm on each part of
      Bodys once. Further, we only calculate the contact forces for the
      mth Body if that body is being updated this time step. */

      #pragma omp single nowait
      { time2 = Get_Time(); }

      /* First, we need to set each particle's contact force to zero. It should
      be noted that we only do this for a particular Body if that body
      is updating it's position this cycle. Otherwise, since the force won't be
      used for anything, there's no reason to waste CPU cycles setting that
      bodies's particle's contact forces to zero. */
      for(b = 0; b < Num_Bodies; b++) {
        unsigned Num_Particles = Bodies[b].Get_Num_Particles();

        #pragma omp for
        for(p = 0; p < Num_Particles; p++) {
          (Bodies[b])[p].Force_Contact = {0,0,0};
          (Bodies[b])[p].Force_Friction = {0,0,0};
        } // for(p = 0; p < (Bodies[b]).Get_Num_Particles(); p++) {
      } // for(b = 0; b < Num_Bodies; b++) {

      /* Now we can apply the contact algorithm. Note that this must be applied
      every time step no matter what (so that bodies that update each step can
      are proprly updated/have the right forces applied each timestpe) */
      for(unsigned b1 = 0; b1 < Num_Bodies - 1; b1++) {
        for(unsigned b2 = b1 + 1; b2 < Num_Bodies; b2++) {
          Body::Contact(Bodies[b2], Bodies[b1]);
        } // for(unsigned b2 = b1 + 1; b2 < Num_Bodies; b2++) {
      } // for(unsinged b1 = 0; b1 < Num_Bodies - 1; b1++) {

      #pragma omp single nowait
      { contact_time += Time_Since(time2); }



      //////////////////////////////////////////////////////////////////////////
      // Update Position (x)

      #pragma omp single nowait
      { time2 = Get_Time(); }

      for(b = 0; b < Num_Bodies; b++) {
        // Note: we don't update P for Bodies that are fixed in place
        if(Bodies[b].Get_Is_Fixed() == true) { continue; }

        else {
          /* First, update the 'F_Index' for the current Body. This
          controls which member of each particle's 'F' array is the 'newest'. */
          #pragma omp single
          { Bodies[b].Increment_F_Index(); }

          // Now update the position of each particle in this body.
          Bodies[b].Update_x(dt);
        } // else {
      } // for(b = 0; b < Num_Bodies; b++) {

      #pragma omp single nowait
      { update_x_time += Time_Since(time2); }
    } // for(t = 0; t < Num_Time_Steps; t++) {



    ////////////////////////////////////////////////////////////////////////////
    // Export the bodies data for the final configuration of the simulation

    #pragma omp single nowait
    { printf("%d time steps complete\n",t); }

    Simulation::Export_Bodies_Data(Bodies, Num_Bodies, t);
  } // #pragma omp parallel


  printf(         "Done!\n");
  simulation_time = Time_Since(time1);

  // If saving is enabled, Dump particle data to file
  if(Save_Simulation_To_File == 1) {
    IO::Save_Simulation(Bodies, Num_Bodies);
  } // if(Save_Simulation_To_File == 1) {

  // Print timing data
  #if defined(_OPENMP)
    printf(       "\nIt took %lf s to perform %u Particle time steps \n",simulation_time, Num_Time_Steps);
    printf(       "%lfs to update BC's\n", update_BC_time);
    printf(       "%lfs to update P\n", update_P_time);
    printf(       "%lfs for Contact\n", contact_time);
    printf(       "%lfs to update x\n", update_x_time);
    printf(       "%lfs to print data to files\n", Print_time);

  #else
    unsigned long MS_Iter,                                           // These are used to store the number of miliseconds that
                  MS_BC,                                             // it took to execute each of the major operations in
                  MS_P,                                              // the code. These are only used the the code is executed
                  MS_Contact,                                        // sequentially
                  MS_x,
                  MS_Print;

    MS_Iter = (unsigned long)((double)time1 / (double)CLOCKS_PER_MS);
    MS_BC = (unsigned long)((double)update_BC_time / (double)CLOCKS_PER_MS);
    MS_P = (unsigned long)((double)update_P_time / (double)CLOCKS_PER_MS);
    MS_Contact = (unsigned long)((double)contact_time / (double)CLOCKS_PER_MS);
    MS_x = (unsigned long)((double)update_x_time / (double)CLOCKS_PER_MS);
    MS_Print = (unsigned long)((double)Print_time / (double)CLOCKS_PER_MS);

    printf(         "\nIt took %lu ms to perform %u Particle time steps \n",MS_Iter, Num_Time_Steps);
    printf(         "%lums to update BC's\n", MS_BC);
    printf(         "%lums to update P\n", MS_P);
    printf(         "%lums for Contact\n", MS_Contact);
    printf(         "%lums to update x\n", MS_x);
    printf(         "%lums to print data to files\n", MS_Print);
  #endif

  delete [] Bodies;
} // void Simulation::Run(void) {





void Simulation::Export_Bodies_Data(Body * Bodies, unsigned Num_Bodies, const unsigned time_steps) {
  /* Function Description:
  This function, as the name implies, exports data for each body in a simulation.
  Position data is always printed. Wheather or not we print Force or
  Net External Force data depends on the simulation paramaters
  Print_Prticle_Force and Print_Next_External_Forces, respectivly.

  Simulation::Run is the only function  that should call this function */
  #pragma omp for nowait
  for(unsigned b = 0; b < Num_Bodies; b++) {
    try {
                                                      Bodies[b].Export_Particle_Positions();
      if(Simulation::Print_Particle_Forces == true) { Bodies[b].Export_Particle_Forces();}
      if(Simulation::Print_Body_Forces == true) {     Bodies[b].Export_Body_Forces(time_steps); }
      if(Simulation::Print_Body_Torques == true) {    Bodies[b].Export_Body_Torques(time_steps); }
    } // try {

    catch(Exception & Error_In) {
      printf("%s\n", Error_In.what());
      abort();
    } // catch(Exception & Error_In) {
  } // for(unsigned b = 0; b < Num_Bodies; b++ ) {
} // void Simulation::Export_Bodies_Data(Body * Bodies, unsigned Num_Bodies, const unsigned time_steps) {
