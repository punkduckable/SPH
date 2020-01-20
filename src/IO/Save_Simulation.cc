#include "Save_Simulation.h"
#include "Body/Body.h"
#include "Particle/Particle.h"
#include "IO_Ops.h"
#include "Errors.h"

void IO::Save_Simulation(const Body * Bodies, const unsigned Num_Bodies) {
  /* This Function is used to save a simulation. This function prints all the
  information that is needed to re-create the state of a simulation. This is
  done by first printing information on the number of bodies, as well
  as their names. Once this is done, each Body gets a new file that
  contains all the information needed to re-create that Body object.
  The intent of this is to allow the user 'save' a simulation, allowing the user
  to - at a future point - pick up where they left off.*/

  /* Open a new file. This file will store the number of Bodys as
  well as their names. */
  FILE * File = fopen("./IO/Saves/Simulation_Data.txt","w");
  if(File == nullptr) {
    throw Cant_Open_File("Can't Open File Exception: Thrown by IO::Save_Simulation\n"
                         "For some reason, /IO/Saves/Simulation_Data.txt could not be opened\n");
  } // if(File == nullptr) {


  // Print number of Bodies to this file
  fprintf(File,   "Number of Bodies:    %u\n\n", Num_Bodies);

  // Now print each Body's name and essential information
  for(unsigned i = 0; i < Num_Bodies; i++) {
    fprintf(File, "Body %3u name:      %s\n", i, Bodies[i].Get_Name().c_str());
    fprintf(File, "     Is a Box:                %u\n",    Bodies[i].Get_Is_Box());

    if(Bodies[i].Get_Is_Box() == true) {
      fprintf(File, "          X_SIDE_LENGTH:      %u\n",  Bodies[i].Get_X_SIDE_LENGTH());
      fprintf(File, "          Y_SIDE_LENGTH:      %u\n",  Bodies[i].Get_Y_SIDE_LENGTH());
      fprintf(File, "          Z_SIDE_LENGTH:      %u\n",  Bodies[i].Get_Z_SIDE_LENGTH());
    } // if(Bodies[i].Get_Is_Box() == true) {

    fprintf(File, "     Is fixed in place:       %u\n",    Bodies[i].Get_Is_Fixed());
    fprintf(File, "     Is Damageable:           %u\n",    Bodies[i].Get_Damagable());
    fprintf(File, "     Number of particles:     %u\n",    Bodies[i].Get_Num_Particles());
    fprintf(File, "     Time steps per update:   %u\n",    Bodies[i].Get_Time_Steps_Between_Updates());

    fprintf(File, "     # times printed net external forces:         %u\n",    Bodies[i].Times_Printed_Net_External_Force);
    fprintf(File, "     # times printed particle forces:             %u\n",    Bodies[i].Times_Printed_Particle_Forces);
    fprintf(File, "     # times printed particle positions:          %u\n\n",  Bodies[i].Times_Printed_Particle_Positions);
  } // for(unsigned i = 0; i < Num_Bodies; i++) {

  // Now print each Body to its own file
  for(unsigned i = 0; i < Num_Bodies; i++) { IO::Save_Body(Bodies[i]); }

  fclose(File);
} // void IO::Save_Simulation(const Body * Bodies, const unsigned Num_Bodies) {



void IO::Save_Body(const Body & Body_In) {
  /* This function prints all data needed to reconstruct a particular
  Body. This information is printed to a file named after the
  Body that it's printing from. This prints the basic material
  properties of the body, as well as all the information that is
  needed to recreate all the particles in that Body. The intent of
  this is to allow the user 'save' a Body object. */

  // Get name + number of particles of passed Body object.
  const std::string Name = Body_In.Get_Name();
  const unsigned Num_Particles = Body_In.Get_Num_Particles();

  // Now make a file path for this Particle's file
  std::string File_Path = "./IO/Saves/";
  File_Path += Name;
  File_Path += ".txt";

  /* Now let's make a file for this Body. Note, if there is already
  a save file for this Body it will be overwritten by this new file */
  FILE * File = fopen(File_Path.c_str(), "w");
  if(File == nullptr) {
    char Buf[500];
    sprintf(Buf,
            "Can't Open File Exception: Thrown by IO::Save_Body\n"
            "For some reason, /IO/Saves/%s.txt couldn't be opened.\n",
            Name.c_str());
    throw Cant_Open_File(Buf);
  } // if(File == nullptr) {


  // Let's begin by printing the Body paramaters
  fprintf(File,   "Name:                         %s\n\n",  Name.c_str());

  fprintf(File,   "Is a Box:                     %u\n",    Body_In.Get_Is_Box());
  if(Body_In.Get_Is_Box() == true) {
    fprintf(File, "     X_SIDE_LENGTH:           %u\n",    Body_In.Get_X_SIDE_LENGTH());
    fprintf(File, "     Y_SIDE_LENGTH:           %u\n",    Body_In.Get_Y_SIDE_LENGTH());
    fprintf(File, "     Z_SIDE_LENGTH:           %u\n",    Body_In.Get_Z_SIDE_LENGTH());
  } //   if(Body_In.Get_Is_Box() == true) {
  fprintf(File,   "Is Fixed in place:            %u\n",    Body_In.Get_Is_Fixed());
  fprintf(File,   "Is Damageable:                %u\n\n",  Body_In.Get_Damagable());

  fprintf(File,   "       -- Kernel Parameters --\n");
  fprintf(File,   "Inter Particle Spacing:       %5lf\n",  Body_In.Get_Inter_Particle_Spacing());
  fprintf(File,   "Support Radius (IPS):         %u\n",    Body_In.Get_Support_Radius());
  fprintf(File,   "Support Radius (mm) aka h:    %5lf\n",  Body_In.Get_h());
  fprintf(File,   "Shape Function Amplitude:     %5lf\n\n",Body_In.Get_Shape_Function_Amplitude());

  fprintf(File,   "       -- Material Parameters --\n");
  fprintf(File,   "Material:                     %s\n",    Body_In.Get_Material().Name.c_str());
  fprintf(File,   "Lame parameter:               %5lf\n",  Body_In.Get_Lame());
  fprintf(File,   "Shear modulus (mu0):          %5lf\n",  Body_In.Get_mu0());
  fprintf(File,   "Viscosity (mu):               %5lf\n",  Body_In.Get_mu());
  fprintf(File,   "F_Index:                      %u\n",    Body_In.Get_F_Index());
  fprintf(File,   "Hourglass Stiffness (E):      %5lf\n",  Body_In.Get_E());
  fprintf(File,   "Material density:             %5lf\n",  Body_In.Get_density());
  fprintf(File,   "alpha (HG parameter):         %5lf\n",  Body_In.Get_alpha());
  fprintf(File,   "Tau (damage parameter):       %5lf\n\n",Body_In.Get_Tau());

  // Now let's print the number of particles
  fprintf(File,   "       -- Particles --\n");
  fprintf(File,   "Number of particles:          %u\n\n",    Body_In.Get_Num_Particles());

  // Finally, let's print the Box paramaters (should be removed if not using
  // a Box)
  //fprintf(File,   "X Side Length:                %u\n",    Simulation::X_SIDE_LENGTH);
  //fprintf(File,   "Y Side Length:                %u\n",    Simulation::Y_SIDE_LENGTH);
  //fprintf(File,   "Z Side Length:                %u\n\n",  Simulation::Z_SIDE_LENGTH);

  // Now let's print all particle data to the file
  for(unsigned i = 0; i < Num_Particles; i++) {
    IO::Save_Particle(Body_In[i], File);
  } // for(unsigned i = 0; i < Num_Particles; i++) {

  // We've now written the 'Particle_Data' file, we can close it.
  fclose(File);
} // void IO::Save_Body(const Body & Body_In) {



void IO::Save_Particle(const Particle & P_In, FILE * File) {
  /* This function prints all the information that is needed to re-create the
  input particle. Notably, this means that we do NOT need to print the first
  Piola Kirchoff stress tensor (P), the deformation gradient (F), any of the
  forces, or most of the neighbor arrays (the ID's are needed, the rest is not).
  P, F, and the forces are not needed because these are all calculated from
  stratch each iteration. Likewise, the neighbor array parameters can be
  recalculated if we know the neighbor IDs.

  Not storing this information in the File makes the file take up less/easier
  to read. This function assumes that the File has already been setup (with
  static particle class paramaters). */

  unsigned i;                                // index variable

  const Vector X = P_In.Get_X();
  const Vector x = P_In.Get_x();
  const Vector V = P_In.Get_V();
  const Tensor F_0 = P_In.Get_F(0);
  const Tensor F_1 = P_In.Get_F(1);

  // Print particle ID, dimensions
  fprintf(File,   "ID:                           %u\n",    P_In.Get_ID());
  fprintf(File,   "Mass:                         %5e\n",   P_In.Get_Mass());
  fprintf(File,   "Volume:                       %5e\n",   P_In.Get_Vol());
  fprintf(File,   "Radius:                       %5lf\n",  P_In.Get_Radius());

  // Print Particle dynamic properties
  fprintf(File,   "X:                            <%6.3lf %6.3lf %6.3lf>\n", X(0), X(1), X(2));
  fprintf(File,   "x:                            <%6.3lf %6.3lf %6.3lf>\n", x(0), x(1), x(2));
  fprintf(File,   "V:                            <%6.3lf %6.3lf %6.3lf>\n", V(0), V(1), V(2));
  fprintf(File,   "F[0]:                         |%6.3lf %6.3lf %6.3lf|\n", F_0(0,0), F_0(0,1), F_0(0,2));
  fprintf(File,   "                              |%6.3lf %6.3lf %6.3lf|\n", F_0(1,0), F_0(1,1), F_0(1,2));
  fprintf(File,   "                              |%6.3lf %6.3lf %6.3lf|\n", F_0(2,0), F_0(2,1), F_0(2,2));
  fprintf(File,   "F[1]:                         |%6.3lf %6.3lf %6.3lf|\n", F_1(0,0), F_1(0,1), F_1(0,2));
  fprintf(File,   "                              |%6.3lf %6.3lf %6.3lf|\n", F_1(1,0), F_1(1,1), F_1(1,2));
  fprintf(File,   "                              |%6.3lf %6.3lf %6.3lf|\n", F_1(2,0), F_1(2,1), F_1(2,2));


  // Damage paramaters
  fprintf(File,   "Stretch_H:                    %5lf\n",  P_In.Get_Stretch_H());
  fprintf(File,   "Stretch_M:                    %5lf\n",  P_In.Get_Stretch_M());
  fprintf(File,   "Stretch_Critical:             %5lf\n",  P_In.Get_Stretch_Critical());
  fprintf(File,   "D:                            %5lf\n",  P_In.Get_D());


  // BC information
  fprintf(File,   "Has Boundary Conditions:      <");
  for(unsigned i = 0; i < 3; i++) {
    if(P_In.Get_Has_BC(i) == true) { fprintf(File, "true   "); }
    else { fprintf(File, "false  "); }
  } // for(unsigned i = 0; i < 3; i++) {
  fprintf(File, ">\n");

  fprintf(File,   "Boundary Conditions:          <");
  for(unsigned i = 0; i < 3; i++) {
    if(P_In.Get_Has_BC(i) == true) { fprintf(File, "%6.3lf ", P_In.Get_BC(i)); }
    else { fprintf(File, "FREE   "); }
  } // for(unsigned i = 0; i < 3; i++) {
  fprintf(File, ">\n");


  // Now, let's figure out how many neighbors this particle has.
  unsigned Num_Neighbors = P_In.Get_Num_Neighbors();

  // Neighbor paramters
  fprintf(File,   "Number of neighbors:          %u\n", P_In.Get_Num_Neighbors());

  // Print neighbor IDs
  fprintf(File,   "Neighbor IDs                  ");
  for(i = 0; i < Num_Neighbors; i++) {
    fprintf(File,"%d ",P_In.Get_Neighbor_IDs(i));
  } // for(i = 0; i < Num_Neighbors; i++) {

  fprintf(File,"\n\n");
} // void IO::Save_Particle(const Particle & P_In, FILE * File) {
