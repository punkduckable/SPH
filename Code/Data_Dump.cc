#if !defined(DATA_DUMP_SOURCE)
#define DATA_DUMP_SOURCE

#include "Data_Dump.h"

void Data_Dump::Print_Particle_Array_To_File(const Particle_Array & Particles) {
  /* This function prints all data needed to reconstruct the Particle's array
  to a file. The intent of this is to allow the user to essentially 'save' the
  state of the program, allowing the user to - at a future point - pick up where
  they last saved/left off. */

  // Get number of particles from passed Particle array.
  const unsigned int Num_Particles = Particles.Get_Num_Particles();

  /* First, let's create the file. Note, if another 'Particle_Data' file already
  exists then this overwrites/deletes that file. */
  FILE * File = fopen("../Files/Particle_Data.txt", "w");

  // Let's begin by printing the 'static' particle class paramaters
  fprintf(File,   "Inter Particle Spacing:       %5lf\n",  Particles.Get_Inter_Particle_Spacing());
  fprintf(File,   "Support Radius (mm):          %5lf\n",  Particles.Get_h());
  fprintf(File,   "Support Radius (IPS):         %u\n",    Particles.Get_Support_Radius());
  fprintf(File,   "Shape Function Amplitude:     %5lf\n",  Particles.Get_Shape_Function_Amplitude());
  fprintf(File,   "Lame parameter:               %5lf\n",  Particles.Get_Lame());
  fprintf(File,   "Shear modulus (mu0):          %5lf\n",  Particles.Get_mu0());
  fprintf(File,   "Viscosity (mu):               %5lf\n",  Particles.Get_mu());
  fprintf(File,   "Hourglass Stiffness (E):      %5lf\n",  Particles.Get_E());
  fprintf(File,   "alpha (HG parameter):         %5lf\n",  Particles.Get_alpha());
  fprintf(File,   "Tau (damage parameter):       %5lf\n",  Particles.Get_Tau());

  // Now let's print the number of particles
  fprintf(File,   "\n");
  fprintf(File,   "Number of particles:          %u\n",    Particles.Get_Num_Particles());

  // Finally, let's print the cuboid paramaters (should be removed if not using
  // a cuboid)
  fprintf(File,   "X Side Length:                %u\n",    Simulation::X_SIDE_LENGTH);
  fprintf(File,   "Y Side Length:                %u\n",    Simulation::Y_SIDE_LENGTH);
  fprintf(File,   "Z Side Length:                %u\n\n",  Simulation::Z_SIDE_LENGTH);

  // Now let's print all particle data to the file
  for(unsigned int i = 0; i < Num_Particles; i++)
    Print_Particle_To_File(Particles[i], File);

  // We've now written the 'Particle_Data' file, we can close it.
  fclose(File);
} // void Data_Dump::Print_Data_To_File(const Particle_Array & Particles) {

void Data_Dump::Print_Particle_To_File(const Particle & P_In, FILE * File) {
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

  unsigned int i;                                // index variable

  const Vector X = P_In.Get_X();
  const Vector x = P_In.Get_x();
  const Vector V = P_In.Get_V();

  // Print particle ID, dimensions
  fprintf(File,   "ID:                           %u\n",    P_In.Get_ID());
  fprintf(File,   "ijk:                          %u %u %u\n", P_In.Get_i(), P_In.Get_j(), P_In.Get_j());
  fprintf(File,   "Mass:                         %5e\n",   P_In.Get_Mass());
  fprintf(File,   "Volume:                       %5e\n",   P_In.Get_Vol());
  fprintf(File,   "Radius:                       %5lf\n",  P_In.Get_Radius());

  // Print Particle dynamic properties
  fprintf(File,   "X:                            <%6.3lf, %6.3lf, %6.3lf>\n", X(0), X(1), X(2));
  fprintf(File,   "x:                            <%6.3lf, %6.3lf, %6.3lf>\n", x(0), x(1), x(2));
  fprintf(File,   "V:                            <%6.3lf, %6.3lf, %6.3lf>\n", V(0), V(1), V(2));

  // Damage paramaters
  fprintf(File,   "Stretch_H:                    %5lf\n",  P_In.Get_Stretch_H());
  fprintf(File,   "Stretch_M:                    %5lf\n",  P_In.Get_Stretch_M());
  fprintf(File,   "Stretch_Critical:             %5lf\n",  P_In.Get_Stretch_Critical());
  fprintf(File,   "D:                            %5lf\n",  P_In.Get_D());

  // Now, let's figure out how many neighbors this particle has.
  unsigned int Num_Neighbors = P_In.Get_Num_Neighbors();

  // Neighbor paramters
  fprintf(File,   "Number of neighbors:          %u\n", P_In.Get_Num_Neighbors());

  // Print neighbor IDs
  fprintf(File,   "Neighbor IDs                  ");
  for(i = 0; i < Num_Neighbors; i++)
    fprintf(File,"%d ",P_In.Get_Neighbor_IDs(i));
  fprintf(File,"\n\n");
} // void Data_Dump::Print_Particle_To_File(const Particle & P_In, FILE * File) {

int Data_Dump::Load_Saved_Data(Particle_Array ** Array_Ptr, unsigned int & Num_Arrays) {
  // First, open up the Particle_data file
  FILE * File = fopen("../Files/Particle_Data.txt","r");

  if(File == NULL) {
    printf("Couldn't find saved data...\n");
    return 1;
  } // if(File == NULL) {

  // First, read how many Particle arrays we need to make.
  unsigned int Buf_Length = 100;
  char Buf[Buf_Length];                                 // Buffer to store text from file
  unsigned int uBuf;
  std::string strBuf;
  fread(Buf, 1, 30, File);   fscanf(File, " %u\n\n", &Num_Arrays);

  // Now use this information to genrate our arrays
  *Array_Ptr = new Particle_Array[Num_Arrays];

  // Now read in each array's name/number of particles
  for(unsigned int i = 0; i < Num_Arrays; i++) {
    /* We want to get each particle's name. The issue here is that there's no
    way to directly scan to a string. Instead, we can scan to a char array, Buf
    in this case. This copies over the contents of the name as well as a null
    terminating character into Buf. We can then assign Buf to a string buffer.
    This works by copying over the characters from Buf into the string until
    the null terminating character is encountered. To ensure that this happens
    (which may not be the case if the string is too long), we manually assign
    the final element of Buf to \0 before doing this (Otherwise, we'd get a segmenetation
    fault if the name was longer than 100 characters). We then allocate a new
    particle array and assign its name to the string. */
    fread(Buf, 1, 30, File); fscanf(File, " \"%s\"\n", Buf);
    Buf[Buf_Length-1] = '\0';
    strBuf = Buf;

    // Allocate particle array
    (*Array_Ptr)[i].Set_Name(strBuf);

    // Now read in number of particles, use it to set the number of particles in
    // the ith particle array.
    fread(Buf, 1, 30, File); fscanf(File," %u\n\n", &uBuf);
    (*Array_Ptr)[i].Set_Num_Particles(uBuf);
  } //   for(unsigned int i = 0; i < Num_Arrays; i++) {

  // CONTINUE WORKING HERE
  // pass particle arrays to read in program, have it read files using names. 

  fclose(File);
  return 0;
} // int Data_Dump::Load_Saved_Data(Particle_Array ** Array_Ptr, unsigned int & Num_Bodies) {

int Data_Dump::Load_Particle_Array_From_File(Particle_Array & Particles) {
  /* This function is designed to read in particle and use it to create a
  particles array. */

  unsigned int Num_Particles = 0;

  // First, open up the file
  FILE * File = fopen("../Files/Particle_Data.txt","r");

  if(File == NULL) {
    printf("Could not open Particle data file\n");
    return 1;
  } // if(File == NULL) {

  // Now read in the static particle members
  char Buf[100];                                 // Buffer to store text from file

  unsigned int uBuf;
  double lfBuf;
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n", &lfBuf);    Particles.Set_Inter_Particle_Spacing(lfBuf);
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n", &lfBuf);
  fread(Buf, 1, 30, File);   fscanf(File, " %u\n", &uBuf);      Particles.Set_Support_Radius(uBuf);
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n", &lfBuf);
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n", &lfBuf);    Particles.Set_Lame(lfBuf);
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n", &lfBuf);    Particles.Set_mu0(lfBuf);
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n", &lfBuf);    Particles.Set_mu(lfBuf);
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n", &lfBuf);    Particles.Set_E(lfBuf);
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n", &lfBuf);    Particles.Set_alpha(lfBuf);
  fread(Buf, 1, 30, File);   fscanf(File, " %lf\n\n", &lfBuf);  Particles.Set_Tau(lfBuf);

  // Now read in number of particles
  fread(Buf, 1, 30, File); fscanf(File, " %u\n", &Num_Particles);
  fread(Buf, 1, 30, File); fscanf(File, " %u\n", &Simulation::X_SIDE_LENGTH);
  fread(Buf, 1, 30, File); fscanf(File, " %u\n", &Simulation::Y_SIDE_LENGTH);
  fread(Buf, 1, 30, File); fscanf(File, " %u\n\n", &Simulation::Z_SIDE_LENGTH);

  // Set up Particles Particle_Array.
  Particles.Set_Num_Particles(Num_Particles);

  // Now read in particles.
  for(unsigned int i = 0; i < Num_Particles; i++)
    Load_Particle_From_File(Particles[i], File);

  //////////////////////////////////////////////////////////////////////////////
  // Now recreate each Particle's neighbor arrays.
  unsigned int Neighbor_ID;
  double V_j;                                    // Volume of jth neighbor               : mm^3
  Tensor A{0,0,0,                                // Shape Tensor (zero initialized)      : unitless Tensor
           0,0,0,
           0,0,0};
  const double Shape_Function_Amp = Particles.Get_Shape_Function_Amplitude();
  const double h = Particles.Get_h();

  for(unsigned int i = 0; i < Num_Particles; i++) {
    // Check that the current particle has Neighbors

    if(Particles[i].Num_Neighbors != 0) {
      // If so, then set up this particle's neighbor arrays.
      A = Tensor(0,0,0,
                 0,0,0,
                 0,0,0);

      for(unsigned int j = 0; j < Particles[i].Num_Neighbors; j++) {
        Neighbor_ID = Particles[i].Neighbor_IDs[j];

        // Calculate displacement vectors
        Particles[i].R[j] = Particles[Neighbor_ID].X - Particles[i].X;
        Particles[i].Mag_R[j] = Particles[i].R[j].Magnitude();

        // Calculate shape function, shape function gradient for jth neighbor
        Particles[i].W[j] = Shape_Function_Amp*(h - Particles[i].Mag_R[j])*
                            (h - Particles[i].Mag_R[j])*
                            (h - Particles[i].Mag_R[j]);
        Particles[i].Grad_W[j] = -3*Shape_Function_Amp*((h - Particles[i].Mag_R[j])*
                                 (h - Particles[i].Mag_R[j]))*
                                 (Particles[i].R[j] / Particles[i].Mag_R[j]);

        // Add in the Current Neighbor's contribution to the Shape tensor
        V_j = Particles[Neighbor_ID].Vol;
        A += Dyadic_Product((V_j*Particles[i].Grad_W[j]), Particles[i].R[j]);
      } // for(unsigned int j = 0; j < Particles[i].Num_Neighbors; i++) {

      // Now we can calculate A^(-1) from A.
      Particles[i].A_Inv = A^(-1);

      // Now that neighbors have been set, we set 'Neighbors_Are_Set' to true
      Particles[i].Neighbors_Are_Set = true;
    } // if(Particles[i].Num_Neighbors != 0) {
  } // for(unsigned int i = 0; i < Num_Particles; i++) {

  // All done, close the file.
  fclose(File);

  return 0;
} // int Data_Dump::Load_Data_From_File(unsigned int & Num_Particles, Particle_Array Particles) {

void Data_Dump::Load_Particle_From_File(Particle & P_In, FILE * File) {
  /* This function reads in the particle data for a specific particle. This
  data is transfered from the File to P_In. */

  char Buf[100];                                 // Buffer to store text from file (discarded)

  // First, read in the particle's ID and dimensions
  fread(Buf, 1, 30, File); fscanf(File, " %u\n", &P_In.ID);
  fread(Buf, 1, 30, File); fscanf(File, " %u %u %u\n", &P_In.i, &P_In.j, &P_In.k);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &P_In.Mass);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &P_In.Vol);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &P_In.Radius);

  // Now read in particle dynamic properties
  fread(Buf, 1, 30, File); fscanf(File, " <%lf, %lf, %lf>\n", &P_In.X(0), &P_In.X(1), &P_In.X(2));
  fread(Buf, 1, 30, File); fscanf(File, " <%lf, %lf, %lf>\n", &P_In.x(0), &P_In.x(1), &P_In.x(2));
  fread(Buf, 1, 30, File); fscanf(File, " <%lf, %lf, %lf>\n", &P_In.V(0), &P_In.V(1), &P_In.V(2));
  P_In.First_Time_Step = false;

  // Damage parameters
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &P_In.Stretch_H);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &P_In.Stretch_M);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &P_In.Stretch_Critical);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &P_In.D);

  // Neighbor paramaters
  P_In.Neighbors_Are_Set = false;
  fread(Buf, 1, 30, File); fscanf(File, " %u\n", &P_In.Num_Neighbors);

  // Now allocate memory for P_In's neighbor arrays
  P_In.Neighbor_IDs = new unsigned int[P_In.Num_Neighbors];
  P_In.R = new Vector[P_In.Num_Neighbors];                                     //        : mm Vector
  P_In.Mag_R = new double[P_In.Num_Neighbors];                                 //        : mm
  P_In.W = new double[P_In.Num_Neighbors];                                     //        : unitless
  P_In.Grad_W = new Vector[P_In.Num_Neighbors];                                //        : 1/mm Vector

  // Now read in neighbor IDs. Before we can do that, however, we need to move
  // the file pointer ahead, past 'Neighbor IDs: '
  fread(Buf, 1, 30, File);
  for(unsigned int i = 0; i < P_In.Num_Neighbors; i++)
    fscanf(File, " %u", &P_In.Neighbor_IDs[i]);
} // void Data_Dump::Load_Particle_From_File(Particle & P_In, FILE * File) {

#endif
