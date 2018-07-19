#if !defined(DATA_DUMP_SOURCE)
#define DATA_DUMP_SOURCE

void Data_Dump::Print_Data_To_File(const Particle * Particles, const unsigned int Num_Particles) {
  /* This function prints all data needed to reconstruct the Particle's array
  to a file. The intent of this is to allow the user to essentially 'save' the
  state of the program, allowing the user to - at a future point - pick up where
  they last saved/left off. */

  /* First, let's create the file. Note, if another 'Particle_Data' file already
  exists then this overwrites/deletes that file. */
  FILE * File = fopen("./Particle_Data.txt", "w");

  // Let's begin by printing the 'static' particle class paramaters
  fprintf(File,   "Inter Particle Spacing:       %5lf\n",Particle::Inter_Particle_Spacing);
  fprintf(File,   "Support Radius (mm):          %5lf\n",Particle::h);
  fprintf(File,   "Support Radius (IPS):         %u\n",Particle::Support_Radius);
  fprintf(File,   "Shape Function Amplitude:     %5lf\n",Particle::Shape_Function_Amp);
  fprintf(File,   "Lame parameter:               %5lf\n",Particle::Lame);
  fprintf(File,   "Shear modulus (mu0):          %5lf\n",Particle::mu0);
  fprintf(File,   "Viscosity (mu):               %5lf\n",Particle::mu);
  fprintf(File,   "Hourglass Stiffness (E):      %5lf\n",Particle::E);
  fprintf(File,   "alpha (HG parameter):         %5lf\n",Particle::alpha);
  fprintf(File,   "Tau (damage parameter):       %5lf\n",Particle::Tau);

  // Now let's print the number of particles
  fprintf(File,   "\n");
  fprintf(File,   "Number of particles:          %u\n",Num_Particles);

  // Finally, let's print the cuboid paramaters (should be removed if not using
  // a cuboid)
  fprintf(File,   "X Side Length:                %u\n", Simulation::X_SIDE_LENGTH);
  fprintf(File,   "Y Side Length:                %u\n", Simulation::Y_SIDE_LENGTH);
  fprintf(File,   "Z Side Length:                %u\n\n", Simulation::Z_SIDE_LENGTH);

  // Now let's print all particle data to the file
  for(unsigned int i = 0; i < Num_Particles; i++)
    Print_Particle_To_File(Particles[i], File);

  // We've now written the 'Particle_Data' file, we can close it.
  fclose(File);
} // void Data_Dump::Print_Data_To_File(const Particle * Particles, const unsigned int Num_Particles) {

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

  // Print particle ID, dimensions
  fprintf(File,   "ID:                           %u\n", P_In.ID);
  fprintf(File,   "ijk:                          %u %u %u\n", P_In.i, P_In.j, P_In.k);
  fprintf(File,   "Mass:                         %5e\n", P_In.Mass);
  fprintf(File,   "Volume:                       %5e\n", P_In.Vol);
  fprintf(File,   "Radius:                       %5lf\n", P_In.Radius);

  // Print Particle dynamic properties
  fprintf(File,   "X:                            <%6.3lf, %6.3lf, %6.3lf>\n", P_In.X(0), P_In.X(1), P_In.X(2));
  fprintf(File,   "x:                            <%6.3lf, %6.3lf, %6.3lf>\n", P_In.x(0), P_In.x(1), P_In.x(2));
  fprintf(File,   "V:                            <%6.3lf, %6.3lf, %6.3lf>\n", P_In.V(0), P_In.V(1), P_In.V(2));

  // Damage paramaters
  fprintf(File,   "Stretch_H:                    %5lf\n", P_In.Stretch_H);
  fprintf(File,   "Stretch_M:                    %5lf\n", P_In.Stretch_M);
  fprintf(File,   "Stretch_Critical:             %5lf\n", P_In.Stretch_Critical);
  fprintf(File,   "D:                            %5lf\n", P_In.D);

  // Now, let's figure out how many neighbors this particle has.
  unsigned int Num_Neighbors = P_In.Num_Neighbors;

  // Neighbor paramters
  fprintf(File,   "Number of neighbors:          %u\n", P_In.Num_Neighbors);

  // Print neighbor IDs
  fprintf(File,   "Neighbor IDs                  ");
  for(i = 0; i < Num_Neighbors; i++)
    fprintf(File,"%d ",P_In.Neighbor_IDs[i]);
  fprintf(File,"\n\n");
} // void Data_Dump::Print_Particle_To_File(const Particle & P_In, FILE * File) {

Particle * Data_Dump::Load_Data_From_File(unsigned int & Num_Particles) {
  /* This function is designed to read in particle and use it to create a
  particles array. */

  // First, open up the file
  FILE * File = fopen("./Particle_Data.txt","r");

  // Now read in the static particle members
  char Buf[100];                                 // Buffer to store text from file

  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &Particle::Inter_Particle_Spacing);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &Particle::h);
  fread(Buf, 1, 30, File); fscanf(File, " %u\n", &Particle::Support_Radius);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &Particle::Shape_Function_Amp);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &Particle::Lame);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &Particle::mu0);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &Particle::mu);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &Particle::E);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n", &Particle::alpha);
  fread(Buf, 1, 30, File); fscanf(File, " %lf\n\n",&Particle::Tau);

  // Now read in number of particles
  fread(Buf, 1, 30, File); fscanf(File, " %u\n", &Num_Particles);
  fread(Buf, 1, 30, File); fscanf(File, " %u\n", &Simulation::X_SIDE_LENGTH);
  fread(Buf, 1, 30, File); fscanf(File, " %u\n", &Simulation::Y_SIDE_LENGTH);
  fread(Buf, 1, 30, File); fscanf(File, " %u\n\n", &Simulation::Z_SIDE_LENGTH);


  // Use this to allocate the particle's array
  Particle * Particles = new Particle[Num_Particles];

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
  const double Shape_Function_Amp = Particle::Shape_Function_Amp;
  const double h = Particle::h;

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

  return Particles;
} // Particle * Data_Dump::Load_Data_From_File(unsigned int & Num_Particles) {

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
