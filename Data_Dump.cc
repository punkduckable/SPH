#if !defined(DATA_DUMP_SOURCE)
#define DATA_DUMP_SOURCE

void Data_Dump::Load_From_File(Particle * Particles, unsigned int & Num_Particles) {

}

void Data_Dump::Print_Data_To_File(const Particle * Particles, const unsigned int Num_Particles) {
  /* This function prints all data needed to reconstruct the Particle's array
  to a file. The intent of this is to allow the user to essentially 'save' the
  state of the program, allowing the user to - at a future point - pick up where
  they last saved/left off. */

  /* First, let's create the file. Note, if another 'Particle_Data' file already
  exists then this overwrites/deletes that file. */
  FILE * File = fopen("./Particle_Data.txt", "w");

  // Let's begin by printing the 'static' particle class paramaters
  fprintf(File,   "Inter Particle Spacing:       %5f\n",Particle::Inter_Particle_Spacing);
  fprintf(File,   "Support radius (h):           %5f\n",Particle::h);
  fprintf(File,   "Shape Function Amplitude:      %5f\n",Particle::Shape_Function_Amp);
  fprintf(File,   "Lame parameter:               %5f\n",Particle::Lame);
  fprintf(File,   "Shear modulus (mu0):          %5f\n",Particle::mu0);
  fprintf(File,   "Viscosity (mu):               %5f\n",Particle::Inter_Particle_Spacing);
  fprintf(File,   "Hourglass Stiffness (E):      %5f\n",Particle::E);
  fprintf(File,   "alpha (HG parameter):         %5f\n",Particle::alpha);
  fprintf(File,   "Tau (damage parameter):       %5f\n",Particle::Tau);

  // Now let's print the number of particles
  fprintf(File,   "\n");
  fprintf(File,   "Number of particles:          %u\n\n",Num_Particles);

  // Now let's print all particle data to the file
  for(unsigned int i = 0; i < Num_Particles; i++)
    Print_Particle_To_File(Particles[i], File);

  // We've now written the 'Particle_Data' file, we can close it.
  fclose(File);
} // void Data_Dump::Print_Data_To_File(const Particle * Particles, const unsigned int Num_Particles) {

void Data_Dump::Print_Particle_To_File(const Particle & P_In, FILE * File) {
  /* This function prints all the information that is needed to re-create the
  input particle. It should be noted that we do NOT print all neighbor arrays.
  this is because these arrays can be recalculated if we know the neighbor IDs.
  Not storing this information in the File makes the file take up less/easier
  to read. This function assumes that the File has already been setup (with
  static particle class paramaters). */

  unsigned int i;                                // index variable

  // Print particle ID, dimensions
  fprintf(File,   "ijk:                          %u %u %u\n", P_In.i, P_In.j, P_In.k);
  fprintf(File,   "ID:                           %u\n", P_In.ID);
  fprintf(File,   "Mass:                         %5e\n", P_In.Mass);
  fprintf(File,   "Volume:                       %5e\n", P_In.Vol);
  fprintf(File,   "Radius:                       %5f\n", P_In.Radius);

  // Print Particle dynamic properties
  fprintf(File,   "X:                            <%6.3f, %6.3f, %6.3f>\n", P_In.X(0), P_In.X(1), P_In.X(2));
  fprintf(File,   "x:                            <%6.3f, %6.3f, %6.3f>\n", P_In.x(0), P_In.x(1), P_In.x(2));
  fprintf(File,   "V:                            <%6.3f, %6.3f, %6.3f>\n", P_In.V(0), P_In.V(1), P_In.V(2));
  fprintf(File,   "First time step:              %i\n",P_In.First_Time_Step);
  fprintf(File,   "PK stress tensor:             | %6.3f %6.3f %6.3f |\n", P_In.P(0,0), P_In.P(0,1), P_In.P(0,2));
  fprintf(File,   "                              | %6.3f %6.3f %6.3f |\n", P_In.P(1,0), P_In.P(1,1), P_In.P(1,2));
  fprintf(File,   "                              | %6.3f %6.3f %6.3f |\n", P_In.P(2,0), P_In.P(2,1), P_In.P(2,2));
  fprintf(File,   "Deformation gradient (F):     | %6.3f %6.3f %6.3f |\n", P_In.P(0,0), P_In.P(0,1), P_In.P(0,2));
  fprintf(File,   "                              | %6.3f %6.3f %6.3f |\n", P_In.P(1,0), P_In.P(1,1), P_In.P(1,2));
  fprintf(File,   "                              | %6.3f %6.3f %6.3f |\n", P_In.P(2,0), P_In.P(2,1), P_In.P(2,2));

  // Forces
  fprintf(File,   "Internal Force:               <%6.3f, %6.3f, %6.3f>\n", P_In.Force_Int(0), P_In.Force_Int(1), P_In.Force_Int(2));
  fprintf(File,   "Contact Force:                <%6.3f, %6.3f, %6.3f>\n", P_In.Force_Contact(0), P_In.Force_Contact(1), P_In.Force_Contact(2));
  fprintf(File,   "HourGlass Force:              <%6.3f, %6.3f, %6.3f>\n", P_In.Force_HG(0), P_In.Force_HG(1), P_In.Force_HG(2));

  // Damage paramaters
  fprintf(File,   "Stretch_H:                    %5f\n", P_In.Stretch_H);
  fprintf(File,   "Stretch_M:                    %5f\n", P_In.Stretch_M);
  fprintf(File,   "Stretch_Critical:             %5f\n", P_In.Stretch_Critical);
  fprintf(File,   "D:                            %5f\n", P_In.D);

  // Now, let's figure out how many neighbors this particle has.
  unsigned int Num_Neighbors = P_In.Num_Neighbors;

  // Print neighbor paramters
  fprintf(File,   "Has Neighbors:                %u\n", P_In.Has_Neighbors);
  fprintf(File,   "Number of neighbors:          %u\n", P_In.Num_Neighbors);

  // Print neighbor IDs
  fprintf(File,   "Neighbor IDs                  ");
  for(i = 0; i < Num_Neighbors; i++)
    fprintf(File,"%d ",P_In.Neighbor_IDs[i]);
  fprintf(File,"\n\n");
} // void Data_Dump::Print_Particle_To_File(const Particle & P_In, FILE * File) {

void Data_Dump::Load_Particle_From_File(Particle * Particles, const FILE * File) {

}

#endif
