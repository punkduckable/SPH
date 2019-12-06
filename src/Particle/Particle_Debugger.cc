#include "Particle_Debugger.h"

std::string Particle_Debugger::Get_File_Name(const std::string & Str) {
  // First, figure out if we've seen this string before. if so, then it'll be
  // in the name list
  unsigned i;

  unsigned Num_Names = Name_List.Node_Count();
  for(i = 0; i < Num_Names; i++) {
    // Compare the ith name in the Name List to the supplied string.
    if(Name_List[i].compare(Str) == 0)
      break;
  } // for(unsigned i = 0; i <= Num_Names; i++) {

  // If the supplied string is NOT in the list, then add Str to the end of the
  // name list and add a zero node to the File_Number list.
  if(i >= Num_Names) {
    Name_List.Add_Back(Str);
    File_Number_List.Add_Back(0);
  } // if(i > Num_Names) {

  // Now print the appropiate file number to the buffer, increment the file number.
  char Buf[6];
  sprintf(Buf,"%05u",File_Number_List[i]);
  File_Number_List[i]++;
  std::string File_Name = Str;

  // Now append file name syntax to file name, return.
  File_Name += "_Forces_";
  File_Name += Buf;
  File_Name += ".txt";

  return File_Name;
} // std::string Get_File_Name(cosnt std::string & Str) {



void Particle_Debugger::Export_Particle_Forces(const Body & Particles) {
  // First, figure out the number of particles
  const unsigned Num_Particles = Particles.Get_Num_Particles();

  // Now create a file path for the new file (based on the Body's name
  // and how many times we've printed this Body)
  std::string File_Path = "../Files/Force_Files/";
  File_Path += Get_File_Name(Particles.Get_Name());
  FILE * File = fopen(File_Path.c_str(), "w");

  // Print header.
  fprintf(File,"  ID  |");
  fprintf(File," Particle Pos  |");
  fprintf(File,"        Internal Force        |");

  #if defined(PARTICLE_DEBUG)
    fprintf(File,"        Viscous Force         |");
  #endif

  fprintf(File,"        Contact Force         |");
  fprintf(File,"        Friction Force        |");
  fprintf(File,"        Hourglass Force       |");
  fprintf(File,"\n");

  // Cycle through particles, print spacial positions of each particle
  for(unsigned i = 0; i < Num_Particles; i++) {
    fprintf(File,"%6u|", Particles[i].Get_ID());
    fprintf(File,"%4.1f,%4.1f,%4.1f | ",    Particles[i].X[0],            Particles[i].X[1],            Particles[i].X[2]);
    fprintf(File,"<%8.1e,%8.1e,%8.1e> | ",  Particles[i].Force_Int[0],    Particles[i].Force_Int[1],    Particles[i].Force_Int[2]);

    #if defined(PARTICLE_DEBUG)
      fprintf(File,"<%8.1e,%8.1e,%8.1e> | ",  Particles[i].Force_Visc[0],   Particles[i].Force_Visc[1],   Particles[i].Force_Visc[2]);
    #endif

    fprintf(File,"<%8.1e,%8.1e,%8.1e> | ",  Particles[i].Force_Contact[0],Particles[i].Force_Contact[1],Particles[i].Force_Contact[2]);
    fprintf(File,"<%8.1e,%8.1e,%8.1e> | ",  Particles[i].Force_Friction[0],Particles[i].Force_Friction[1],Particles[i].Force_Friction[2]);
    fprintf(File,"<%8.1e,%8.1e,%8.1e>\n",   Particles[i].Force_HG[0],     Particles[i].Force_HG[1],     Particles[i].Force_HG[2]);
  } // for(unsigned i = 0; i < Num_Particles; i++) {

  fclose(File);
} // void Particle_Debugger::Export_Particle_Forces(const Body & Particles) {
