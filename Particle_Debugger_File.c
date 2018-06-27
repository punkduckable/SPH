#if !defined(PARTICLE_FILE_SOURCE)
#define PARTICLE_FILE_SOURCE

#include "Particle_Debugger_File.h"
#include "Particle.h"

void Particle_Debugger_File::Get_File_Name(string & Str) {
  char Buf[6];
  sprintf(Buf,"%05d",File_Number);
  File_Number++;

  Str += "_variables_";
  Str += Buf;
  Str += ".txt";
} // void Get_File_Name(string & Str) {

void Particle_Debugger_File::Export_Pariticle_Properties(const unsigned int Num_Particles, const Particle * Particles) {
  string File_Name = "Particle";
  Get_File_Name(File_Name);

  string File_Path = "./Particle_Files/";
  File_Path += File_Name;

  FILE * File = fopen(File_Path.c_str(), "w");

  // Print header.
  fprintf(File  , "   Particle ID   |          Internal Force          |          External Force          |          Hourglass Force\n");

  // Cycle through particles, print spacial positions of each particle
  for(unsigned int i = 0; i < Num_Particles; i++) {
    fprintf(File,"<%04.1f,%4.1f,%4.1f>: ",Particles[i].X[0], Particles[i].X[1], Particles[i].X[2]);
    fprintf(File,"<%10.5f,%10.5f,%10.5f> ",Particles[i].Force_Int[0], Particles[i].Force_Int[1], Particles[i].Force_Int[2]);
    fprintf(File,"<%10.5f,%10.5f,%10.5f> ",Particles[i].Force_Ext[0], Particles[i].Force_Ext[1], Particles[i].Force_Ext[2]);
    fprintf(File,"<%10.5f,%10.5f,%10.5f>\n",Particles[i].Force_Hg[0], Particles[i].Force_Hg[1], Particles[i].Force_Hg[2]);
  } // for(unsigned int i = 0; i < Num_Particles; i++) {

  fclose(File);
} // void Export_Pariticle_Properties(const unsigned int Num_Particles, const Particle * Particles) {

#endif
