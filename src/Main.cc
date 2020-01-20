//#define NDEBUG
#include "Simulation/Simulation.h"
#include "Errors.h"

int main(void) {
  try {
    Simulation::Run_Simulation();
  } // try {
  catch(const Cant_Open_File & Error_In) {
    printf("%s\n", Error_In.what());
  } // catch(Exception & Error_In) {

  return 0;
} // catch.hpp
