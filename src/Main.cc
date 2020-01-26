//#define NDEBUG
#include "Simulation/Simulation.h"
#include "Errors.h"

int main(void) {
  try {
    Simulation::Load_Setup_File();
  } // try {
  catch(const Exception & Error_In) {
    printf("%s\n", Error_In.what());
  } // catch(Exception & Error_In) {

  return 0;
} // catch.hpp
