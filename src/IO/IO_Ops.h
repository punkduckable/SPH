#if !defined(IO_OPS_HEADER)
#define IO_OPS_HEADER

#include <string>
#include <stdio.h>
#include <vector>
#include "Vector/Vector.h"

/* File description:
This file houses several functions that make reading and writing from a file
easier. */

std::string read_line_after_char(FILE * File, const char delim);
Vector read_box_BC(std::string BC_Str);
//void write_box_BC(FILE * File, const Vector & BC_In);

namespace String_Ops {
  bool Contains(const char* Buffer,                                            // Intent: Read
                const char* Word,                                              // Intent: Read
                unsigned Search_At = 0);                                       // Intent: Read

  std::vector<std::string> Split(std::string & S,                              // Intent: Read
                                 const char Delim = ',');                      // Intent: Read
  std::vector<std::string> Split(const char* S,                                // Intent: Read
                                 const char Delim = ',');                      // Intent: Read
} // namespace String_Ops {

#endif
