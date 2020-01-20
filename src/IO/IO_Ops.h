#if !defined(IO_OPS_HEADER)
#define IO_OPS_HEADER

#include <string>
#include <stdio.h>
#include <vector>
#include "Classes.h"

/* File description:
This file houses several functions that make reading and writing from a file
easier. */

std::string read_line_after_char(FILE * File, const char delim);

namespace String_Ops {
  bool Contains(const char* Buffer,                                            // Intent: Read
                const char* Word,                                              // Intent: Read
                unsigned Search_At = 0);                                       // Intent: Read
} // namespace String_Ops {

#endif
