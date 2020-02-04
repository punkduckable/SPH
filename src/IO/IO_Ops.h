#if !defined(IO_OPS_HEADER)
#define IO_OPS_HEADER

#include <string>
#include <stdio.h>
#include <vector>
#include <fstream>
#include "Classes.h"

/* File description:
This file houses several functions that make reading and writing from a file
easier. */

namespace IO  {
  std::string read_line_after(std::ifstream & File,
                              const char* Word,
                              bool Case_Sensitive = true);

  namespace String_Ops {
    bool Contains(const char* Buffer,
                  const char* Word,
                  bool Case_Sensitive = true,
                  unsigned Start_At = 0);

    int Index_After_Word(const char* Buffer,
                         const char* Word,
                         bool Case_Sensitive = true,
                         unsigned Start_At = 0);
  } // namespace String_Ops {
}

#endif
