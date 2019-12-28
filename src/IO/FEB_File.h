#if !defined(FEB_HEADER)
#define FEB_HEADER

#include "Classes.h"
#include <string>

namespace FEB_File {
  using std::string;

  int Read_FEB_File(const string & File_Name, Vector ** X_Ptr, unsigned & Num_Nodes);
}

#endif
