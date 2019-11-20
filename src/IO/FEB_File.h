#if !defined(FEB_HEADER)
#define FEB_HEADER

namespace FEB_File {
  using std::string;

  int Read_FEB_File(const string & File_Name, Vector ** X_Ptr, unsigned int & Num_Nodes);
}

#endif
