#if !defined(ARRAY_HEADER)
#define ARRAY_HEADER

#include "List.h"
#include "Errors.h"
#include <assert.h>

template <typename T>
class Array {
  private:
    T* Internal_Array;
    unsigned Length;
    bool Length_Set;

    void Private_Setup_From_List(List<T> & List_In) {
      /* This function exists to reduce code redundancy.
      the List constructor and "Setup_From_List" function basically do the same
      thing; most of their code is the same. Thus, to eliminate repeating the
      same code twice, I made this function (both of the aformentioned functions
      call this one). */

      unsigned Num_Nodes = List_In.Get_Num_Nodes();
      (*this).Length = Num_Nodes;
      (*this).Internal_Array = new T[Num_Nodes];
      (*this).Length_Set = true;

      for(unsigned i = 0; i < Num_Nodes; i++) {
        (*this).Internal_Array[i] = List_In.Pop_Front();
      } // for(unsigned i = 0; i < Num_Nodes; i++)
    }

  public:
    ////////////////////////////////////////////////////////////////////////////
    // Constructors, destructor.
    Array(void) : Internal_Array(nullptr), Length(0), Length_Set(false) {} // Default constructor

    Array(const Array & Ar_In) = delete;         // Copy constructor
    Array(const unsigned Length) : Internal_Array(new T[Length]), Length(Length), Length_Set(true) {}
    Array(List<T> & List_In) {                   // Make Array from List.
      /* This function exists to convert Lists into Arrays.
      The array has the same length as the input list. The list is destroyed in
      the process. */
      Private_Setup_From_List(List_In);
    } // Array(List<T> & List_In) {
    ~Array(void) { delete [] Internal_Array; }



    ////////////////////////////////////////////////////////////////////////////
    // Operator Overloading

    Array & operator=(const Array & Ar_In) = delete;

    T operator[](const unsigned index) const {
      if((*this).Length_Set == false) {
        throw Array_Length_Not_Set("Array Length Not Set Exception: Thrown by Array::operator[] (access)\n"
                                   "An array must have a length before you can access its elements\n");
      } // if((*this).Length_Set == false) {
      assert(index < (*this).Length);

      return Internal_Array[index];
    } // T operator[](const unsigned index) const {

    T & operator[](const unsigned index) {
      if((*this).Length_Set == false) {
        throw Array_Length_Not_Set("Array Length Not Set Exception: Thrown by Array::operator[] (modify)\n"
                                   "An array must have a length before you can access its elements\n");
      } // if((*this).Length_Set == false) {
      assert(index < (*this).Length);

      return Internal_Array[index];
    } // T & operator[](const unsigned index) {



    ////////////////////////////////////////////////////////////////////////////
    // Getters

    unsigned Get_Length(void) const { return (*this).Length; }
    bool Get_Length_Set(void) const { return (*this).Length_Set; }



    ////////////////////////////////////////////////////////////////////////////
    // Setters

    void Set_Lenght(const unsigned Length_In) {
      (*this).Length = Length_In;
      (*this).Length_Set = true;
    } // void Set_Lenght(const unsigned Length_In) {



    ////////////////////////////////////////////////////////////////////////////
    // Other methods

    void Setup_From_List(List<T> & List_In) {
      /* This function exists to convert a List into an Array.
      The array has the same length as the input list. The list is destroyed in
      the process. */
      if((*this).Length_Set == true) {
        throw Array_Length_Already_Set("List Length Already Set Exception: thrown by Array::Setup_From_List\n"
                                       "An Array can not be setup from a list if its length has already been set\n");
      } // if((*this).Length_Set == true) {

      (*this).Private_Setup_From_List(List_In);
    } // void Setup_From_List(List<T> & List_In) {
}; // class Array {

#endif
