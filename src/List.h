#if !defined(LIST_HEADER)
#define LIST_HEADER

#include "classes.h"
#include "Errors.h"
#include <stdio.h>

template <typename T>
class List {
  private:
    struct Node {
      T Value;
      Node *Next_Node;
      Node *Previous_Node;
    }; // struct Node {

    Node* Start;                                 // Start of the list
    Node* End;                                   // End of the list
    unsigned Num_Nodes;                          // Number of nodes in the list

  public:
    // Constructors, destructor
    List();                                      // Default constructor
    List(const List & L_In) = delete;            // Copy constructor (to prevent inadvertent shallow copying)
    ~List();                                     // Destructor

    // Operator overloading
    List & operator=(const List &) = delete;     // Explicit = operator overload (to prevent inadvertent shallow copying)

    // Methods to modify the list.
    void Push_Back(const T ID_In);               // Add an item to the start of the list
    void Push_Front(const T ID_In);              // Add an item to the end of the list
    T Pop_Back(void);                            // Remove an item from the end of the list
    T Pop_Front(void);                           // Add an item to the start of the list
    T & Front(void);                             // Returns a reference to the value of the first item of the list.
    T & Back(void);                              // Returns a reference to the value of the last item of the list.

    // Methods to get list info
    unsigned Get_Length(void) const;             // Returns number of nodes
    void Print_Node_Info(void) const;            // For testing
};



////////////////////////////////////////////////////////////////////////////////
// Constructors, destructor
template <typename T>
List<T>::List(void) {
  Num_Nodes = 0;
  Start = nullptr;
  End = nullptr;
} // List<T>::List() {



template <typename T>
List<T>::~List(void) {
  /* To delete this node, we need to free all memory that has bee nallocated by
  our list. To do this, we need to delete our nodes one by one. */

  Node * temp;
  while(Start != nullptr) {
    temp = Start;                                // Store current Start in temp (so we can free it!)
    Start = Start->Next_Node;                    // Move Start forward one node
    delete temp;                                 // free the old Start node

    Num_Nodes--;
    //printf("%d nodes remaining\n",Num_Nodes);
  } // while(Start != nullptr) {
} // List<T>::~List(void) {





////////////////////////////////////////////////////////////////////////////////
// Method that modify the list.

template <typename T>
void List<T>::Push_Back(const T ID_In) {
  /* This method is used to add a node onto the end of our list. We do this by
  dynamically allocating a node, adding it into our list (old End, if there is
  one, points to the new node, make End point to new node). */
  Node * add = new Node;                       // Dynamically allocate the new node
  add->Value = ID_In;                    // New Node's ID is the input ID
  add->Next_Node = nullptr;                       // New node is End, so there is no next node

  /* Check if there are any nodes in our list.

  If the list is empty (no nodes), then the added element will be both the
  Start and the End. However, this means that there is no previous node.

  If our list already has elements, then add is the new End. The old End
  is the Previous_Node for add.*/
  if(Num_Nodes == 0) {
    Start = add;
    add->Previous_Node = nullptr;                // Add is new Start (and End), so there is no previous node.
  } // if(Num_Nodes == 0) {
  else {
    add->Previous_Node = End;                    // add is new End. Previous node previous node for add (new end)
    End->Next_Node = add;                        // old end's new node is add (new end)
  } // else {

  End = add;                                     // add is new End
  Num_Nodes++;                                   // We added a new node, incremenet number of nodes
} // void List<T>::Push_Back(const T ID_In) {



template <typename T>
void List<T>::Push_Front(const T ID_In) {
  /* this method adds a node onto the Start of our list. We do this by
  dynamically allocating a node and adding it to the Start of our list (this
  new node will become the new Start/will point to the old Start) */

  Node *add = new Node;                          // Dynamically allocate new node
  add->Value = ID_In;                            // The new node's value is the input ID
  add->Previous_Node = nullptr;                  // add will be new Start, so there is no node before add.

  /* Check if there are any nodes in our list:

  If the list is empty, then add is both the Start and the End. This means that
  there is no next node to point to.

  If the list is not empty, then add is the new Start. The old Start is the
  Next_Node for add */
  if(Num_Nodes == 0) {
    End = add;
    add->Next_Node = nullptr;
  } // if(Num_Nodes == 0) {
  else {
    add->Next_Node = Start;                      // Old Start is next_node for add (new Start).
    Start->Previous_Node = add;                  // add (new Start) is previous node of old Start.
  } // else {

  Start = add;                                   // add is new Start
  Num_Nodes++;                                   // Increment number of nodes
} // void List<T>::Push_Front(const T ID_In) {



template <typename T>
T List<T>::Pop_Back(void) {
  /* This method removes the last node on our list (if there is one). This is
  done by having End point to the Old End's Previous node and then freeing the
  old End. */

  // Check if there is a node to remove.
  if((*this).Num_Nodes == 0) {
    throw Empty_List("Empty List Exception: Thrown by List<T>::Pop_Back\n"
                     "This list is empty. You can not pop an item from an empty list!\n");
  } // if((*this).Num_Nodes == 0) {

  T ID_Out = End->Value;

  Node * temp = End;                  // temp node points to old End
  End = temp->Previous_Node;          // Move End backward (to old End's previous node)
  delete temp;                         // free memory allocated by old End

  Num_Nodes--;                         // decrement number of nodes
  //printf("Node removed from end. %d nodes remaining\n",Num_Nodes);      // For testing

  /* check if there are any nodes left in the list. If there are none, then
  Start should point to nullptr */
  if(Num_Nodes == 0) { Start = nullptr; }

  return ID_Out;
} // T List<T>::Pop_Back(void) {



template <typename T>
T List<T>::Pop_Front(void) {
  /* This method removes the first node on our list (if there is one). This is
  done by having Start point to the Old Start's Next_Node, and then freeing the
  old Start. */

  // Check if there is a node to remove
  if((*this).Num_Nodes == 0) {
    throw Empty_List("Empty List Exception: Thrown by List<T>::Pop_Front\n"
                     "This list is empty. You can not pop an item from an empty list!\n");
  } // if((*this).Num_Nodes == 0)

  T ID_Out = Start->Value;

  Node * temp = Start;                           // temp points to old Start.
  Start = temp->Next_Node;                       // Move Start forward
  delete temp;                                   // free memory allocated by old Start

  Num_Nodes--;                                   // decrement number of nodes
  //printf("Node remove from Start. %d nodes remaining\n",Num_Nodes);       // For testing

  /* Check if there are any nodes left in the list. If there are none, then End
  should point to nullptr. */
  if(Num_Nodes == 0) { End = nullptr; }

  return ID_Out;
} // T List<t>::Pop_Front(void) {



template <typename T>
T & List<T>::Front(void) {
  if((*this).Num_Nodes == 0) {
    throw Empty_List("Empty List Exception: Thrown by List<T>::Front\n"
                     "This list is empty. You can not pop an item from an empty list!\n");
  } // if((*this).Num_Nodes == 0) {

  return (*this).Start->Value;
} // T & List<T>::Front(void) {



template <typename T>
T & List<T>::Back(void) {
  if((*this).Num_Nodes == 0) {
    throw Empty_List("Empty List Exception: Thrown by List<T>::Back\n"
                     "This list is empty. You can not pop an item from an empty list!\n");
  } // if((*this).Num_Nodes == 0) {

  return (*this).End->Value;
} // T & List<T>::Back(void) {






////////////////////////////////////////////////////////////////////////////////
// Methods to get list info

template <typename T>
unsigned int List<T>::Get_Length(void) const {
  return Num_Nodes;
} // void List<T>::Get_Length(void) const {



template <typename T>
void List<T>::Print_Node_Info(void) const {
  Node * Current_Node = Start;
  int Node_Num = 0;

  // Print info on list
  printf("Num Nodes: %d\n",Num_Nodes);
  printf("Start:     %p\n",Start);
  printf("End:       %p\n\n",End);

  while(Current_Node != nullptr) {
    // Print info on current node
    printf("Node #%d\n",Node_Num);
    printf("Value:     %d\n",Current_Node->Value);
    printf("Current    %p\n",Current_Node);
    printf("Previous:  %p\n",Current_Node->Previous_Node);
    printf("Next:      %p\n\n",Current_Node->Next_Node);

    // Go to next node, cycle node count
    Current_Node = Current_Node->Next_Node;
    Node_Num++;
  } //   while(Current_Node != nullptr) {
} // void List<T>::Print_Node_Info(void) const {

#endif
