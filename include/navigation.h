#include <cstdlib>
#include <cmath>
#include "par_io.h"

#ifndef _NAVIGATION_H_
#define _NAVIGATION_H_

//This Class handles the navigation on each timeslice. Two persistent 2d-arrays are
//allocated with the number of 3d points on the lattice. These lookup tables are 3dimensional,
//one holds the indices in the upwards, one in the downwards direction. Every
//spacepoint is encoded in an integer respecting periodic boundary conditions.
//Map from 3d-spatial coordinates to lookup table integers is:
//ind(x,y,z) = x*L2*L1 + y*L1 + z
//second dimension of lookup tables gives index in x- (0), y- (1) and z- (2)
//direction.  
//typedef boost::multi_array<int,2> look;
typedef int** look;
class Nav {

  private:
    //data members
    look iup;
    look idown;
    int vol;
    // Constructor taking a pointer to a parameters object as argument. The
    // corresponding variables are used to construct the lookup tables.
    Nav();
    //copy assignment
    ~Nav();
  protected:

  public:
    
    //Get index in distance. org stands for the origin, dir for the direction,
    //dist for the distance in lattice units.
    static Nav* getInstance();
    void init();
    int dist_ind(const int org, const int dir, const int dist);
    int get_up(const int pos, const int dir);
    int get_dn(const int pos, const int dir);

};

#endif //_NAVIGATION_H_
