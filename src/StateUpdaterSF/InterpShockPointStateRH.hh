// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_InterpShockPointStateRH_hh
#define ShockFitting_InterpShockPointStateRH_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include <algorithm>
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"
#include "Framework/FileLogManip.hh"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <boost/iterator/counting_iterator.hpp>
#include <fstream>
#include <vector>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/utils_classes.h>
#include <CGAL/Tetrahedron_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel_2;
typedef Kernel_2::Point_3 Point_3;
typedef CGAL::Tetrahedron_3<Kernel_2> Tetrahedron_3;



//definition of a non-mutable lvalue property map,
//with the get function as a friend function to give it
//access to the private member

class My_point_property_map{
  const std::vector<Point_3>& points;
public:
  typedef Point_3 value_type;
  typedef const value_type& reference;
  typedef std::size_t key_type;
  typedef boost::lvalue_property_map_tag category;

  My_point_property_map(const std::vector<Point_3>& pts):points(pts){}

  reference operator[](key_type k) const {return points[k];}

  friend reference get(const My_point_property_map& ppmap,key_type i)
  {return ppmap[i];}
};

typedef CGAL::Search_traits_3<Kernel_2>                                        Traits_base;
typedef CGAL::Search_traits_adapter<std::size_t,My_point_property_map,Traits_base> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits>                      K_Neighbor_search;
typedef K_Neighbor_search::Tree                                         Tree;
typedef Tree::Splitter                                                  Splitter;
typedef K_Neighbor_search::Distance                                     Distance;

//vincent 28/04/2017
//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a InterpShockPointStateRH,

class InterpShockPointStateRH {
public:

  //Constructor
  InterpShockPointStateRH();

  /// Destructor
  ~InterpShockPointStateRH();

  void InterpUpState(unsigned);

  void InterpDownState(unsigned);


private: // helper functions

  std::string getClassName() const { return "InterpShockPointStateRH"; }

  /// assign variables used in InterpShockPointStateRH to MeshData pattern
  void setMeshData();

  /// assign variables used in InterpShockPointStateRH to PhysicsData pattern
  void setPhysicsData();

  /// set starting pointers of array 2D and 3D
  void setAddress();

  // save the new interpolated upstream state
  void saveNewState(unsigned);

  /// de-allocate dynamic arrays
  void freeArray();

private: // data

  /// number of vertices
  unsigned* nvt;

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of shocks
  unsigned* nShocks;

  /// number of shock points
  std::vector<unsigned>* nShockPoints;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of mesh elements
  std::vector<unsigned>* nelem;

  /// mesh points state (assignable to MeshData)
  std::vector<double>* zroeVect;

  /// mesh points coordinates (assignable to MeshData)
  std::vector<double>* coorVect;

  /// vector characterizing nodes elements
  std::vector<int>* celnodVect;

  /// code characterizing mesh points
  std::vector <int>* nodcod;

  /// mesh points coordinates referring the shocked mesh
  Array2D <double>* XYZ;

  /// mesh points state referring to the shocked mesh
  Array2D <double>* zroe;

  /// upstream status
  Array3D <double>* ZRoeShu;

  /// downstream status
  Array3D <double>* ZRoeShd;

  /// celnod(0)(i-elem) 1° node of i-element
  /// celnod(1)(i-elem) 2° node of i-element
  /// celnod(2)(i-elem) 3° node of i-element
  Array2D <int>* celnod;

  /// shock points coordinates
  Array3D<double>* XYZSh;

  /// shock points coordinates belonging to upstream zone
  Array3D <double>* XYZShu;

  /// shock points coordinates belonging to downstream zone
  Array3D <double>* XYZShd;

  /// file storing info
     FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_InterpShockPointStateRH_hh
