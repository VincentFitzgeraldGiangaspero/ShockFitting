// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_Interp3D_hh
#define ShockFitting_Interp3D_hh

//--------------------------------------------------------------------------//

#include <algorithm>
#include <cmath>
#include "Framework/StateUpdater.hh"
#include "Framework/FileLogManip.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"


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



//--------------------------------------------------------------------------//

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
//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a Interp3D, whose task is to update values in the phantom
/// nodes of the background mesh (referred to index (0)) using values in the
/// shocked mesh (referred to index (1))

class Interp3D : public StateUpdater {
public:

  /// Constructor
  /// @param objectName the concrete class name
  Interp3D(const std::string& objectName);

  /// Destructor
  virtual ~Interp3D();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// update phantom nodes values
  virtual void update();


private: // helper functions

  /// return class name
  std::string getClassName() const {return "Interp3D";}

  /// assign variables used in Interp3D to MeshData pattern
  void setMeshData();

  /// assign variables used in Interp3D to PhysicsData pattern
  void setPhysicsData();

  /// set starting pointers of array 2D and 3D
  void setAddress();

  /// set vectors for CGAL routine
  void setCGALvectors();

  /// find the cell which the phantom node belongs to
  void finder(unsigned);

  /// return cell which the phantom node belongs to
  unsigned getCell() const { return ielem; }

  /// return boolean variable checking if the phantom node is found
  unsigned getIfound() const { return ifound; }

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

  /// map vector
  std::vector<int>* M02M1;

  /// map vector -> M02M1+NPOIN(0)
  std::vector<int>* M02M12;

  /// mesh points coordinates referring the background mesh
  Array2D <double>* XYZBkg;

  /// mesh points state referring to the background mesh
  Array2D <double>* zBkg;

  /// mesh points coordinates referring the shocked mesh
  Array2D <double>* XYZ;

  /// mesh points state referring to the shocked mesh
  Array2D <double>* zroe;

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

  /// i-cell including the phantom node
  unsigned ielem;

  /// shock faces
  Array3D <unsigned>* ShFaces;

  /// number of shock faces for each shock
  std::vector <unsigned>* nShockFaces;

  // cgal vector
  std::vector<Point_3> mesh_points;

  /// variable checking if the phantom node is been found
  /// @param ifound = 1 phantom node not found
  /// @param ifound = 0 phantom node found
  unsigned ifound;

  /// store log file infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_Interp3D_hh

