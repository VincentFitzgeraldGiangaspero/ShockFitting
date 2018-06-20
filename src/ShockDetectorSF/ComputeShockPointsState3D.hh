// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ComputeShockPointsState3D_hh
#define ShockFitting_ComputeShockPointsState3D_hh

//--------------------------------------------------------------------------//

#include <algorithm>
#include <cmath>
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

/// This class define a ComputeShockPointsState3D, whose task is to evaluate
/// the downstream and upstream state of the shock points by interpolating
/// them in the directions of the shock points normal vectors

/// @author Vincent

/// vincent 31/01/17

class ComputeShockPointsState3D {
public:

  /// Constructor
  ComputeShockPointsState3D();

  /// Destructor
  ~ComputeShockPointsState3D();

  /// setup all the allocated vectors and arrays
  /// @param shockLayerThickess    distance used to extract the upstream
  ///                              and downstream coordinates
  void setup(double shockLayerThickess);

  /// unsetup all the allocated vectors and arrays
  void unsetup();

  /// extract upstream and downstream points even if the downstream
  /// upstream state will be evaluated after
  /// @param normals   array storing normal vector to the shock points
  void extractDownstreamAndUpstreamPoints(Array3D<double> normals);

  /// evaluate the downstream and upstream state for each shock point
  /// by interpolating the state of the surrounding nodes
  void interpDownstreamAndUpstreamState();

  /// assign upstream and downstream state according to the
  /// interpolated values
  void assignDownstreamAndUpstreamState();

private: // helper functions

   /// return class name
   std::string getClassName() const { return "ComputeShockPointsState3D"; }

   /// set the shock layer thickness
   void setShockLayerThickness(double shLayerThick)
    { m_shockLayerThickness = shLayerThick; }

   /// resize vectors and array
   void setSize();

   /// assign strating pointers to array
   void setAddress();

   /// assign variables used in ComputeShockPointsState3D to MeshData
   void setMeshData();

   /// assign variables used in ComputeShockPointsState3D to PhysicsData
   void setPhysicsData();

   /// de-allocate dynamic array
   void freeArray();

private: // data

   /// distance used to extract the upstream and downstream coordinates
   double m_shockLayerThickness;

   /// number of chemical species
   unsigned* nsp;

   /// number of vertices for each mesh element
   unsigned* nvt;

   /// number of degrees of freedom
   unsigned* ndof;

   /// number of shock points
   unsigned* nShocks;

   /// number of points for each shock
   std::vector<unsigned>* nShockPoints;

   /// number of shock edges for each shock
   std::vector<unsigned>* nShockEdges;

   /// number of mesh points
   std::vector<unsigned>* npoin;

   /// number of mesh elements
   std::vector<unsigned>* nelem;

   /// mesh points state (assigbale to MeshData)
   std::vector<double>* zroeVect;

   /// mesh points coordinates(assignable to MeshData)
   std::vector<double>* coorVect;

   /// vector characterizing nodes elements (assignable to MeshData)
   std::vector<int>* celnodVect;

   /// mesh points state (in array storing)
   Array2D<double>* zroe;

   /// mesh points coordinates (in array storing)
   Array2D<double>* XYZ;

   /// celnod(0)(i-elem) 1° node of i-element
  /// celnod(1)(i-elem) 2° node of i-element
  /// celnod(2)(i-elem) 3° node of i-element
  /// celnod(3)(i-elem) 4° node of i-element
   Array2D<int>* celnod;

   /// shock points coordinates
   Array3D<double>* XYZSh;

   /// upstream shock points state
   Array3D<double>* ZRoeShu;

   /// downstream shock points state
   Array3D<double>* ZRoeShd;

   /// old upstream status
   Array3D <double>* ZRoeShuOld;

   /// old downstream status
   Array3D <double>* ZRoeShdOld;

   /// working array storing upstream points coordinates
   Array3D<double> XYZShu;

   /// working array storing downstream points coordinates
   Array3D<double> XYZShd;

   /// file storing info
   FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_ComputeShockPointsState3D_hh
