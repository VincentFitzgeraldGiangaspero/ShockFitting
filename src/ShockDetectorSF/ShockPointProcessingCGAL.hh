// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ShockPointProcessingCGAL_hh
#define ShockFitting_ShockPointProcessingCGAL_hh


//--------------------------------------------------------------------------//
#include "ShockDetectorSF/ShockPointProcessingCGAL.hh"
#include "MathTools/MinMax.hh"
#include "Framework/ChemicalInfo.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "SConfig/ObjectProvider.hh"

#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

#include <CGAL/Kernel_traits.h>
#include <CGAL/Origin.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/edge_aware_upsample_point_set.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/bilateral_smooth_point_set.h>
#include <CGAL/tags.h>


// from surface reconstr
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/tuple.h>
#include <boost/lexical_cast.hpp>
#include <algorithm>

// vector from wlop
#include <vector>

#include <utility> // defines std::pair
#include <list>
#include <fstream>
#include <iostream>


// Types Definition
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel_pred; //from normals
typedef CGAL::Simple_cartesian<double> Kernel; //from jet and wlop
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::cpp11::array<std::size_t,3> Facet;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;


// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

//--------------------------------------------------------------------------//


//--------------------------------------------------------------------------//

struct Perimeter {
	double bound;
	Perimeter(double bound)
	: bound(bound)
	{}
	template <typename AdvancingFront, typename Cell_handle>
	double operator() (const AdvancingFront& adv, Cell_handle& c,
					 const int& index) const
	{
	// bound == 0 is better than bound < infinity
	// as it avoids the distance computations
	if(bound == 0){
	  return adv.smallest_radius_delaunay_sphere (c, index);
	}
	// If perimeter > bound, return infinity so that facet is not used
	double d  = 0;
	d = sqrt(squared_distance(c->vertex((index+1)%4)->point(),
							  c->vertex((index+2)%4)->point()));
	if(d>bound) return adv.infinity();
	d += sqrt(squared_distance(c->vertex((index+2)%4)->point(),
							   c->vertex((index+3)%4)->point()));
	if(d>bound) return adv.infinity();
	d += sqrt(squared_distance(c->vertex((index+1)%4)->point(),
							   c->vertex((index+3)%4)->point()));
	if(d>bound) return adv.infinity();
	// Otherwise, return usual priority value: smallest radius of
	// delaunay sphere
	return adv.smallest_radius_delaunay_sphere (c, index);
	}
	};
//--------------------------------------------------------------------------//


//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

class ShockPointProcessingCGAL {
public:

  /// Constructor
  ShockPointProcessingCGAL();

  /// Destructor
  ~ShockPointProcessingCGAL();

//--------------------------------------------------------------------------//

  virtual void PointReconstruction();

  void Working();

  void Wlop();

  void Normals();

  void Smoothing();

  void AutoReconstruction();

  // void SurfaceMeshingOLD();


//--------------------------------------------------------------------------//

private: // helper functions


 /// assign strating pointers to array
 void setAddress();

 /// assign variables used in ComputeShockPointsState to MeshData
 void setMeshData();

 /// assign variables used in ComputeShockPointsState to PhysicsData
 void setPhysicsData();

 /// de-allocate dynamic array
 void freeArray();



 //--------------------------------------------------------------------------//

private: // data

  std::vector<Point> points_shock;


  std::vector<Point> m_points_shock;


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

	/// number of shock edges for each shock
  std::vector<unsigned>* nShockFaces;

  // Shock Faces
  Array3D<unsigned int>* ShFaces;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of mesh elements
  std::vector<unsigned>* nelem;

  /// mesh points coordinates(assignable to MeshData)
  std::vector<double>* coorVect;

  /// vector characterizing nodes elements (assignable to MeshData)
  std::vector<int>* celnodVect;

  /// mesh points coordinates (in array storing)
  Array2D<double>* XYZ;

  /// celnod(0)(i-elem) 1° node of i-element
  /// celnod(1)(i-elem) 2° node of i-element
  /// celnod(2)(i-elem) 3° node of i-element
  /// celnod(3)(i-elem) 4° node of i-element
  Array2D<int>* celnod;

  /// shock points coordinates
  Array3D<double>* XYZSh;

  /// number of shock points
  unsigned m_nbShPoints;


public:
	
	double minDomainXCoordinate;

	double maxDomainXCoordinate;

	double minDomainYCoordinate;

	double maxDomainYCoordinate;

	double minDomainZCoordinate;

	double maxDomainZCoordinate;

  /// file storing info
  FileLogManip logfile;

};

//--------------------------------------------------------------------------//

}

//--------------------------------------------------------------------------//
//

#endif
