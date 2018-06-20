// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ComputeStateDps4Pg3DAssignUpstream_hh
#define ShockFitting_ComputeStateDps4Pg3DAssignUpstream_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include "Framework/FileLogManip.hh"
#include "Framework/StateUpdater.hh"
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

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel_2;
typedef Kernel_2::Point_3 Point_3;
typedef CGAL::Tetrahedron_3<Kernel_2> Tetrahedron_3;



//definition of a non-mutable lvalue property map,
//with the get function as a friend function to give it
//access to the private member


typedef CGAL::Search_traits_3<Kernel_2>                                        Traits_base;
typedef CGAL::Search_traits_adapter<std::size_t,My_point_property_map,Traits_base> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits>                      K_Neighbor_search;
typedef K_Neighbor_search::Tree                                         Tree;
typedef Tree::Splitter                                                  Splitter;
typedef K_Neighbor_search::Distance                                     Distance;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a ComputeStateDps4Pg3DAssignUpstream, whose task is to update
/// the solution within the grid-points located on the discontinuities
/// by enforcing Rankine-Hugoniot relations between the upstream and
/// the downstream states for a perfect gas model

class ComputeStateDps4Pg3DAssignUpstream : public StateUpdater {
public:

  /// Constructor
  /// @param objectName the concrete class name
  ComputeStateDps4Pg3DAssignUpstream(const std::string& objectName);

  /// Destructor
  ~ComputeStateDps4Pg3DAssignUpstream();

  /// Set up this object before its first use
  void setup();

  /// Unset up this object before its first use
  void unsetup();

  /// Update solution
  void update();

private: // functions

  /// assign start pointers of Array2D and 3D
  void setAddress();

  /// resize discontinuity speed array
  void setDiscSpeedSize();

  /// assign values used in ComputeState to MeshData pattern
  void setMeshData();

  /// assign values used in ComputeState to PhysicsData pattern
  void setPhysicsData();

  /// de-allocate the dynamic arrays
  void freeArray();

private: // helper functions

  /// return class name
  std::string getClassName () const {return std::string("ComputeStateDpsPg");}

  /// upload downstream status
  void recoverDownState(unsigned, unsigned);

  // upload old downstream status
  void recoverDownOlDState(unsigned, unsigned);

  /// upload upstream status
  void recoverUpState(unsigned, unsigned);

  /// save old downstream status
  void saveDownState(unsigned, unsigned);

  /// compute new upstream or downstream status
  /// and store solution in ZRoe arrays
  void computeDownState(unsigned, unsigned);

  void computeUpState(unsigned, unsigned);

  void InterpShockVelocity(unsigned, unsigned);

  void InterpDownState(unsigned);

  void InterpUpState(unsigned);

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of shocks
  unsigned* nShocks;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of shock points for each shock
  std::vector<unsigned>* nShockPoints;

  /// number of shock faces for each shock
  std::vector<unsigned>* nShockFaces;

  /// type of shock
  std::vector <std::string>* typeSh;

  /// mesh points status
  std::vector<double>* zroeVect;

  /// shock points normal vectors
  Array3D <double>* vShNor;

  /// shock/discontinuity speed
  Array3D <double>* WSh;

  /// upstream status
  Array3D <double>* ZRoeShu;

  /// downstream status
  Array3D <double>* ZRoeShd;

  /// shock points coordinates
  Array3D <double>* XYZSh;

  /// old upstream status
  Array3D <double>* ZRoeShuOld;

  /// old downstream status
  Array3D <double>* ZRoeShdOld;

  /// total number of shock points
  unsigned TotnbShockPoints;

  /// dummy variables storing normal vector values
  double dx; double dy; double dz;

  /// dummy variables storing tangential vector values
  double dxt; double dyt; double dzt;

  /// discontinuity speed
  double WS;

  /// working variable storing the riemann invariant
  double R2;

  /// helping variables
  double help;

  /// work vector  used to store upstream values
  std::vector<double> xu;

  /// working vector used to store downstream values
  std::vector<double> xd;

  /// store log file infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_ComputeStateDps4Pg3DAssignUpstream_hh
