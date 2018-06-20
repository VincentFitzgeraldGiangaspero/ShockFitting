// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_LaplacianSmoothing3D_hh
#define ShockFitting_LaplacianSmoothing3D_hh

//--------------------------------------------------------------------------//

#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a LaplacianSmoothing3D, whose task is to re-distribute shock points
/// in each shock with uniform spacing

class LaplacianSmoothing3D {
public:

  /// Constructor
  /// @param objectName the concrete class name
  LaplacianSmoothing3D();

  /// Destructor
  ~LaplacianSmoothing3D();

  //--------------------------------------------------------------------------//


  /// re-distribute shock points
  void Smooth();

  void WriteShockFile();

  /// find Shock Points on the shock surface edges
  void findShockEdges();

  /// Manual LaplacianSmooth with 3 steps of the algorithm
  void LaplacianSmooth();

private: // helper function


  ///assign PhysicsData values to LaplacianSmoothing3D
  void setPhysicsData();

  ///assign MeshData values to LaplacianSmoothing3D
  void setMeshData();

  /// assign start pointers of Array2D and 3D
  void setAddress();

  /// resize vectors and arrays
  void setSize();

  /// de-allocate dynamic arrays
  void freeArray();

private: // data

  //working variables
  double verta, vertb, vertc;

  //working variables
  std::vector <double> v1,v2,v3;

  /// dummy variable used as index
  unsigned ISH;

  /// shock edge length
  double Sh_Edge_length;

  /// number of shock points of redistribution
  unsigned nShockPoints_new;

  /// array stores new shock points
  Array2D <double> XYZSh_New;

  /// array stores new zroe upstream values
  Array2D <double> ZRoeShu_New;

  /// array stores new zroe downstream values
  Array2D <double> ZRoeShd_New;

  /// shock length edge
  std::vector <double> Sh_ABSC;

  /// new shock length edge
  std::vector <double> Sh_ABSC_New;

  /// nof degree of freedom
  unsigned* ndof;

  /// number of shocks
  unsigned* nShocks;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of shock edges for each shock
  std::vector <unsigned>* nShockEdges;

  /// number of shock points for each shock
  std::vector <unsigned>* nShockPoints;

  /// shock points coordinates
  Array3D <double>* XYZSh;

  /// shock points coordinates
  Array3D <unsigned>* ShEdgePoints;

  // Shock Faces connectivity
  Array3D<unsigned int>* ShFaces;

  /// number of shock edges for each shock
  std::vector <unsigned>* nShockFaces;

  /// upstream state
  Array3D <double>* ZRoeShu;

  /// downstream state
  Array3D <double>* ZRoeShd;

  /// mesh points status
  std::vector <double>* zroeVect;


};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_LaplacianSmoothing3D_hh
