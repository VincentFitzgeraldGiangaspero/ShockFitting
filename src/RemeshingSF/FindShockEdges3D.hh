// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_FindShockEdges3D_hh
#define ShockFitting_FindShockEdges3D_hh

//--------------------------------------------------------------------------//

#include "Framework/Remeshing.hh"
#include "Framework/FileLogManip.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

// vincent 04/05/2017

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a FindShockEdges3D, whose task is to find the Shock Points
/// on the edges of the surface

class FindShockEdges3D {
public:

  /// Constructor
  /// @param objectName the concrete class name
  FindShockEdges3D(const std::string& objectName);

  /// Destructor
  virtual ~FindShockEdges3D();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// re-distribute shock points
  virtual void remesh();

private: // helper functions

  /// return class name
  std::string getClassName() const {return "FindShockEdges3D";}

  ///assign PhysicsData values to FindShockEdges3D
  void setPhysicsData();

  ///assign MeshData values to FindShockEdges3D
  void setMeshData();

  /// assign start pointers of Array2D and 3D
  void setAddress();

  /// resize vectors and arrays
  void setSize();

  /// de-allocate dynamic arrays
  void freeArray();

private: // data

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

  /// upstream state
  Array3D <double>* ZRoeShu;

  /// downstream state
  Array3D <double>* ZRoeShd;

  /// mesh points status
  std::vector <double>* zroeVect;

  /// store file log infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_FindShockEdges3D_hh
