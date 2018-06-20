// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_FindShockSurfaceNormals3DCGAL_hh
#define ShockFitting_FindShockSurfaceNormals3DCGAL_hh

//--------------------------------------------------------------------------//

#include "Framework/FileLogManip.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a FindShockSurfaceNormals3DCGAL, whose task is to
/// compute the normal unit vector to the shock surface

/// @author VIncent

class FindShockSurfaceNormals3DCGAL {
public:

  /// Constructor
  FindShockSurfaceNormals3DCGAL();

  /// Destructor
  ~FindShockSurfaceNormals3DCGAL();

  /// return shock points normal vector
  Array3D<double> getUnitNormalsVector3D() { return normals; }

  /// compute shock points normals
  void computeUnitNormals3D();

private: // helper functions

  /// assign variables used in FindShockSurfaceNormals3DCGAL to PhysicsData
  void setPhysicsData();

  /// compute shock points normals in the special points
  void computeUnitNormalsSpecialPoints3D();

private: // data

  /// vector storing components of the normal unit vector
  Array3D<double> normals;

  /// number of shocks
  unsigned* nShocks;

  /// number of special points
  unsigned* nSpecPoints;

  /// number of points for each shock
  std::vector<unsigned>* nShockPoints;

  /// type of special points
  std::vector <std::string>* typeSpecPoints;

  /// shock points coordinates
  Array3D <double>* XYZSh;

  /// array characterizing special points
  Array3D <unsigned>* SHinSPPs;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_FindShockSurfaceNormals3DCGAL_hh
