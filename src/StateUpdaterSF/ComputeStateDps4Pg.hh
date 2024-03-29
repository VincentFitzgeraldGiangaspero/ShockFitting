// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
// 
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ComputeStateDps4Pg_hh
#define ShockFitting_ComputeStateDps4Pg_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include "Framework/FileLogManip.hh"
#include "Framework/StateUpdater.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a ComputeStateDps4Pg, whose task is to update 
/// the solution within the grid-points located on the discontinuities
/// by enforcing Rankine-Hugoniot relations between the upstream and 
/// the downstream states for a perfect gas model

class ComputeStateDps4Pg : public StateUpdater {
public:

  /// Constructor
  /// @param objectName the concrete class name
  ComputeStateDps4Pg(const std::string& objectName);

  /// Destructor
  ~ComputeStateDps4Pg();

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

  /// upload upstream status
  void recoverUpState(unsigned, unsigned);

  /// save old downstream status
  void saveDownState(unsigned, unsigned);

  /// compute new upstream or downstream status
  /// and store solution in Zroe arrays
  void computeDownState(unsigned, unsigned);

  void computeUpState(unsigned, unsigned);

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of shocks
  unsigned* nShocks;

  /// number of mesh points
  std::vector<unsigned>* npoin;
  
  /// number of shock points for each shock
  std::vector<unsigned>* nShockPoints;

  /// number of shock edges for each shock
  std::vector<unsigned>* nShockEdges;

  /// type of shock
  std::vector <std::string>* typeSh;

  /// mesh points status
  std::vector<double>* zroeVect;

  /// shock points normal vectors
  Array3D <double>* vShNor;

  /// shock/discontinuity speed
  Array3D <double>* WSh;

  /// upstream status
  Array3D <double>* ZroeShu;

  /// downstream status
  Array3D <double>* ZroeShd;

  /// old upstream status
  Array3D <double>* ZroeShuOld;

  /// old downstream status
  Array3D <double>* ZroeShdOld;

  /// total number of shock points
  unsigned TotnbShockPoints;

  /// dummy variables storing normal vector values
  double dx; double dy;

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

#endif // ShockFitting_ComputeStateDps4Pg_hh
