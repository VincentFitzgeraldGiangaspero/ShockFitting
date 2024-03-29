// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_CoNorm4Pg3D_hh
#define ShockFitting_CoNorm4Pg3D_hh

// vincent 9/3/17

//--------------------------------------------------------------------------//

#include <fstream>
#include "Framework/Remeshing.hh"
#include "Framework/FileLogManip.hh"
#include "MathTools/Array3D.hh"


//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a CoNorm4Pg3D, whose task is to compute the normal unit
/// vectors to the shocks and discontinuities in shock/discontinuity points
/// for NDOF = 5 && MODEL = PG

class CoNorm4Pg3D : public Remeshing {
public:

  /// Constructor
  /// @param objectName the concrete class name
  CoNorm4Pg3D(const std::string& objectName);

  /// Destructor
  ~CoNorm4Pg3D();

  /// Set up this object before its first use
  void setup();

  /// Unset up this object before its first use
  void unsetup();

  /// Compute normal vectors
  void remesh();

private: // functions

  /// assign start pointers of Array2D and 3D
  void setAddress();

  /// resize vectors and arrays
  void setSize();

  /// assign values used in CoNorm to MeshData pattern
  void setMeshData();

  /// assign values used in CoNorm to PhysicsData pattern
  void setPhysicsData();

  /// de-allocate dynamic arrays
  void freeArray();

private: // helper functions

  /// return class name
  std::string getClassName () const {return std::string("CoNorm4Pg3D");}

  /// compute tangetial vectors
  void computeTau(unsigned, unsigned);

  void compute3Dnormals();

  void compute3Dnormals_dependency();

  /// recover shock point status
  void recoverState(std::string, unsigned, unsigned);

  /// fix normal vector for typeShock = S
  void setVShNorForStype();

  /// fix normal vector for WPNRX special point
  void setVShNorForWPNRX(unsigned);

  /// fix normal vector for C special point
  void setVShNorForC(unsigned);

  /// fix normal vector for TP special point
  void setVShNorForTP(unsigned);

  /// write tecplot file
  void writeTecPlotFile();

  /// set indeces characterizing shock
  void setShockIndeces(unsigned, unsigned);

  /// take one point forward coordinates
  void onePointForward(unsigned, unsigned);

  /// take two points forward coordinates
  void twoPointsForward(unsigned, unsigned);

  /// take one point backward coordinates
  void onePointBackward(unsigned, unsigned);

  /// take two points backward coordinates
  void twoPointsBackward(unsigned, unsigned);

  /// set variables used to compute normal vectors
  void setLm();
  void setLp();
  void setTauIp1ToZero();
  void setTauIp2ToZero();
  void setTauIm1ToZero();
  void setTauIm2ToZero();

private: // data

  /// dummy variables for the shock speed
  double ush,vsh;
  double ui, vi;

  /// dummy variables for the shock points coordinates
  double xi,yi,xj,yj,xj2,yj2;

  /// dummy variables for the normal vectors computation
  double taux, tauy, tau;
  double tauxip1, tauyip1, tauxip2, tauyip2;
  double tauxim1, tauyim1, tauxim2, tauyim2;
  double lp1, lp2, lm1, lm2, lp12, lm12, lm22, lp22;
  double dum, nx1, nx2, ny1, ny2, nx4, ny4;

  /// dummy variables for the status recovering
  double uj, vj, roj, help, aj, pj;

  /// dummy variables for the shock dependencies
  int depip1, depim1;

  /// dummy variables for shock indeces
  unsigned I;
  std::vector<unsigned> ISH;
  std::vector<unsigned> IP;

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of shocks
  unsigned* nShocks;

  /// number of special points
  unsigned* nSpecPoints;

  /// number of shock points for each shock
  std::vector<unsigned>* nShockPoints;

  /// vector storing components of the normal unit vector
  Array3D<double> shock_normals;

  // Shock Faces connectivity
  Array3D<unsigned int>* ShFaces;

  /// number of shock edges for each shock
  std::vector <unsigned>* nShockFaces;

  /// type of shock
  std::vector<std::string>* typeSh;

  /// type of special points
  std::vector<std::string>* typeSpecPoints;

  /// mesh points status
  std::vector <double>* zroe;

  /// store log file infos
  FileLogManip logfile;

  /// shock points coordinates
  Array3D <double>* XYZSh;

  /// downstream status
  Array3D <double>* ZRoeShd;

  /// shock points normal vectors
  Array3D <double>* vShNor;

  /// array characterizing special points
  Array3D <unsigned>* SHinSPPs;
};

//----------------------------------------------------------------------------//

} // namespace ShockFitting

#endif // ShockFitting_CoNorm4Pg3D
