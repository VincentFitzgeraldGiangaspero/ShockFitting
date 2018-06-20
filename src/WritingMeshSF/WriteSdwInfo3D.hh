// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_WriteSdwInfo3D_hh
#define ShockFitting_WriteSdwInfo3D_hh

// vincent 29/04/2017

//--------------------------------------------------------------------------//

#include <vector>
#include "Framework/FileLogManip.hh"
#include "Framework/WritingMesh.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a WriteSdwInfo3D, whose task is to write file
/// sh99.dat containg informations about shocks and or discontinuities

class WriteSdwInfo3D : public WritingMesh {
public:

  /// Constructor
  /// @param objectName the concrete class name
  WriteSdwInfo3D(const std::string& objectName);

  /// Destructor
  virtual ~WriteSdwInfo3D();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its last use
  virtual void unsetup();

  /// write the sh99 file
  virtual void write();

 /// write the sh99 file
  virtual void writeShockTecplot(unsigned);

private: // helper functions

  /// return class name
  std::string getClassName() const { return std::string("WriteSdwInfo3D"); }

  /// assign variables used in WriteSdwInfo3D to MeshData pattern
  void setMeshData();

  /// assign variables used in WriteSdwInfo3D to PhysicsData pattern
  void setPhysicsData();

  /// assign starting pointers to Array3D
  void setAddress();

  /// write SHinSPPs elements of sh99 file
  void writeSHinSPPs(unsigned, unsigned);

  /// de-allocate dynamic arrays
  void freeArray();

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// number oer shocks
  unsigned* nShocks;

  /// number of special points
  unsigned* nSpecPoints;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of shock points for each shock
  std::vector<unsigned>* nShockPoints;

  /// number of shock points
  std::vector<unsigned>* nShockFaces;

  /// type of shock
  std::vector<std::string>* typeSh;

  /// type of special points
  std::vector <std::string>* typeSpecPoints;

  /// mesh points state
  std::vector<double>* zroe;

  /// shock points upstream state
  Array3D <double>* ZroeShu;

  /// shock points downstream state
  Array3D <double>* ZroeShd;

  /// shock points coordinates
  Array3D <double>* XYZSh;

  /// shock points coordinates
  Array3D <double>* WSh;

  /// shock points coordinates
  Array3D <unsigned int>* ShFaces;

  /// array characterizing special points
  Array3D <unsigned>* SHinSPPs;

  /// variable writing on sh99 file
  FILE* file;

  FILE* off_file;

  /// variable writing on upstream file
  FILE* shplt_up;

  /// variable writing on upstream file
  FILE* shplt_down;

  /// store log file infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_WriteSdwInfo3D.hh
