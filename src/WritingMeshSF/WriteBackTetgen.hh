// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_WriteBackTetgen_hh
#define ShockFitting_WriteBackTetgen_hh

//--------------------------------------------------------------------------//

#include <stdio.h>
#include <vector>
#include "Framework/WritingMesh.hh"
#include "MathTools/Array2D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a WriteBackTriangle, whose task is to write the
/// updates value of ZROE on triangle format file .node
/// File node format:
/// First line: <# of vertices> <dimension (must be 2)> <# of
/// attributes> <# of boundary markers (0 or 1)>
/// Remaining lines: <vertex #> <x> <y> [attributes]
/// [boundary marker] 

class WriteBackTetgen : public WritingMesh {
public:

  /// Constructor
  /// @param objectName the concrete class name
  WriteBackTetgen(const std::string& objectName);

  /// Destructor
  virtual ~WriteBackTetgen();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its last use
  virtual void unsetup();

  /// Write Tetgen fomat file
  virtual void write();

  /// write one or more output files
  virtual void writeShockTecplot(unsigned);

private: // helper functions

  /// assign variables used in WriteBackTetgen to MeshData pattern
  void setMeshData();

  /// assign variables used in WriteBackTetgen to PhysicsData pattern
  void setPhysicsData();

  /// assign starting pointers for Array2D
  void setAddress();

  /// de-allocate dynamic arrays
  void freeArray();

private: // data

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of shock points
  std::vector<unsigned>* npoin;

  /// number of shock points
    std::vector<unsigned>* nelem;

  /// mesh points state
  std::vector<double>* zroeVect;

  /// mesh points coordinates
  std::vector<double>* coorVect;

  /// mesh points coordinates
  std::vector<int>* celnodVect;


  /// code characterizing mesh points
  std::vector<int>* nodcod;

  /// name of the current output file
  std::string* fnameBack;

  /// mesh points state (in array storing)
  Array2D<double>* Zroe;

  /// mesh points coordinates (in array storing)
  Array2D<double>* XYZ;

    Array2D<int>* celnod;

  /// variable writing on node file
  FILE* file;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_WriteBackTetgen
