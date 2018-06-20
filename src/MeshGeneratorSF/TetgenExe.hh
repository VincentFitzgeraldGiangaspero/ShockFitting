// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_TetgenExe_hh
#define ShockFitting_TetgenExe_hh

//--------------------------------------------------------------------------//

#include <vector>
#include <string>
#include "Framework/MeshGenerator.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a TetgenExe, whose task is to call Tetgen
/// Mesh Generator.

class TetgenExe : public MeshGenerator {
public:

  /// Constructor
  TetgenExe(const std::string& objectName);

  /// Destructor
  virtual ~TetgenExe();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its last use
  virtual void unsetup();

  /// generate the new mesh
  virtual void generate();

  /// generate the new mesh from a given processing file
  virtual void generate(std::string);

private:

  /// dummy string
  std::string command;

  /// name of current file
  std::stringstream* fname;

};

//--------------------------------------------------------------------------//

} // namespace ShockFitting_TetgenExe_hh

//--------------------------------------------------------------------------//

#endif // ShockFitting_TetgenExe_hh
