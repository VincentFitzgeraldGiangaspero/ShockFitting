// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_BndryNodePtr_hh
#define ShockFitting_BndryNodePtr_hh

//--------------------------------------------------------------------------//

#include <vector>
#include "Framework/Remeshing.hh"
#include "MathTools/Array2D.hh"
#include "Framework/FileLogManip.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a BndryNodePtr, whose task is to create an array
/// NODPTR which owns informations about boundary nodes.
///

class BndryNodePtr: public Remeshing {
public:

  /// Constructor
  /// @param objectName the concrete class name
  BndryNodePtr(const std::string& objectName);

  /// Destructor
  virtual ~BndryNodePtr();

  /// Set up this object before its first use
  virtual void setup();

  /// Unset up this object before its first use
  virtual void unsetup();

  /// create array NODPTR
  virtual void remesh();

private: // helper functions

  /// return class name
  std::string getClassName() const {return "BndryNodePtr";}

  /// local function for NODPTR computation
  void myroutine();

  /// count number of boundary points
  void getnbBndryPoints();
 
  /// local function computes vector NODPTR
  void setBndryNodePtr();

  /// assign array starting pointers
  void setAddress();

  /// assign values used in BndryNodePtr to MeshData pattern
  void setMeshData();

  /// de-allocate dynamic arrays
  void freeArray();

private: //data

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of boundary points
  std::vector<unsigned>* nbpoin;

  /// number of boundary faces in the mesh
  std::vector<unsigned>* nbfac;

  /// code charactering mesh points
  std::vector <int>* nodcod;

  /// vector characterizing boundary faces (assignable to MeshData)
  std::vector<int>* bndfacVect;

  /// vector characterizing boundary nodes
  std::vector<int>* nodptrVect;

  /// nodptr(i-poin)(0) global code number				// io non ho due edges, ho un numero di face
  /// nodptr(i-poin)(1) one of the edges it belongs to
  /// nodptr(i-poin)(2) other edge
  Array2D <int>* nodptr;

  /// bndfac(0)(i-face) 1° endpoint of i-boundary face
  /// bndfac(1)(i-face) 2° endpoint of i-boundary face
  /// bndfac(2)(i-face) boundary marker of i-boundary face
  Array2D <int>* bndfac;

  /// local vector
  /// It is used to store NBPOINTS code number
  std::vector <int> iwork_;

  /// local vector
  /// Its is used to store NBPOINTS code number indeces
  std::vector <int> iworkRank_;

  /// store log infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_BndryNodePtr_hh
