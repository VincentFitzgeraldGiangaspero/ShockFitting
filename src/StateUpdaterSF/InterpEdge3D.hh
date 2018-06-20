// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_InterpEdge3D_hh
#define ShockFitting_InterpEdge3D_hh

//--------------------------------------------------------------------------//

#include <cmath>
#include "Framework/StateUpdater.hh"
#include "Framework/FileLogManip.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

/// This class defines a InterpEdge3D, whose task is to update values in the phantom
/// nodes of the background mesh (referred to index (0)) using values in the
/// shocked mesh (referred to index (1))

class InterpEdge3D {
public:

  /// Constructor
  /// @param objectName the concrete class name
  InterpEdge3D();

  /// Destructor
  ~InterpEdge3D();

/// update phantom nodes values
  void updateEdge();

  Array2D<unsigned> getVertFaceInterp() { return vertfaceinterp;}

  void interpEdge();

private: // helper functions

  /// return class name
  std::string getClassName() const {return "InterpEdge3D";}

  /// assign variables used in InterpEdge3D to MeshData pattern
  void setMeshData();

  /// assign variables used in InterpEdge3D to PhysicsData pattern
  void setPhysicsData();

  /// set starting pointers of array 2D and 3D
  void setAddress();

  /// de-allocate dynamic arrays
  void freeArray();








private: // data

  std::vector <double> a,b,c,d,n,a_b,a_c,d_a;

  /// number of vertices
  unsigned* nvt;

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of shocks
  unsigned* nShocks;

  /// number of shock points
  std::vector<unsigned>* nShockPoints;

  /// number of mesh points
  std::vector<unsigned>* npoin;

  /// number of mesh elements
  std::vector<unsigned>* nelem;

  /// number of shock Faces
  std::vector<unsigned>* nShockFaces;


  /// mesh points state (assignable to MeshData)
  std::vector<double>* zroeVect;

  /// mesh points coordinates (assignable to MeshData)
  std::vector<double>* coorVect;

  /// vector characterizing nodes elements
  std::vector<int>* celnodVect;

  /// code characterizing mesh points
  std::vector <int>* nodcod;

  /// map vector
  std::vector<int>* M02M1;

  /// map vector -> M02M1+NPOIN(0)
  std::vector<int>* M02M12;



  /// mesh points state referring to the background mesh
  Array2D <double>* zBkg;

  /// mesh points coordinates referring the shocked mesh
  Array2D <double>* XYZ;

  /// mesh points state referring to the shocked mesh
  Array2D <double>* zroe;

  /// celnod(0)(i-elem) 1° node of i-element
  /// celnod(1)(i-elem) 2° node of i-element
  /// celnod(2)(i-elem) 3° node of i-element
  Array2D <int>* celnod;
  
  Array2D <unsigned> vertfaceinterp;


  /// shock points coordinates
  Array3D<double>* XYZSh;


  Array3D <unsigned>* ShEdgePoints;

  Array3D <unsigned>* ShFaces;

  /// store log file infos
  FileLogManip logfile;
};

//--------------------------------------------------------------------------//

} // namespace ShockFitting

//--------------------------------------------------------------------------//

#endif // ShockFitting_InterpEdge3D_hh
