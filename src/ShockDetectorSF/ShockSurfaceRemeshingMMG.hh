// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef ShockFitting_ShockSurfaceRemeshingMMG_hh
#define ShockFitting_ShockSurfaceRemeshingMMG_hh

//--------------------------------------------------------------------------//

#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"

#include <iostream>
#include <vector>
#include <fstream>
#include <utility> // defines std::pair



//--------------------------------------------------------------------------//



namespace ShockFitting {

//--------------------------------------------------------------------------//

class ShockSurfaceRemeshingMMG {
public:

  /// Constructor
  ShockSurfaceRemeshingMMG();

  /// Destructor
  ~ShockSurfaceRemeshingMMG();

//--------------------------------------------------------------------------//


  void Remeshing();

  void detectEdges();

  void createMMGfile();

  void readMMGfile();



//--------------------------------------------------------------------------//

private: // helper functions


 void setPhysicsData();

 //--------------------------------------------------------------------------//

private: // data

 //working variables
  double verta, vertb, vertc;

  //working variables
  std::vector <double> v1,v2,v3;

  /// number of vertices for each mesh element
  unsigned* nvt;

  /// number of degrees of freedom
  unsigned* ndof;

  /// number of shock points
  unsigned* nShocks;

  /// number of points for each shock
  std::vector<unsigned>* nShockPoints;

	/// number of shock edges for each shock
  std::vector<unsigned>* nShockFaces;

  /// shock points coordinates
  Array3D<double>* XYZSh;

	// Shock Faces
  Array3D<unsigned int>* ShFaces;

  Array3D <unsigned>* ShEdgePoints;


  /// number of shock points
  unsigned m_nbShPoints;


private:
	//
	std::string command;

  std::string command2;

  std::string command3;



};

//--------------------------------------------------------------------------//

}

//--------------------------------------------------------------------------//
//

#endif
