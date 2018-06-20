// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShockDetectorSF/ShockSurfaceMeshingCGAL.hh"
#include "Framework/ChemicalInfo.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"



#include <iostream>
#include <fstream>
#include <algorithm>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/tuple.h>
#include <boost/lexical_cast.hpp>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Origin.h>

//--------------------------------------------------------------------------//


typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::cpp11::array<std::size_t,3> Facet;


//--------------------------------------------------------------------------//

namespace std {
  std::ostream&
  operator<<(std::ostream& os, const Facet& f)
  {
    os << "3 " << f[0] << " " << f[1] << " " << f[2];
    return os;
  }
}


//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//


namespace ShockFitting {

//--------------------------------------------------------------------------//

/// Constructor
ShockSurfaceMeshingCGAL::ShockSurfaceMeshingCGAL(){
 }

//--------------------------------------------------------------------------//

/// Destructor
ShockSurfaceMeshingCGAL::~ShockSurfaceMeshingCGAL(){
 }

 //--------------------------------------------------------------------------//

void ShockSurfaceMeshingCGAL::SurfaceMeshing(){

  // std::cout << "     => ShockSurfaceMeshingCGAL::SurfaceMeshing()";
  //
  // setMeshData();
  // setPhysicsData();
  // setAddress();
  //
  // // unsigned ISH=0;
  // //
  // // for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
  // //  points_shock.push_back(Point((*XYZSh)(0,ISHPOIN,ISH),(*XYZSh)(1,ISHPOIN,ISH),(*XYZSh)(2,ISHPOIN,ISH)));
  // // }
  // //
  // // cout << "points_shock.size() in SurfaceMeshing = " << points_shock.size() << " \n";
  // //
  // //
  // // double per = 0;
  // //
  // // std::vector<Facet> facets;
  // //
  // //
  // // Perimeter perimeter(per);
  // // CGAL::advancing_front_surface_reconstruction(points_shock.begin(),
  // // 									   points_shock.end(),
  // // 									   std::back_inserter(facets),
  // // 									   perimeter);
  // //
  // //
  // // nShockFaces->at(0) = facets.size();
  // //
  // // ShFaces->resize(3,
  // //               facets.size(),
  // //               PhysicsInfo::getnbShMax());
  // //
  // // cout << "nShockFaces in SurfaceMeshing = " << nShockFaces->at(0) << " \n";
  // //
  // // // .off file of surface meshed
  // // std::ofstream out_surf("shock.off");
  // // out_surf << "OFF\n" << points_shock.size() << " " << facets.size() << " 0\n";
  // // std::copy(points_shock.begin(),
  // // 	points_shock.end(),
  // // 	std::ostream_iterator<Point>(out_surf, "\n"));
  // // std::copy(facets.begin(),
  // // 	facets.end(),
  // // 	std::ostream_iterator<Facet>(out_surf, "\n"));
  // //
  // //
  // //
  // //
  // // // .node file
  // // std::ofstream out_node("log/shock_surface.node");
  // //
  // // out_node << points_shock.size() << "\n";
  // //
  // // for (std::size_t i=0; i<points_shock.size(); ++i){
  // //
  // // out_node << i+33898+1 << " " << points_shock[i] << " 7420.3797440538401133  2586.9173518748098104 " <<
  // //
  // // out_node << "509.3701808266359876   31.6751022599442997  2603.7394822446799481 0\n";
  // //
  // // }
  // // //7420.3797440538401133  2586.9173518748098104  509.3701808266359876   31.6751022599442997  2603.7394822446799481  0
  // //
  // // std::ofstream out_poly("log/shock_surface.poly");
  // //
  // // out_poly << facets.size() << "\n";
  // //
  // // for( unsigned IFACE=0; IFACE<nShockFaces->at(0); IFACE++){
  // //   (*ShFaces)(0,IFACE,ISH)=facets[IFACE][0];
  // //   (*ShFaces)(1,IFACE,ISH)=facets[IFACE][1];
  // //   (*ShFaces)(2,IFACE,ISH)=facets[IFACE][2];
  // // }
  // //
  // // // std::cout << " (*ShFaces)(0,0,0) " << (*ShFaces)(0,0,0) << " \n";
  // // // std::cout << " (*ShFaces)(1,0,0) " << (*ShFaces)(1,0,0) << " \n";
  // // // std::cout << " (*ShFaces)(2,0,0) " << (*ShFaces)(2,0,0) << " \n";
  // // //
  // // // std::cout << " (*XYZSh)(0,408,0) " << (*XYZSh)(0,(*ShFaces)(0,0,0),0) << " \n";
  // // // std::cout << " (*XYZSh)(1,408,0) " << (*XYZSh)(1,(*ShFaces)(0,0,0),0) << " \n";
  // // // std::cout << " (*XYZSh)(2,408,0) " << (*XYZSh)(2,(*ShFaces)(0,0,0),0) << " \n";
  // //
  // // for (std::size_t j=0; j<facets.size(); j++){
  // // out_poly << "1     0   5 \n";
  // // out_poly  << "3 " << facets[j][0]+33899 << " " <<facets[j][1]+33899 << " " << facets[j][2]+33899  << "\n";
  // // }
  //
  //
  // freeArray();
}

//--------------------------------------------------------------------------//

void ShockSurfaceMeshingCGAL::setAddress()
{
  unsigned totsize = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                                 PhysicsInfo::getnbShPointsMax();
  XYZ = new Array2D<double> (PhysicsInfo::getnbDim(),
                            totsize,
                            &coorVect->at(0));



}

//--------------------------------------------------------------------------//

void ShockSurfaceMeshingCGAL::freeArray()
{
  delete XYZ;
}

//--------------------------------------------------------------------------//

void ShockSurfaceMeshingCGAL::setMeshData()
{
  nvt = MeshData::getInstance().getData<unsigned>("NVT");
  nelem = MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  // zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  // celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
}

//--------------------------------------------------------------------------//

void ShockSurfaceMeshingCGAL::setPhysicsData()
{
  // nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
   PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  // nShockEdges =
  //  PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");		// decidere che fare con questo
  nShockFaces =
   PhysicsData::getInstance().getData <vector <unsigned> > ("nShockFaces");		// aggiunto

  XYZSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYZSH");
  ShFaces =
   PhysicsData::getInstance().getData <Array3D<unsigned> > ("ShFaces");

}

//--------------------------------------------------------------------------//


}//namespace ShockFitting
