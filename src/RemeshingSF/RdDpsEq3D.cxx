// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cmath>
#include "RemeshingSF/RdDpsEq3D.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/Remeshing.hh"
#include "Framework/Log.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/MeshData.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<RdDpsEq3D, Remeshing> rdDpsEq3DProv("RdDpsEq3D");

//--------------------------------------------------------------------------//

RdDpsEq3D::RdDpsEq3D(const std::string& objectName) :
  Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

RdDpsEq3D::~RdDpsEq3D()
{
}

//--------------------------------------------------------------------------//

void RdDpsEq3D::setup()
{
  LogToScreen(VERBOSE, "RdDpsEq3D::setup() => start\n");

  LogToScreen(VERBOSE, "RdDpsEq3D::setup() => end\n");
}

//--------------------------------------------------------------------------//

void RdDpsEq3D::unsetup()
{
  LogToScreen(VERBOSE, "RdDpsEq3D::unsetup()\n");
}

//--------------------------------------------------------------------------//

void RdDpsEq3D::remesh()
{
  LogToScreen(INFO, "RdDpsEq3D::remesh()\n");

  setPhysicsData();
  setMeshData();

  logfile.Open(getClassName());
  // assign start pointers of Array2D and 3D
  setAddress();

  // resize vectors and arrays
  setSize();

  // de-allocate dynamic arrays
  freeArray();

  logfile.Close();
}

//--------------------------------------------------------------------------//



//--------------------------------------------------------------------------//

void RdDpsEq3D::setAddress()
{
  unsigned start;
  start = npoin->at(0)*PhysicsInfo::getnbDofMax();
  ZRoeShu = new Array3D <double> (PhysicsInfo::getnbDofMax(),
                                  PhysicsInfo::getnbShPointsMax(),
                                  PhysicsInfo::getnbShMax(),
                                  &zroeVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDofMax() +
          PhysicsInfo::getnbShPointsMax() *
          PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDofMax();
  ZRoeShd = new Array3D <double> (PhysicsInfo::getnbDofMax(),
                                  PhysicsInfo::getnbShPointsMax(),
                                  PhysicsInfo::getnbShMax(),
                                  &zroeVect->at(start));
}

//--------------------------------------------------------------------------//

void RdDpsEq3D::setSize()
{
  XYZSh_New.resize(PhysicsInfo::getnbDim(),
                  PhysicsInfo::getnbShPointsMax() );
  ZRoeShu_New.resize(PhysicsInfo::getnbDofMax(),
                     PhysicsInfo::getnbShPointsMax() );
  ZRoeShd_New.resize(PhysicsInfo::getnbDofMax(),
                     PhysicsInfo::getnbShPointsMax() );

}

//--------------------------------------------------------------------------//

void RdDpsEq3D::freeArray()
{
  delete ZRoeShu; delete ZRoeShd;
}

//--------------------------------------------------------------------------//

void RdDpsEq3D::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints = PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  nShockEdges = PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  XYZSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
}

//--------------------------------------------------------------------------//

void RdDpsEq3D::setMeshData()
{
  zroeVect = MeshData::getInstance().getData <vector <double> >("ZROE");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
