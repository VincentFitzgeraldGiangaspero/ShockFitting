// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/MoveDps4Pg3D.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "SConfig/ObjectProvider.hh"

// vincent 19/04/2017


//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<MoveDps4Pg3D, StateUpdater> moveDpsPg3DProv("MoveDps4Pg3D");

//----------------------------------------------------------------------------//

MoveDps4Pg3D::MoveDps4Pg3D(const std::string& objectName) :
 StateUpdater(objectName)
{
}

//----------------------------------------------------------------------------//

MoveDps4Pg3D::~MoveDps4Pg3D()
{
}

//----------------------------------------------------------------------------//

void MoveDps4Pg3D::setup()
{
  LogToScreen(VERBOSE,"MoveDps4Pg3D::setup() => start\n");

  LogToScreen(VERBOSE,"MoveDps4Pg3D::setup() => end\n");
}

//----------------------------------------------------------------------------//

void MoveDps4Pg3D::unsetup()
{
  LogToScreen(VERBOSE,"MoveDps4Pg3D::unsetup()\n");
}

//----------------------------------------------------------------------------//

void MoveDps4Pg3D::update()
{
  LogToScreen(INFO,"MoveDps4Pg3D::update()\n");

  setMeshData();

  setPhysicsData();

  setAddress();

  logfile.Open(getClassName().c_str());

  WShMax = 0;
  
  double roref;
  roref=ReferenceInfo::getpref()/(ReferenceInfo::getRgas()*ReferenceInfo::getTref());

  // compute the dt max
  dt = 1e39;
  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   unsigned iShock = ISH+1;
   logfile("Shock n. ", iShock, "\n");
   for(unsigned IV=0; IV<nShockPoints->at(ISH); IV++) {
    unsigned iShockPoint = IV+1;

     ro = (*ZroeSh)(0,IV,ISH)*(*ZroeSh)(0,IV,ISH);
     help = pow((*ZroeSh)(2,IV,ISH),2) + pow((*ZroeSh)(3,IV,ISH),2)+pow((*ZroeSh)(4,IV,ISH),2);
     p = PhysicsInfo::getGm1()/PhysicsInfo::getGam() *
         ((*ZroeSh)(0,IV,ISH) * (*ZroeSh)(1,IV,ISH)-0.5*help);


    //ro = ((*ZroeSh)(0,IV,ISH)/((*ZroeSh)(4,IV,ISH)*ReferenceInfo::getRgas()))/roref;
    //p = (*ZroeSh)(0,IV,ISH)/ReferenceInfo::getpref();
    
    //ro = ((*ZroeSh)(0,IV,ISH)/((*ZroeSh)(4,IV,ISH)*ReferenceInfo::getRgas()));
    //p = (*ZroeSh)(0,IV,ISH);

    a = sqrt(PhysicsInfo::getGam()*p/ro);

    dum = 0;
    for(unsigned K=0; K<3; K++) { dum = dum + pow((*WSh)(K,IV,ISH),2); }
    WShMod = sqrt(dum);

    logfile("Shock point n. ", iShockPoint," speed: ", WShMod);
    logfile( " a_sound: ", a, "\n");

    if (WShMod>WShMax) { WShMax = WShMod; }

    dum = (MeshData::getInstance().getSHRELAX()) *
          (MeshData::getInstance().getDXCELL())  *
          (MeshData::getInstance().getSNDMIN()) /(a+WShMod);

   // dum = (MeshData::getInstance().getSHRELAX()) *
   //       (MeshData::getInstance().getSNDMIN()) /(a+WShMod);


    if(dt>dum) { dt = dum; }
   }
  }

  logfile("DT max: ", dt, "\n");

  std::cout << "-DT max: " << dt <<'\n';
  std::cout << "-WShMax: " << WShMax <<'\n';



  // if the connectivity freezed option is true and
  // if the freezed connectivity ratio is less than the value
  // chosen in the input.case, the freezed connectivity bool variable
  // is set to true
  // if(MeshData::getInstance().freezedConnectivityOption()) {
  //
  // //  // open file storing freezedConnectivityRatioValues
  //  ofstream memFCR("FreezedConnectivityRatioValues.dat",ios::app);
  //  memFCR << WShMax*dt/MeshData::getInstance().getSNDMIN() << endl;
  //  memFCR.close();
  //
  //  if(WShMax*dt/MeshData::getInstance().getSNDMIN() <
  //   MeshData::getInstance().getFreezedAdimConnectivityRatio())
  //  {
  //   MeshData::getInstance().setFreezedConnectivity(true);
  //   logfile("\n (!) Freezed connectivity\n");
  //   logfile("Freezed connectivity ratio = ",
  //           WShMax*dt/MeshData::getInstance().getSNDMIN(),"\n\n");
  //   LogToScreen(INFO,"MoveDps4Pg3D::warning => freezed connectivity\n");
  //  }
  // }



  // move the discontinuity
  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   unsigned iShock = ISH+1;
   logfile("Shock/Discontinuity n. ", iShock, "\n");
   for(unsigned IV=0; IV<nShockPoints->at(ISH); IV++) {
    for(unsigned I=0; I<3; I++) {
     (*XYZSh)(I,IV,ISH) = (*XYZSh)(I,IV,ISH) + (*WSh)(I,IV,ISH) * dt;
    }
    logfile(IV+1,") : ",(*XYZSh)(0,IV,ISH), " ", (*XYZSh)(1,IV,ISH), " ");
    logfile((*XYZSh)(2,IV,ISH), "\n");
   }
  }

  // de-allocate ZRoeSh
  freeArray();

  logfile.Close();
}

//----------------------------------------------------------------------------//

void MoveDps4Pg3D::setAddress()
{
  unsigned start;
  start = npoin->at(0) * PhysicsInfo::getnbDofMax() +
          PhysicsInfo::getnbShPointsMax() * PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDofMax();
  ZroeSh = new Array3D <double> (PhysicsInfo::getnbDofMax(),
                                 PhysicsInfo::getnbShPointsMax(),
                                 PhysicsInfo::getnbShMax(),
                                 &zroeVect->at(start));
}

//----------------------------------------------------------------------------//

void MoveDps4Pg3D::freeArray()
{
  delete ZroeSh;
}

//----------------------------------------------------------------------------//

void MoveDps4Pg3D::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector <double> > ("ZROE");
}

//----------------------------------------------------------------------------//

void MoveDps4Pg3D::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  gref = PhysicsData::getInstance().getData <double> ("GREF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
    PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  // nShockEdges =
  //   PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  nShockFaces =
    PhysicsData::getInstance().getData <vector <unsigned> > ("nShockFaces");
  typeSh =
    PhysicsData::getInstance().getData <vector <string> > ("TYPESH");
  WSh = PhysicsData::getInstance().getData <Array3D<double> > ("WSH");
  XYZSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYZSH");
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
