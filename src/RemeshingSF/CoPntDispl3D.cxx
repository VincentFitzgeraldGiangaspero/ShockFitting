// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cmath>
#include "RemeshingSF/CoPntDispl3D.hh"
#include "RemeshingSF/CoIntrPnt.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "SConfig/ObjectProvider.hh"

//vincent

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<CoPntDispl3D, Remeshing> coPointdispl3DProv("CoPntDispl3D");

//--------------------------------------------------------------------------//

CoPntDispl3D::CoPntDispl3D (const std::string& objectName) :
  Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

CoPntDispl3D::~CoPntDispl3D()
{
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::setup()
{
  LogToScreen(VERBOSE, "CoPntDispl3D::setup() => start\n");

  LogToScreen(VERBOSE, "CoPntDispl3D::setup() => end\n");
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::unsetup()
{
  LogToScreen(VERBOSE, "CoPntDispl3D::unsetup()\n");
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::remesh()
{
  LogToScreen(INFO, "CoPntDispl3D::remesh()\n");

  setMeshData();
  setPhysicsData();
  setAddress();

  logfile.Open(getClassName());

  // initialize nodcodsh values to -99
  setNodCodSh();

  // add second layer of shock nodes
  // place extra shock nodes on the shock normal
  addSecondShockLayer();

  // // correct the displacement in the boundary special points
  // for(unsigned ISPPNTS=0; ISPPNTS<(*nSpecPoints); ISPPNTS++) {
  //
  //  if (typeSpecPoints->at(ISPPNTS) == "IPX" ||
  //      typeSpecPoints->at(ISPPNTS) == "OPX" ||
  //      typeSpecPoints->at(ISPPNTS) == "WPNRX")
  //   { setCoorForIPXorOPXorWPNRX(ISPPNTS); }
  //
  //  else if (typeSpecPoints->at(ISPPNTS) == "IPY" ||
  //           typeSpecPoints->at(ISPPNTS) == "OPY" ||
  //           typeSpecPoints->at(ISPPNTS) == "WPNRY")
  //   { setCoorForIPYorOPYorWPNRY(ISPPNTS); }
  //
  //  else if (typeSpecPoints->at(ISPPNTS) == "TP")
  //   { setCoorForTP(ISPPNTS); }
  //
  //  else if (typeSpecPoints->at(ISPPNTS) == "RRX")
  //   { setCoorForRRX(ISPPNTS); }
  //
  //  else if (typeSpecPoints->at(ISPPNTS) == "QP")
  //   { setCoorForQP(ISPPNTS); }
  //
  //  else if (typeSpecPoints->at(ISPPNTS) == "EP")
  //   { setCoorForEP(ISPPNTS); }
  //
  //  else if (typeSpecPoints->at(ISPPNTS) == "C")
  //   { setCoorForC(ISPPNTS); }
  //
  //  else
  //   { logfile ("Condition not implemented\n");
  //     cout << "Condition not implemented\n";
  //     exit(1); }
  // }

  // de-allocate dynamic arrays
  freeArray();

  logfile.Close();
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::setNodCodSh()
{
  for(unsigned ISH=0; ISH<PhysicsInfo::getnbShMax(); ISH++) {
   for(unsigned I=0; I<nShockPoints->at(0); I++) {
    (*NodCodSh)(I,ISH) = -99;
   }
  }
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::addSecondShockLayer()
{
  for(unsigned ish=0; ish<(*nShocks); ish++) {
   for(unsigned i=0; i<nShockPoints->at(ish); i++) {
    dx = (*vShNor)(0,i,ish);
    dy = (*vShNor)(1,i,ish);
    dz = (*vShNor)(2,i,ish);
    (*XYZShu)(0,i,ish) = (*XYZSh)(0,i,ish) +
                        0.5 * (MeshData::getInstance().getEPS()) * dx;
    (*XYZShu)(1,i,ish) = (*XYZSh)(1,i,ish) +
                        0.5 * (MeshData::getInstance().getEPS()) * dy;
    (*XYZShu)(2,i,ish) = (*XYZSh)(2,i,ish) +
                        0.5 * (MeshData::getInstance().getEPS()) * dz;
    (*XYZShd)(0,i,ish) = (*XYZSh)(0,i,ish) -
                        0.5 * (MeshData::getInstance().getEPS()) * dx;
    (*XYZShd)(1,i,ish) = (*XYZSh)(1,i,ish) -
                        0.5 * (MeshData::getInstance().getEPS()) * dy;
    (*XYZShd)(2,i,ish) = (*XYZSh)(2,i,ish) -
                        0.5 * (MeshData::getInstance().getEPS()) * dz;
    (*NodCodSh)(i,ish) = 10;
   }
  }
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::setCoorForIPXorOPXorWPNRX(unsigned ISPPNTS)
{
  setShockIndeces(1,ISPPNTS);

  dx = (*vShNor)(0,IP.at(0),ISH.at(0));
  dy = (*vShNor)(1,IP.at(0),ISH.at(0));

  (*XYZShu)(0,IP.at(0),ISH.at(0)) = (*XYZSh)(0,IP.at(0),ISH.at(0))
                                    + 0.5 * (MeshData::getInstance().getEPS())/dx;
  (*XYZShu)(1,IP.at(0),ISH.at(0)) = (*XYZSh)(1,IP.at(0),ISH.at(0));
  (*XYZShd)(0,IP.at(0),ISH.at(0)) = (*XYZSh)(0,IP.at(0),ISH.at(0))
                                    - 0.5 * (MeshData::getInstance().getEPS())/dx;
  (*XYZShd)(1,IP.at(0),ISH.at(0)) = (*XYZSh)(1,IP.at(0),ISH.at(0));
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::setCoorForIPYorOPYorWPNRY(unsigned ISPPNTS)
{
  setShockIndeces(1,ISPPNTS);

  dx = (*vShNor)(0,IP.at(0),ISH.at(0));
  dy = (*vShNor)(1,IP.at(0),ISH.at(0));

  (*XYZShu)(0,IP.at(0),ISH.at(0)) = (*XYZSh)(0,IP.at(0),ISH.at(0));
  (*XYZShu)(1,IP.at(0),ISH.at(0)) = (*XYZSh)(1,IP.at(0),ISH.at(0))
                                   + 0.5 * (MeshData::getInstance().getEPS())/dy;
  (*XYZShd)(0,IP.at(0),ISH.at(0)) = (*XYZSh)(0,IP.at(0),ISH.at(0));
  (*XYZShd)(1,IP.at(0),ISH.at(0)) = (*XYZSh)(1,IP.at(0),ISH.at(0))
                                   - 0.5 * (MeshData::getInstance().getEPS())/dy;
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::setCoorForTP(unsigned ISPPNTS)
{
  // ISH.at(0) incident shock
  // ISH.at(1) reflected shock
  // ISH.at(2) mach stem
  // ISH.at(3) contact discontinuity
  setShockIndeces(4,ISPPNTS);

  // move incident shock
  moveDiscontinuity(IP.at(0),ISH.at(0));

  // move reflected shock
  moveDiscontinuity(IP.at(1),ISH.at(1));

  // move mach stem
  moveDiscontinuity(IP.at(2),ISH.at(2));

  // move contact discontinuity
  moveDiscontinuity(IP.at(3),ISH.at(3));

  // superimpose the upstream point of the incident shock to
  // the upstream point of the mach stem
  superimposeDiscPoints("up",IP.at(0),ISH.at(0),
                        "up",IP.at(2),ISH.at(2));

  // superimpose the downstream point of the incident shock to
  // the upstream point of the reflected shock
  superimposeDiscPoints("down",IP.at(0),ISH.at(0),
			"up",  IP.at(1),ISH.at(1));

  // superimpose the downstream point of the mach stem to
  // the downstream point of the contact discontinuity
  // (pay attention to the c.d. direction)
  superimposeDiscPoints("down",IP.at(2),ISH.at(2),
                        "down",IP.at(3),ISH.at(3));

  // superimpose downstream point of the reflected shock to
  // upstream point of the contact discontinuity
  // (pay attention to the c.d. direction)
  superimposeDiscPoints("down",IP.at(1),ISH.at(1),
		        "up",  IP.at(3),ISH.at(3));
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::superimposeDiscPoints(string discState1, unsigned IP1,
                                       unsigned ISH1,
				       string discState2, unsigned IP2,
                                       unsigned ISH2)
{
  Array3D <double>* XYZSh1 =
    new Array3D <double>(PhysicsInfo::getnbDim(),
                         nShockPoints->at(0),
                         PhysicsInfo::getnbShMax());
  Array3D <double>* XYZSh2 =
    new Array3D <double>(PhysicsInfo::getnbDim(),
                         nShockPoints->at(0),
                         PhysicsInfo::getnbShMax());

  if (discState1=="up")   {XYZSh1 = XYZShu;}
  if (discState1=="down") {XYZSh1 = XYZShd;}
  if (discState2=="up")   {XYZSh2 = XYZShu;}
  if (discState2=="down") {XYZSh2 = XYZShd;}

  for(unsigned K=0; K<PhysicsInfo::getnbDim(); K++) {
   dum = 0.5 * ((*XYZSh1)(K,IP1,ISH1)+(*XYZSh2)(K,IP2,ISH2));
   (*XYZSh1)(K,IP1,ISH1) = dum;
   (*XYZSh2)(K,IP2,ISH2) = dum;
  }
  delete XYZSh1; delete XYZSh2;
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::setCoorForRRX(unsigned ISPPNTS)
{
  // ISH.at(0) incident shock
  // ISH.at(1) reflected shock
  setShockIndeces(2,ISPPNTS);

  for (unsigned i=0; i<2; i++) {
   (*XYZShu)(0,IP.at(i),ISH.at(i)) = (*XYZSh)(0,IP.at(i),ISH.at(i))
				    + 0.5 * (MeshData::getInstance().getEPS())/dx;
   (*XYZShu)(1,IP.at(i),ISH.at(i)) = (*XYZSh)(1,IP.at(i),ISH.at(i));

   (*XYZShd)(0,IP.at(i),ISH.at(i)) = (*XYZSh)(0,IP.at(i),ISH.at(i))
                                    - 0.5 * (MeshData::getInstance().getEPS())/dx;
   (*XYZShd)(1,IP.at(i),ISH.at(i)) = (*XYZSh)(1,IP.at(i),ISH.at(i));
  }

  xc = getDownXVect(IP.at(0),ISH.at(0));
  yc = getDownYVect(IP.at(0),ISH.at(0));

  xs = getUpXVect(IP.at(1),ISH.at(1));
  ys = getUpYVect(IP.at(1),ISH.at(1));

  CoIntrPnt findNewCoor;
  findNewCoor.callCoIntrPnt(xc,yc,xs,ys);

  (*XYZShd)(0,IP.at(0),ISH.at(0)) = findNewCoor.getX();
  (*XYZShd)(1,IP.at(0),ISH.at(0)) = findNewCoor.getY();

  (*XYZShu)(0,IP.at(0),ISH.at(0)) = findNewCoor.getX();
  (*XYZShu)(1,IP.at(0),ISH.at(0)) = findNewCoor.getY();
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::setCoorForQP(unsigned ISPPNTS)
{
  // ISH.at(0) first incident shock
  // ISH.at(1) first reflected shock
  // ISH.at(2) second incident shock
  // ISH.at(3) second reflected shock
  // ISH.at(4) contact discontinuity
  setShockIndeces(5,ISPPNTS);

  // set the family which the first incident shock belongs to
  f1 = (*vShNor)(0,IP.at(0),ISH.at(0))* (*ZRoeShu)(3,IP.at(0),ISH.at(0))
      -(*vShNor)(1,IP.at(0),ISH.at(0))* (*ZRoeShu)(2,IP.at(0),ISH.at(0));
  f1 = copysign(1,f1);

  // set the family which the second incident shock belongs to
  f3 = (*vShNor)(0,IP.at(2),ISH.at(2))* (*ZRoeShu)(3,IP.at(2),ISH.at(2))
      -(*vShNor)(1,IP.at(2),ISH.at(2))* (*ZRoeShu)(2,IP.at(2),ISH.at(2));
  f3 = copysign(1,f1);

  // move point belongs to the first incident shock
  moveDiscontinuity(IP.at(0),ISH.at(0));

  // move point belongs to the first incident shock
  moveDiscontinuity(IP.at(1),ISH.at(1));

  // move point belongs to the second incident shock
  moveDiscontinuity(IP.at(2),ISH.at(2));

  // move point belongs to the second reflected shock
  moveDiscontinuity(IP.at(3),ISH.at(3));

  // move point belongs to contact discontinuity
  moveDiscontinuity(IP.at(4),ISH.at(4));

  CoIntrPnt findNewCoor;

  // superimpose point belongs to upstream state of
  // the first incident shock to point belongs to
  // downstream state of the second incident shock
  // if the incident shocks belong to different families
  if (f1>0) { xc = getUpXVect(IP.at(0),ISH.at(0));
   	      yc = getUpYVect(IP.at(0),ISH.at(0)); }

  else 	    { xc = getDownXVect(IP.at(0),ISH.at(0));
              yc = getDownYVect(IP.at(0),ISH.at(0)); }

  if (f3<0) { xs = getUpXVect(IP.at(2),ISH.at(2));
              ys = getUpYVect(IP.at(2),ISH.at(2)); }

  else      { xs = getDownXVect(IP.at(2),ISH.at(2));
              ys = getDownYVect(IP.at(2),ISH.at(2)); }

  findNewCoor.callCoIntrPnt(xc,yc,xs,ys);

  if (f1>0) { (*XYZShu)(0,IP.at(0),ISH.at(0)) = findNewCoor.getX();
              (*XYZShu)(1,IP.at(0),ISH.at(0)) = findNewCoor.getY(); }

  else      { (*XYZShd)(0,IP.at(0),ISH.at(0)) = findNewCoor.getX();
              (*XYZShd)(1,IP.at(0),ISH.at(0)) = findNewCoor.getY(); }

  if (f3<0) { (*XYZShu)(0,IP.at(2),ISH.at(2)) = findNewCoor.getX();
              (*XYZShu)(1,IP.at(2),ISH.at(2)) = findNewCoor.getY(); }

  else      { (*XYZShd)(0,IP.at(2),ISH.at(2)) = findNewCoor.getX();
  	      (*XYZShd)(1,IP.at(2),ISH.at(2)) = findNewCoor.getY(); }

  // superimpose point belongs to downstream state of first incident
  // shock to point belongs to downstream state of first
  // reflected shock
  if (f1>0) { xc = getDownXVect(IP.at(0),ISH.at(0));
              yc = getDownYVect(IP.at(0),ISH.at(0)); }

  else      { xc = getUpXVect(IP.at(0),ISH.at(0));
   	      yc = getUpYVect(IP.at(0),ISH.at(0)); }

  xs = getUpXVect(IP.at(1),ISH.at(1));
  ys = getUpYVect(IP.at(1),ISH.at(1));

  findNewCoor.callCoIntrPnt(xc,yc,xs,ys);

  if(f1>0) { (*XYZShd)(0,IP.at(0),ISH.at(0)) = findNewCoor.getX();
  	     (*XYZShd)(1,IP.at(0),ISH.at(0)) = findNewCoor.getY(); }

  else 	   { (*XYZShu)(0,IP.at(0),ISH.at(0)) = findNewCoor.getX();
             (*XYZShu)(1,IP.at(0),ISH.at(0)) = findNewCoor.getY(); }

  (*XYZShu)(0,IP.at(1),ISH.at(1)) = findNewCoor.getX();
  (*XYZShu)(1,IP.at(1),ISH.at(1)) = findNewCoor.getY();

  // superimpose point belongs to downstream state of second
  // incident shock to point belongs to upstream state of
  // second reflected shock
  if (f3>0) { xc = getUpXVect(IP.at(2),ISH.at(2));
   	      yc = getUpYVect(IP.at(2),ISH.at(2)); }

  else      { xc = getDownXVect(IP.at(2),ISH.at(2));
 	      yc = getDownYVect(IP.at(2),ISH.at(2)); }

  xs = getUpXVect(IP.at(3),ISH.at(3));
  ys = getUpYVect(IP.at(3),ISH.at(3));

  findNewCoor.callCoIntrPnt(xc,yc,xs,ys);

  if (f3<0) { (*XYZShd)(0,IP.at(2),ISH.at(2)) = findNewCoor.getX();
   	      (*XYZShd)(1,IP.at(2),ISH.at(2)) = findNewCoor.getY(); }

  else 	    { (*XYZShu)(0,IP.at(2),ISH.at(2)) = findNewCoor.getX();
   	      (*XYZShu)(1,IP.at(2),ISH.at(2)) = findNewCoor.getY(); }

  (*XYZShu)(0,IP.at(3),ISH.at(3)) = findNewCoor.getX();
  (*XYZShu)(1,IP.at(3),ISH.at(3)) = findNewCoor.getY();

  // superimpose point belongs to downstream state of first
  // reflected shock to point belongs to upstream state of
  // contact discontinuity
  xc = getDownXVect(IP.at(1),ISH.at(1));
  yc = getDownYVect(IP.at(1),ISH.at(1));

  xs = getUpXVect(IP.at(4),ISH.at(4));
  ys = getUpYVect(IP.at(4),ISH.at(4));

  findNewCoor.callCoIntrPnt(xc,yc,xs,ys);

  (*XYZShd)(0,IP.at(1),ISH.at(1)) = findNewCoor.getX();
  (*XYZShd)(1,IP.at(1),ISH.at(1)) = findNewCoor.getY();

  (*XYZShu)(0,IP.at(4),ISH.at(4)) = findNewCoor.getX();
  (*XYZShu)(1,IP.at(4),ISH.at(4)) = findNewCoor.getY();

  // superimpose point belongs to downstream state of
  // second reflected shock to point belong to downstream
  // state of contact discontinuity
  xc = getDownXVect(IP.at(3),ISH.at(3));
  yc = getDownYVect(IP.at(3),ISH.at(3));

  xs = getDownXVect(IP.at(4),ISH.at(4));
  ys = getDownYVect(IP.at(4),ISH.at(4));

  findNewCoor.callCoIntrPnt(xc,yc,xs,ys);

  (*XYZShd)(0,IP.at(3),ISH.at(3)) = findNewCoor.getX();
  (*XYZShd)(1,IP.at(4),ISH.at(4)) = findNewCoor.getY();
}

//--------------------------------------------------------------------------/e {

void CoPntDispl3D::setCoorForEP(unsigned ISPPNTS)
{
  setShockIndeces(1,ISPPNTS);

  // superimpose point belongs to downstream state to point
  // belongs to upstream state
  (*XYZShd)(0,IP.at(0),ISH.at(0)) = (*XYZSh)(0,IP.at(0),ISH.at(0));
  (*XYZShd)(1,IP.at(0),ISH.at(0)) = (*XYZSh)(1,IP.at(0),ISH.at(0));
  (*XYZShu)(0,IP.at(0),ISH.at(0)) = (*XYZSh)(0,IP.at(0),ISH.at(0));
  (*XYZShu)(1,IP.at(0),ISH.at(0)) = (*XYZSh)(1,IP.at(0),ISH.at(0));
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::setCoorForC(unsigned ISPPNTS)
{
  // ISH.at(0) incident shock
  // ISH.at(1) reflected shock
  setShockIndeces(2,ISPPNTS);

  // superimpose point belongs to downstream state to point
  // belongs to upstream state of connection point
  (*XYZShd)(0,IP.at(1),ISH.at(1)) = (*XYZShd)(0,IP.at(0),ISH.at(0));
  (*XYZShd)(1,IP.at(1),ISH.at(1)) = (*XYZShd)(1,IP.at(0),ISH.at(0));
  (*XYZShu)(0,IP.at(1),ISH.at(1)) = (*XYZShu)(0,IP.at(0),ISH.at(0));
  (*XYZShu)(1,IP.at(1),ISH.at(1)) = (*XYZShu)(1,IP.at(0),ISH.at(0));
}

//--------------------------------------------------------------------------//

std::vector <double> CoPntDispl3D::getUpXVect(unsigned IP, unsigned ISH)
{
  std::vector <double> x(2);
  x.at(0) = (*XYZShu)(0,IP,ISH);
  if (IP==0) { x.at(1) = (*XYZShu)(0,IP+1,ISH); }
  else       { x.at(1) = (*XYZShu)(0,IP-1,ISH); }
  return x;
}

//--------------------------------------------------------------------------//

std::vector <double> CoPntDispl3D::getUpYVect(unsigned IP, unsigned ISH)
{
  std::vector <double> y(2);
  y.at(0) = (*XYZShu)(1,IP,ISH);
  if (IP==0) { y.at(1) = (*XYZShu)(1,IP+1,ISH); }
  else       { y.at(1) = (*XYZShu)(1,IP-1,ISH); }
  return y;
}

//--------------------------------------------------------------------------//

std::vector <double> CoPntDispl3D::getDownXVect(unsigned IP, unsigned ISH)
{
  std::vector <double> x(2);
  x.at(0) = (*XYZShd)(0,IP,ISH);
  if (IP==0) { x.at(1) = (*XYZShd)(0,IP+1,ISH); }
  else       { x.at(1) = (*XYZShd)(0,IP-1,ISH); }
  return x;
}

//--------------------------------------------------------------------------//

std::vector <double> CoPntDispl3D::getDownYVect(unsigned IP, unsigned ISH)
{
  std::vector <double> y(2);
  y.at(0) = (*XYZShd)(1,IP,ISH);
  if (IP==0) { y.at(1) = (*XYZShd)(1,IP+1,ISH); }
  else       { y.at(1) = (*XYZShd)(1,IP-1,ISH); }
  return y;
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::moveDiscontinuity(unsigned IP, unsigned ISH)
{
  if(IP==1) {
   tx = (*XYZSh)(0,IP,ISH) - (*XYZSh)(0,IP-1,ISH);
   ty = (*XYZSh)(1,IP,ISH) - (*XYZSh)(1,IP-1,ISH);
  }
  else {
   tx = (*XYZSh)(0,IP-1,ISH) - (*XYZSh)(0,IP,ISH);
   ty = (*XYZSh)(1,IP-1,ISH) - (*XYZSh)(1,IP,ISH);
  }
  dum = sqrt(pow(tx,2)+pow(ty,2));
  tx = tx/dum;
  ty = ty/dum;

  (*XYZShu)(0,IP,ISH) = (*XYZShu)(0,IP,ISH) + (MeshData::getInstance().getEPS())*tx;
  (*XYZShu)(1,IP,ISH) = (*XYZShu)(1,IP,ISH) + (MeshData::getInstance().getEPS())*ty;

  (*XYZShd)(0,IP,ISH) = (*XYZShd)(0,IP,ISH) + (MeshData::getInstance().getEPS())*tx;
  (*XYZShd)(1,IP,ISH) = (*XYZShd)(1,IP,ISH) + (MeshData::getInstance().getEPS())*ty;
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::setShockIndeces(unsigned nbDiscontinuities, unsigned ISPPNTS)
{
  ISH.resize(nbDiscontinuities);
  IP.resize(nbDiscontinuities);
  for(unsigned i=0; i<nbDiscontinuities; i++) {
   ISH.at(i) = (*SHinSPPs)(0,i,ISPPNTS)-1; // c++ indeces start from 0
   I = (*SHinSPPs)(1,i,ISPPNTS) - 1;
   IP.at(i) = I * (nShockPoints->at(ISH.at(i))-1); // c++ indeces start from 0
  }
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::setAddress()
{
  unsigned start;
  start = npoin->at(0);
  NodCodSh = new Array2D <int> (nShockPoints->at(0),
                                PhysicsInfo::getnbShMax(),
                                &nodcod->at(start));
  start = npoin->at(0)*(*ndof);
  ZRoeShu =
    new Array3D <double> ((*ndof),
                          nShockPoints->at(0),
                          PhysicsInfo::getnbShMax(),
                          &zroe->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDim();
  XYZShu =
    new Array3D <double> (PhysicsInfo::getnbDim(),
                          nShockPoints->at(0),
                          PhysicsInfo::getnbShMax(),
                          &coor->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDim() +
          nShockPoints->at(0) *
          PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDim();
  XYZShd =
    new Array3D <double> (PhysicsInfo::getnbDim(),
                          nShockPoints->at(0),
                          PhysicsInfo::getnbShMax(),
                          &coor->at(start));
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::freeArray()
{
  delete ZRoeShu; delete NodCodSh;
  delete XYZShu; delete XYZShd;
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> >("NPOIN");
  nodcod = MeshData::getInstance().getData <vector<int> > ("NODCOD");
  zroe = MeshData::getInstance().getData <vector<double> > ("ZROE");
  coor = MeshData::getInstance().getData <vector<double> > ("COOR");
}

//--------------------------------------------------------------------------//

void CoPntDispl3D::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nSpecPoints =
        PhysicsData::getInstance().getData <unsigned> ("nSpecPoints");
  nShockPoints =
        PhysicsData::getInstance().getData <vector<unsigned> > ("nShockPoints");
  nShockEdges =
        PhysicsData::getInstance().getData <vector<unsigned> > ("nShockEdges");
  nShockFaces =
        PhysicsData::getInstance().getData <vector<unsigned> > ("nShockFaces");
  ShFaces =
         PhysicsData::getInstance().getData <Array3D<unsigned> > ("ShFaces");
  typeSh =
        PhysicsData::getInstance().getData <vector<string> > ("TYPESH");
  typeSpecPoints =
        PhysicsData::getInstance().getData <vector<string> > ("TypeSpecPoints");
  XYZSh =
        PhysicsData::getInstance().getData <Array3D<double> > ("XYZSH");
  SHinSPPs =
        PhysicsData::getInstance().getData <Array3D<unsigned> > ("SHinSPPs");
  vShNor =
        PhysicsData::getInstance().getData <Array3D<double> > ("VSHNOR");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
