// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <iomanip>
#include "RemeshingSF/FixMshSps3D.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/Log.hh"
#include "Framework/Remeshing.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Array3D.hh"
#include "RemeshingSF/FindBEdg.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<FixMshSps3D, Remeshing> FixMshSps3DProv("FixMshSps3D");

//--------------------------------------------------------------------------//

FixMshSps3D::FixMshSps3D (const std::string& objectName)
 : Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

FixMshSps3D::~FixMshSps3D()
{
}

//--------------------------------------------------------------------------//

void FixMshSps3D::setup()
{
  LogToScreen(VERBOSE,"FixMshSps3D::setup() => start\n");

  LogToScreen(VERBOSE,"FixMshSps3D::setup() => end\n");
}

//--------------------------------------------------------------------------//

void FixMshSps3D::unsetup()
{
  LogToScreen(VERBOSE,"FixMshSps3D::unsetup()\n");
}

//--------------------------------------------------------------------------//

void FixMshSps3D::remesh()
{
  LogToScreen(INFO, "FixMshSps3D::remesh()\n");

  setMeshData();
  setPhysicsData();

  setAddress();

  logfile.Open(getClassName().c_str());

  createShockNodesList();

  // create downstream shock edges
  createDownShEdges();

  // for(unsigned ISPPNTS=0; ISPPNTS<(*nSpecPoints); ISPPNTS++) {
  //
  //  if       (typeSpecPoints->at(ISPPNTS) == "IPX"   ||
  //            typeSpecPoints->at(ISPPNTS) == "IPY"   ||
  //            typeSpecPoints->at(ISPPNTS) == "OPX"   ||
  //            typeSpecPoints->at(ISPPNTS) == "OPY"   ||
  //            typeSpecPoints->at(ISPPNTS) == "WPNRX" ||
  //            typeSpecPoints->at(ISPPNTS) == "WPNRY")
  //     { fixMeshForIPorOPorWPNR(ISPPNTS); }
  //
  //   else if (typeSpecPoints->at(ISPPNTS) == "RRX")
  //     { fixMeshForRRX(ISPPNTS); }
  //
  //   else if (typeSpecPoints->at(ISPPNTS) == "TP" ||
  //      typeSpecPoints->at(ISPPNTS) == "QP" ||
  //            typeSpecPoints->at(ISPPNTS) == "EP" ||
  //            typeSpecPoints->at(ISPPNTS) == "C" )
  //      {} // nothing to do
  //   else
  //      { cout << "Condition not implemented\n"; exit(1);}
  // }

  (*nbfacSh) = ibfac;

  // de-allocate dynamic arrays
  freeArray();

  logfile.Close();
}

//--------------------------------------------------------------------------//

void FixMshSps3D::createShockNodesList()
{
  ISHPlistu.resize( nShockPoints->at(0),
    // nShockPoints->at(0),
                    PhysicsInfo::getnbShMax());

  ISHPlistd.resize(nShockPoints->at(0),
    // nShockPoints->at(0),
                    PhysicsInfo::getnbShMax());

  unsigned ilist = npoin->at(0);

  for(unsigned ISH=0; ISH<PhysicsInfo::getnbShMax(); ISH++) {
  //  for(unsigned I=0; I<nShockPoints->at(0); I++) {
   for(unsigned I=0; I<nShockPoints->at(0); I++) {
    ++ilist;
    ISHPlistu(I,ISH) = ilist;
   }
  }

  for(unsigned ISH=0; ISH<PhysicsInfo::getnbShMax(); ISH++) {
  //  for(unsigned I=0; I<nShockPoints->at(0); I++) {
   for(unsigned I=0; I<nShockPoints->at(0); I++) {
    ++ilist;
    ISHPlistd(I,ISH) = ilist;
   }
  }
}

//--------------------------------------------------------------------------//

void FixMshSps3D::createDownShEdges()
{

  ibfac = nbfac->at(0);

  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
    unsigned n_counter_d=ISHPlistd(0,ISH);
  //  nShockEdges->at(ISH) = nShockPoints->at(ISH)-1;
   for (unsigned I=0; I<(nShockFaces->at(ISH)); I++) {
    (*bndfac)(0,ibfac) = (*ShFaces)(0,I,ISH)+n_counter_d;
    (*bndfac)(1,ibfac) = (*ShFaces)(1,I,ISH)+n_counter_d;
    (*bndfac)(2,ibfac) = (*ShFaces)(2,I,ISH)+n_counter_d;
    (*bndfac)(3,ibfac) = 10;
    ++ibfac;
   }

  unsigned n_counter_u=ISHPlistu(0,ISH);

   for (unsigned I=0; I<(nShockFaces->at(ISH)); I++) {
     (*bndfac)(0,ibfac) = (*ShFaces)(0,I,ISH)+n_counter_u;
     (*bndfac)(1,ibfac) = (*ShFaces)(1,I,ISH)+n_counter_u;
     (*bndfac)(2,ibfac) = (*ShFaces)(2,I,ISH)+n_counter_u;
     (*bndfac)(3,ibfac) = 10;
    ++ibfac;
   }

  }

}

//--------------------------------------------------------------------------//

void FixMshSps3D::fixMeshForIPorOPorWPNR(unsigned ISPPNTS)
{
  setShockIndeces(1,ISPPNTS);
  unsigned ishock = ISH.at(0)+1 ;

  FindBEdg lookForBE;

  // look for boundary edge for downstream point
  setShockPointCoor("Down",IP.at(0),ISH.at(0));
  logfile("\nTypeSpecPoints: ", typeSpecPoints->at(ISPPNTS), "\n");

  iedg1 = lookForBE.getBEdg(xsh, ysh);

  s1 = lookForBE.getS();

  if (iedg1==-1) {
   cout << "FixMshSps3D::error => failed matching shock points\n";
   cout << "                    look in the FixMshSps3D.log file for more info\n";
   logfile("Failed matching 1st shock point of the shock n.",ishock,"\n");
   logfile("Downstream point: xsh ",xsh,"ysh ",ysh,"\n");
   exit(1);}
  int iedg = iedg1+1; // only for log file printing, c++ indeces start from 0
  logfile(" s(1): ", s1, " ");
  logfile(xsh," ", ysh, " ", iedg);
  logfile ("\nShock point (1)", xsh, " ", ysh, " ");
  logfile ("falls within ", (*bndfac)(0,iedg1), " ", (*bndfac)(1,iedg1),"\n");

  // look for boundary edge for upstream point
  setShockPointCoor("Up",IP.at(0),ISH.at(0));
  iedg2 = lookForBE.getBEdg(xsh, ysh);
  s2 = lookForBE.getS();

  if (iedg2==-1) {
   cout << "FixMshSps3D::error => failed matching the shock points\n";
   cout << "                    look in the FixMshSps3D.log file for more info\n";
   logfile("Failed matching 1st shock point of the shock n.",ishock,"\n");
   logfile("Upstream point: xsh ",xsh,"ysh ",ysh,"\n");
   exit(1);}
  iedg = iedg2+1; // only for log file printing, c++ indeces start from 0
  logfile(" s(2): ", s2, " ");
  logfile(xsh," ", ysh, " ", iedg);
  logfile ("\nShock point (1)", xsh, " ", ysh, " ");
  logfile ("falls within ", (*bndfac)(0,iedg2), " ", (*bndfac)(1,iedg2),"\n");

  // check if the shock crosses a boundary point
  checkShockBndryEdgeCrossing();

  iedg = iedg1+1; // only for log file printing, c++ indeces start from 0
  logfile ("****************************");
  logfile ("    Shock: ", ishock, "\n");
  logfile ("    ", iedg, "\n");
  logfile ("****************************");


  // Split the existing edges of the background mesh
  // and create 2 new edges at  boundary
  // if the first point of the shock  is at the boundary
  if (iedg>0) { splitEdges();
    createNewEdges(IP.at(0), ISH.at(0)); }
}

//--------------------------------------------------------------------------//

void FixMshSps3D::fixMeshForRRX(unsigned ISPPNTS)
{
  // ISH.at(0) incident shock
  // ISH.at(1) reflected shock
  setShockIndeces(2,ISPPNTS);
  unsigned ishock = ISH.at(1)+1;

  FindBEdg lookForBE;

  // look for boundary edge for downstream point
  setShockPointCoor("Down",IP.at(1),ISH.at(1));
  iedg1 = lookForBE.getBEdg(xsh, ysh);
  s1 = lookForBE.getS();
  if (iedg1==-1) {
   logfile("Failed matching 1st shock point of the shock n.",ishock,"\n");
   exit(1);}
  int iedg = iedg1+1; // only for log file printing, c++ indeces start from 0
  logfile(" s(1): ", s1, " ");
  logfile(xsh," ", ysh, " ", iedg);
  logfile ("\nShock point (1)", xsh, " ", ysh, " ");
  logfile ("falls within ", (*bndfac)(0,iedg1), " ", (*bndfac)(1,iedg1),"\n");

  // look for boundary edge for upstream point
  setShockPointCoor("Up",IP.at(0),ISH.at(0));
  iedg2 = lookForBE.getBEdg(xsh, ysh);
  s2 = lookForBE.getS();
  if (iedg2==-1) {
   ishock = ISH.at(0)+1;
   logfile("Failed matching 1st shock point of the shock n.",ishock,"\n");
   exit(1);}
  iedg = iedg2+1; // only for log file printing, c++ indeces start from 0
  logfile(" s(2): ", s2, " ");
  logfile(xsh," ", ysh, " ", iedg);
  logfile ("\nShock point (1)", xsh, " ", ysh, " ");
  logfile ("falls within ", (*bndfac)(0,iedg2), " ", (*bndfac)(1,iedg2),"\n");

  // check if the shock crosses a boundary point
  checkShockBndryEdgeCrossing();

  iedg = iedg1+1; // c++ indeces start from 0
  logfile ("****************************");
  logfile ("    Shock: ", ishock, "\n");
  logfile ("    ", iedg, "\n");
  logfile ("****************************");

  // Split the existing edges of the background mesh
  // and create two new edges at  boundary
  // if the first point of the shock  is at the boundary
  if (iedg>0) { splitEdges();
                createNewEdges(IP.at(0), ISH.at(0),IP.at(1), ISH.at(1)); }
}

//--------------------------------------------------------------------------//

void FixMshSps3D::setShockIndeces(unsigned nbDiscontinuities,unsigned ISPPNTS)
{
  unsigned I;
  ISH.resize(nbDiscontinuities);
  IP.resize(nbDiscontinuities);
  for(unsigned i=0; i<nbDiscontinuities; i++) {
   ISH.at(i) = (*SHinSPPs)(0,i,ISPPNTS)-1; // c++ indeces start from 0
   I = (*SHinSPPs)(1,i,ISPPNTS) - 1;
   IP.at(i) = I * (nShockPoints->at(ISH.at(i))-1); // c++ indeces start from 0
  }
}

//--------------------------------------------------------------------------//

void FixMshSps3D::setShockPointCoor(string zone, unsigned IP, unsigned ISH)
{
  if (zone=="Up")        { xsh = (*XYZShu)(0,IP,ISH);
                           ysh = (*XYZShu)(1,IP,ISH); }
  else if (zone=="Down") { xsh = (*XYZShd)(0,IP,ISH);
                           ysh = (*XYZShd)(1,IP,ISH); }
  logfile ("xsh: ", xsh, " ysh: ", ysh,"\n");
}

//--------------------------------------------------------------------------//

void FixMshSps3D::checkShockBndryEdgeCrossing()
{
  if(s1<0 || s1>1 || s2 <0 || s2>1) {
   logfile ("s(1): ", s1, " s(2): ", s2, "out of bonds\n");
   exit(1);
  }
  int iedge1 = iedg1+1; // only for file log printing, c++ indeces start from 0
  int iedge2 = iedg2+1; // only for file log printing, c++ indeces start from 0
  if (iedg2 != iedg1) {
   logfile ("Shock points (1) (2) NOT on the same boundary edge\n");
   logfile (iedge1, " ", iedge2, "\n");
   logfile (s1, " ", s2, "\n");
  }
}

//--------------------------------------------------------------------------//

void FixMshSps3D::splitEdges()
{
  I1 = (*bndfac)(0,iedg1); // c++ indeces start from 0
  I2 = (*bndfac)(1,iedg1); // c++ indeces start from 0
  IBC = (*bndfac)(2,iedg1); // c++ indeces start from 0

  (*bndfac)(2,iedg1) = -IBC;
  int iedg = iedg1+1; // only for file log printing, c++ indeces start from 0
  logfile("Removing background edge ", iedg, " ");
  logfile(I1, " ", I2, "\n\n");
}

//--------------------------------------------------------------------------//

void FixMshSps3D::createNewEdges(unsigned IP, unsigned ISH)
{
  (*bndfac)(2,ibfac) = IBC;
  if (s1<s2) {
   (*bndfac)(0,ibfac)=I1;
   (*bndfac)(1,ibfac) = ISHPlistd(IP,ISH);
  }
  else {
   (*bndfac)(0,ibfac)=I1;
   (*bndfac)(1,ibfac) = ISHPlistu(IP,ISH);
  }

  ibfac++;
  (*bndfac)(2,ibfac) = IBC;
  if (s1<s2) {
   (*bndfac)(0,ibfac)=ISHPlistu(IP,ISH);
   (*bndfac)(1,ibfac) = I2;
  }
  else {
   (*bndfac)(0,ibfac)=ISHPlistd(IP,ISH);
   (*bndfac)(1,ibfac) = I2;
  }
  ibfac++;
}

//--------------------------------------------------------------------------//

void FixMshSps3D::createNewEdges(unsigned IP1, unsigned ISH1,
                               unsigned IP2, unsigned ISH2)
{
  (*bndfac)(2,ibfac) = IBC;
  if (s1<s2) {
   (*bndfac)(0,ibfac)=I1;
   (*bndfac)(1,ibfac) = ISHPlistd(IP2,ISH2);
  }
  else {
   (*bndfac)(0,ibfac)=I1;
   (*bndfac)(1,ibfac) = ISHPlistu(IP1,ISH1);
  }

  ibfac++;
  (*bndfac)(2,ibfac) = IBC;
  if (s1<s2) {
   (*bndfac)(0,ibfac)=ISHPlistu(IP1,ISH1);
   (*bndfac)(1,ibfac) = I2;
  }
  else {
   (*bndfac)(0,ibfac)=ISHPlistd(IP2,ISH2);
   (*bndfac)(1,ibfac) = I2;
  }
  ibfac++;
}

//--------------------------------------------------------------------------//

void FixMshSps3D::setAddress()
{
  unsigned start;
  start = 0;
  XYZ = new Array2D <double> (PhysicsInfo::getnbDim(), npoin->at(0),
                             &coorVect->at(start));
  // unsigned totsize = nbfac->at(0) + 2 *
  //                    PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShEdgesMax();
  // bndfac = new Array2D<int> (4,totsize,&bndfacVect->at(start));
  unsigned totsize = nbfac->at(0) + 2 *
                     PhysicsInfo::getnbShMax() * nShockFaces->at(0);
  bndfac = new Array2D<int> (4,totsize,&bndfacVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDim();
  XYZShu =
    new Array3D <double> (PhysicsInfo::getnbDim(),
                          nShockPoints->at(0),
                          PhysicsInfo::getnbShMax(),
                          &coorVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDim() +
          nShockPoints->at(0) *
          PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDim();
  XYZShd =
    new Array3D <double> (PhysicsInfo::getnbDim(),
                          nShockPoints->at(0),
                          PhysicsInfo::getnbShMax(),
                          &coorVect->at(start));
}

//--------------------------------------------------------------------------//

void FixMshSps3D::freeArray()
{
  delete bndfac; delete XYZ;
  delete XYZShu; delete XYZShd;
}

//--------------------------------------------------------------------------//

void FixMshSps3D::setMeshData()
{
  nbfac = MeshData::getInstance().getData <vector<unsigned> > ("NBFAC");
  nbfacSh = MeshData::getInstance().getData <unsigned> ("NBFACSH");
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  nelem = MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  coorVect = MeshData::getInstance().getData <vector<double> > ("COOR");
  bndfacVect = MeshData::getInstance().getData <vector <int> > ("BNDFAC");
}

//--------------------------------------------------------------------------//

void FixMshSps3D::setPhysicsData()
{
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nSpecPoints = PhysicsData::getInstance().getData <unsigned> ("nSpecPoints");
  nShockPoints =
      PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  nShockEdges =
      PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  typeSpecPoints =
      PhysicsData::getInstance().getData <vector <string> > ("TypeSpecPoints");
  SHinSPPs =
      PhysicsData::getInstance().getData <Array3D <unsigned> > ("SHinSPPs");
  ShFaces =
     PhysicsData::getInstance().getData <Array3D<unsigned> > ("ShFaces");
  nShockFaces =
     PhysicsData::getInstance().getData <vector<unsigned> > ("nShockFaces");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
