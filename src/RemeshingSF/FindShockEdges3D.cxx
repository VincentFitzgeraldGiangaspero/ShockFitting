// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cmath>
#include "RemeshingSF/FindShockEdges3D.hh"
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
ObjectProvider<FindShockEdges3D, Remeshing> rdDpsEqProv("FindShockEdges3D");

//--------------------------------------------------------------------------//

FindShockEdges3D::FindShockEdges3D(const std::string& objectName) :
  Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

FindShockEdges3D::~FindShockEdges3D()
{
}

//--------------------------------------------------------------------------//

void FindShockEdges3D::setup()
{
  LogToScreen(VERBOSE, "FindShockEdges3D::setup() => start\n");

  LogToScreen(VERBOSE, "FindShockEdges3D::setup() => end\n");
}

//--------------------------------------------------------------------------//

void FindShockEdges3D::unsetup()
{
  LogToScreen(VERBOSE, "FindShockEdges3D::unsetup()\n");
}

//--------------------------------------------------------------------------//

void FindShockEdges3D::remesh()
{
  LogToScreen(INFO, "FindShockEdges3D::remesh()\n");

  setPhysicsData();
  setMeshData();

  logfile.Open(getClassName());
  // assign start pointers of Array2D and 3D
  setAddress();

  // resize vectors and arrays
  setSize();

  for (unsigned I=0; I<(*nShocks); I++){

   // compute length of shock edge
   computeShEdgeLength(I);

   // compute number of shock points of redistribution
   // and check the new number of shock points
   computeNbDistrShPoints(I);

   // compute distribution step
   computeDistrStep(I);

   // interpolation of new shock points
   // set the first and the last point
   newPointsInterp(I);

   // computation of internal points
   computeInterPoints(I);

   // rewrite the state arrays and indices
   rewriteValues(I);
  }

  // de-allocate dynamic arrays
  freeArray();

  logfile.Close();
}

//--------------------------------------------------------------------------//

void FindShockEdges3D::computeShEdgeLength(unsigned index)
{
  ISH=index;
  unsigned ish = ISH+1;
  logfile("Shock n. :", ish,"\n");
  Sh_ABSC.at(1) = 0;
  for (unsigned IV=1; IV<=nShockEdges->at(ISH); IV++) {
   unsigned I = IV+1;
   Sh_Edge_length = 0;
   for (unsigned K=0; K<PhysicsInfo::getnbDim(); K++) {
    Sh_Edge_length = Sh_Edge_length +
                     pow(((*XYZSh)(K,IV-1,ISH)-(*XYZSh)(K,IV,ISH)),2);
   }
   Sh_Edge_length = sqrt(Sh_Edge_length);
   Sh_ABSC.at(I) = Sh_ABSC.at(I-1) + Sh_Edge_length;
  }

  logfile("Number of points old distribution:", nShockPoints->at(ISH),"\n");
}

//--------------------------------------------------------------------------//

void FindShockEdges3D::computeNbDistrShPoints(unsigned index)
{
  ISH=index;
  nShockPoints_new = Sh_ABSC.at(nShockPoints->at(ISH)) /
                     (MeshData::getInstance().getDXCELL())+1;

  if (nShockPoints_new > PhysicsInfo::getnbShPointsMax()) {
   cout << "FindShockEdges3D::error => Too many shock points! Increase NPSHMAX in input.case\n";
   exit(1);}
}

//--------------------------------------------------------------------------//

void FindShockEdges3D::computeDistrStep(unsigned index)
{
  ISH=index;
  double dx = Sh_ABSC.at(nShockPoints->at(ISH))/(nShockPoints_new-1);
  Sh_ABSC_New.at(1) = 0;
  for (unsigned IV=1; IV<=nShockPoints_new; IV++) {
   unsigned I = IV+1;
   Sh_ABSC_New.at(I) = Sh_ABSC_New.at(I-1) + dx;
  }
  logfile ("Number of points new distribution:",nShockPoints_new, "\n");
}

//--------------------------------------------------------------------------//

void FindShockEdges3D::newPointsInterp(unsigned index)
{
  ISH = index;

   for (unsigned K=0; K<PhysicsInfo::getnbDim(); K++) {
    XYZSh_New(K,0) = (*XYZSh)(K,0,ISH);
   }

   for (unsigned K=0; K<(*ndof); K++) {
    ZRoeShu_New(K,0) = (*ZRoeShu)(K,0,ISH);
    ZRoeShd_New(K,0) = (*ZRoeShd)(K,0,ISH);
   }

   for (unsigned K=0; K<PhysicsInfo::getnbDim(); K++) {
    XYZSh_New(K,nShockPoints_new-1) = (*XYZSh)(K,nShockPoints->at(ISH)-1,ISH);

   }

   for (unsigned K=0; K<(*ndof); K++) {
    ZRoeShu_New(K,nShockPoints_new-1) =
                              (*ZRoeShu)(K,nShockPoints->at(ISH)-1,ISH);
    ZRoeShd_New(K,nShockPoints_new-1) =
                              (*ZRoeShd)(K,nShockPoints->at(ISH)-1,ISH);
   }
}

//--------------------------------------------------------------------------//

void FindShockEdges3D::computeInterPoints(unsigned index)
{
  double ds, dsj, alpha, beta;
  unsigned j1, i1;
  ISH=index;
  for(unsigned i=2; i<=nShockPoints_new-1; i++) {
   for (unsigned j=2; j<=nShockPoints->at(ISH); j++) {
    if ( (Sh_ABSC_New.at(i) >= Sh_ABSC.at(j-1)) &&
         (Sh_ABSC_New.at(i) <  Sh_ABSC.at(j)) ) {
     ds = Sh_ABSC.at(j)-Sh_ABSC.at(j-1);
     dsj = Sh_ABSC_New.at(i)-Sh_ABSC.at(j-1);
     alpha = dsj/ds;
     beta = 1-alpha;
     j1 = j-1; i1 = i-1;
     logfile("node = ",j1,"\n");
     logfile("node_new = ",i1,"\n");
     for (unsigned K=0; K<PhysicsInfo::getnbDim(); K++) {
      // XYZSh has the indeces that start from 0, therefore "i-1" index
      // is used
      XYZSh_New(K,i-1) = beta * (*XYZSh)(K,j-2,ISH) + alpha * (*XYZSh)(K,j-1,ISH);
     }
     for (unsigned K=0; K<(*ndof); K++) {
      // ZRoeShu_New has the indeces that start from 0, therefore "i-1" index
      // is used
      ZRoeShu_New(K,i-1) = beta * (*ZRoeShu)(K,j-2,ISH) +
                           alpha * (*ZRoeShu)(K,j-1,ISH);
      // ZRoeShd_New has the indeces that start from 0, therefore "i-1" index
      // is used
      ZRoeShd_New(K,i-1) = beta * (*ZRoeShd)(K,j-2,ISH) +
                           alpha * (*ZRoeShd)(K,j-1,ISH);
      unsigned m_K = K+1;
      logfile("ZROESHd_j-1(",m_K,"=",(*ZRoeShd)(K,j-2,ISH), "\n");
     } // K
    } // if
   } // j
  } // i
}

//--------------------------------------------------------------------------//

void FindShockEdges3D::rewriteValues(unsigned index)
{
  ISH = index;
  logfile("new shock point coordinates and quantities");

  for (unsigned I=0; I<nShockPoints_new; I++) {
   for(unsigned K=0; K<PhysicsInfo::getnbDim(); K++) {
    unsigned m_dim = K+1;
    (*XYZSh)(K,I,ISH)= XYZSh_New(K,I);
    logfile("XYZSh(",m_dim,") ",(*XYZSh)(K,I,ISH), "\n");
   }
   logfile("ZRoeShu(1) ",ZRoeShu_New(0,I), " ");
   logfile("ZRoeShd(1) ",ZRoeShd_New(0,I), " ");
   for(unsigned IV=0;IV<(*ndof);IV++) {
    (*ZRoeShu)(IV,I,ISH) = ZRoeShu_New(IV,I);
    (*ZRoeShd)(IV,I,ISH) = ZRoeShd_New(IV,I);
   }
   logfile("\n");
  }
  nShockPoints->at(ISH)=nShockPoints_new;
  nShockEdges->at(ISH) = nShockPoints_new-1;
}

//--------------------------------------------------------------------------//

void FindShockEdges3D::setAddress()
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

void FindShockEdges3D::setSize()
{
  XYZSh_New.resize(PhysicsInfo::getnbDim(),
                  PhysicsInfo::getnbShPointsMax() );
  ZRoeShu_New.resize(PhysicsInfo::getnbDofMax(),
                     PhysicsInfo::getnbShPointsMax() );
  ZRoeShd_New.resize(PhysicsInfo::getnbDofMax(),
                     PhysicsInfo::getnbShPointsMax() );
  // Sh_ABSC is filled with the indeces that start from 1
  // Sh_ABSC(1:NSHMAX+1)
  Sh_ABSC.resize(PhysicsInfo::getnbShPointsMax()+1);
  // Sh_ABSC_New is filled with the indeces that start from 1
  // Sh_ABSC_New(1:NSHMAX+1)
  Sh_ABSC_New.resize(PhysicsInfo::getnbShPointsMax()+1);
}

//--------------------------------------------------------------------------//

void FindShockEdges3D::freeArray()
{
  delete ZRoeShu; delete ZRoeShd;
}

//--------------------------------------------------------------------------//

void FindShockEdges3D::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints = PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  nShockEdges = PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  XYZSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYSH");
}

//--------------------------------------------------------------------------//

void FindShockEdges3D::setMeshData()
{
  zroeVect = MeshData::getInstance().getData <vector <double> >("ZROE");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
