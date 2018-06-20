// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShockDetectorSF/LaplacianSmoothing3D.hh"

#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/MeshData.hh"

//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

LaplacianSmoothing3D::LaplacianSmoothing3D(){

}


//--------------------------------------------------------------------------//

LaplacianSmoothing3D::~LaplacianSmoothing3D()
{
}

//--------------------------------------------------------------------------//


//--------------------------------------------------------------------------//

void LaplacianSmoothing3D::Smooth()
{
  LogToScreen(INFO, "LaplacianSmoothing3D::remesh()\n");

  setPhysicsData();
  setMeshData();

  // assign start pointers of Array2D and 3D
  setAddress();

  // resize vectors and arrays
  setSize();

  for (unsigned I=0; I<(*nShocks); I++){

  findShockEdges();

  LaplacianSmooth();

  WriteShockFile();

  }

  // de-allocate dynamic arrays
  freeArray();

}

//--------------------------------------------------------------------------//

void LaplacianSmoothing3D::findShockEdges()
{

  std::cout << "     => LaplacianSmoothing3D::findShockEdges() \n";

  unsigned ISH=0;

  std::vector <double> alphat(nShockPoints->at(ISH));

  std::vector <double> n1(3), n2(3), n3(3), diff(3);

  v1.resize(3);
  v2.resize(3);
  v3.resize(3);


  double det;

  const double pi = 3.1415926535897;

  double npoinedge=0;

  for (unsigned IFACE=0; IFACE < nShockFaces->at(ISH); IFACE++){

    verta=(*ShFaces)(0,IFACE,ISH);
    vertb=(*ShFaces)(1,IFACE,ISH);
    vertc=(*ShFaces)(2,IFACE,ISH);


    for (unsigned IV=0; IV<3; IV++){
      v1.at(IV)=(*XYZSh)(IV,verta,ISH);
      v2.at(IV)=(*XYZSh)(IV,vertb,ISH);
      v3.at(IV)=(*XYZSh)(IV,vertc,ISH);
    }

    // compute normals  n1
    for (unsigned IV=0; IV<3; IV++){
      diff.at(IV)=v2.at(IV)-v1.at(IV);
    }
    det =sqrt(diff.at(0)*diff.at(0)+diff.at(1)*diff.at(1)+diff.at(2)*diff.at(2));

    for (unsigned IV=0; IV<3; IV++){
      n1.at(IV)=diff.at(IV)/det;
    }

    // compute normals n2
    for (unsigned IV=0; IV<3; IV++){
      diff.at(IV)=v3.at(IV)-v1.at(IV);
    }
    det =sqrt(diff.at(0)*diff.at(0)+diff.at(1)*diff.at(1)+diff.at(2)*diff.at(2));

    for (unsigned IV=0; IV<3; IV++){
      n2.at(IV)=diff.at(IV)/det;
    }

    // compute normals n3
    for (unsigned IV=0; IV<3; IV++){
      diff.at(IV)=v3.at(IV)-v2.at(IV);
    }
    det =sqrt(diff.at(0)*diff.at(0)+diff.at(1)*diff.at(1)+diff.at(2)*diff.at(2));

    for (unsigned IV=0; IV<3; IV++){
      n3.at(IV)=diff.at(IV)/det;
    }

    // compute alphat
    alphat.at(verta)=alphat.at(verta)+acos(n1.at(0)*n2.at(0)+n1.at(1)*n2.at(1)+n1.at(2)*n2.at(2));
    alphat.at(vertb)=alphat.at(vertb)+acos(-n1.at(0)*n3.at(0)-n1.at(1)*n3.at(1)-n1.at(2)*n3.at(2));
    alphat.at(vertc)=alphat.at(vertc)+acos(n2.at(0)*n3.at(0)+n2.at(1)*n3.at(1)+n2.at(2)*n3.at(2));

  }// for IFACE

  for (unsigned NPOIN=0; NPOIN < nShockPoints->at(ISH); NPOIN++){
    if (alphat.at(NPOIN)<1.5*pi){
      npoinedge=npoinedge+1;
      (*ShEdgePoints)(0,NPOIN,ISH)=1;
    }
  }


  std::cout << "-nShockEdgePoints: " << npoinedge << " \n";

}

//--------------------------------------------------------------------------//

void LaplacianSmoothing3D::LaplacianSmooth(){

  std::cout << "     => LaplacianSmoothing3D::LaplacianSmooth \n";

  unsigned ISH=0;

  std::vector <double> d(3);

  std::vector <unsigned int> counter(nShockPoints->at(ISH));

  for (unsigned IFACE=0; IFACE<nShockFaces->at(ISH); IFACE++){

    verta=(*ShFaces)(0,IFACE,ISH);
    vertb=(*ShFaces)(1,IFACE,ISH);
    vertc=(*ShFaces)(2,IFACE,ISH);

    for (unsigned IV=0; IV<3; IV++){
      v1.at(IV)=(*XYZSh)(IV,verta,ISH);
      v2.at(IV)=(*XYZSh)(IV,vertb,ISH);
      v3.at(IV)=(*XYZSh)(IV,vertc,ISH);
    }

    if ((*ShEdgePoints)(0,verta,ISH)==1 & (*ShEdgePoints)(0,vertb,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=(1.0/2)*(v1.at(IV)+v2.at(IV));
      }
    }

    else if ((*ShEdgePoints)(0,verta,ISH)==1 & (*ShEdgePoints)(0,vertc,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=(1.0/2)*(v1.at(IV)+v3.at(IV));
      }
    }

    else if ((*ShEdgePoints)(0,vertb,ISH)==1 & (*ShEdgePoints)(0,vertc,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=(1.0/2)*(v2.at(IV)+v3.at(IV));
      }
    }

    else if ((*ShEdgePoints)(0,verta,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=v1.at(IV);
      }
    }

    else if ((*ShEdgePoints)(0,vertb,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=v2.at(IV);
      }
    }

    else if ((*ShEdgePoints)(0,vertc,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=v3.at(IV);
      }
    }

    else {
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=(1.0/3)*(v1.at(IV)+v2.at(IV)+v3.at(IV));
      }
    }

    for(unsigned IV=0; IV<3; IV++){
      XYZSh_New(IV,verta)=d.at(IV)+XYZSh_New(IV,verta);
      XYZSh_New(IV,vertb)=d.at(IV)+XYZSh_New(IV,vertb);
      XYZSh_New(IV,vertc)=d.at(IV)+XYZSh_New(IV,vertc);
    }
    // only in the first step
    counter.at(verta)=counter.at(verta)+1;
    counter.at(vertb)=counter.at(vertb)+1;
    counter.at(vertc)=counter.at(vertc)+1;

  }// for IFACE


  for (unsigned NPOIN=0; NPOIN < nShockPoints->at(ISH); NPOIN++){
    for (unsigned IV=0; IV<3; IV++){
      (*XYZSh)(IV,NPOIN,ISH)=XYZSh_New(IV,NPOIN)/counter.at(NPOIN);
    }
  }


  //2 step

  XYZSh_New.resize(PhysicsInfo::getnbDim(),
  PhysicsInfo::getnbShPointsMax());

  for (unsigned IFACE=0; IFACE<nShockFaces->at(ISH); IFACE++){

    verta=(*ShFaces)(0,IFACE,ISH);
    vertb=(*ShFaces)(1,IFACE,ISH);
    vertc=(*ShFaces)(2,IFACE,ISH);

    for (unsigned IV=0; IV<3; IV++){
      v1.at(IV)=(*XYZSh)(IV,verta,ISH);
      v2.at(IV)=(*XYZSh)(IV,vertb,ISH);
      v3.at(IV)=(*XYZSh)(IV,vertc,ISH);
    }

    if ((*ShEdgePoints)(0,verta,ISH)==1 & (*ShEdgePoints)(0,vertb,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=(1.0/2)*(v1.at(IV)+v2.at(IV));
      }
    }

    else if ((*ShEdgePoints)(0,verta,ISH)==1 & (*ShEdgePoints)(0,vertc,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=(1.0/2)*(v1.at(IV)+v3.at(IV));
      }
    }

    else if ((*ShEdgePoints)(0,vertb,ISH)==1 & (*ShEdgePoints)(0,vertc,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=(1.0/2)*(v2.at(IV)+v3.at(IV));
      }
    }

    else if ((*ShEdgePoints)(0,verta,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=v1.at(IV);
      }
    }

    else if ((*ShEdgePoints)(0,vertb,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=v2.at(IV);
      }
    }

    else if ((*ShEdgePoints)(0,vertc,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=v3.at(IV);
      }
    }

    else {
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=(1.0/3)*(v1.at(IV)+v2.at(IV)+v3.at(IV));
      }
    }

    for(unsigned IV=0; IV<3; IV++){
      XYZSh_New(IV,verta)=d.at(IV)+XYZSh_New(IV,verta);
      XYZSh_New(IV,vertb)=d.at(IV)+XYZSh_New(IV,vertb);
      XYZSh_New(IV,vertc)=d.at(IV)+XYZSh_New(IV,vertc);
    }
  }// for IFACE


  for (unsigned NPOIN=0; NPOIN < nShockPoints->at(ISH); NPOIN++){

    for (unsigned IV=0; IV<3; IV++){
      (*XYZSh)(IV,NPOIN,ISH)=XYZSh_New(IV,NPOIN)/counter.at(NPOIN);
    }

  }

  // 3 step

  XYZSh_New.resize(PhysicsInfo::getnbDim(),
  PhysicsInfo::getnbShPointsMax());

  for (unsigned IFACE=0; IFACE<nShockFaces->at(ISH); IFACE++){

    verta=(*ShFaces)(0,IFACE,ISH);
    vertb=(*ShFaces)(1,IFACE,ISH);
    vertc=(*ShFaces)(2,IFACE,ISH);

    for (unsigned IV=0; IV<3; IV++){
      v1.at(IV)=(*XYZSh)(IV,verta,ISH);
      v2.at(IV)=(*XYZSh)(IV,vertb,ISH);
      v3.at(IV)=(*XYZSh)(IV,vertc,ISH);
    }

    if ((*ShEdgePoints)(0,verta,ISH)==1 & (*ShEdgePoints)(0,vertb,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=(1.0/2)*(v1.at(IV)+v2.at(IV));
      }
    }

    else if ((*ShEdgePoints)(0,verta,ISH)==1 & (*ShEdgePoints)(0,vertc,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=(1.0/2)*(v1.at(IV)+v3.at(IV));
      }
    }

    else if ((*ShEdgePoints)(0,vertb,ISH)==1 & (*ShEdgePoints)(0,vertc,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=(1.0/2)*(v2.at(IV)+v3.at(IV));
      }
    }

    else if ((*ShEdgePoints)(0,verta,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=v1.at(IV);
      }
    }

    else if ((*ShEdgePoints)(0,vertb,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=v2.at(IV);
      }
    }

    else if ((*ShEdgePoints)(0,vertc,ISH)==1){
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=v3.at(IV);
      }
    }

    else {
      for (unsigned IV=0; IV<3; IV++){
        d.at(IV)=(1.0/3)*(v1.at(IV)+v2.at(IV)+v3.at(IV));
      }
    }

    for(unsigned IV=0; IV<3; IV++){
      XYZSh_New(IV,verta)=d.at(IV)+XYZSh_New(IV,verta);
      XYZSh_New(IV,vertb)=d.at(IV)+XYZSh_New(IV,vertb);
      XYZSh_New(IV,vertc)=d.at(IV)+XYZSh_New(IV,vertc);
    }
  }// for IFACE

  for (unsigned NPOIN=0; NPOIN < nShockPoints->at(ISH); NPOIN++){
    for (unsigned IV=0; IV<3; IV++){
      (*XYZSh)(IV,NPOIN,ISH)=XYZSh_New(IV,NPOIN)/counter.at(NPOIN);
    }
  }

  XYZSh_New.resize(PhysicsInfo::getnbDim(),
  PhysicsInfo::getnbShPointsMax());
}

//--------------------------------------------------------------------------//

void LaplacianSmoothing3D::WriteShockFile(){

  std::ofstream out_smooth("shock_smooth.off");
  out_smooth << "OFF\n" << nShockPoints->at(0) << " " << nShockFaces->at(0) << " 0\n";
  for (unsigned NPOIN=0; NPOIN<nShockPoints->at(0); NPOIN ++){
    out_smooth << (*XYZSh)(0,NPOIN,0) << " " << (*XYZSh)(1,NPOIN,0) << " " <<(*XYZSh)(2,NPOIN,0) << "\n";
  }
  for (unsigned IFACE=0; IFACE<nShockFaces->at(0); IFACE ++){
    out_smooth << "3 " << (*ShFaces)(0,IFACE,0) << " " << (*ShFaces)(1,IFACE,0) <<  " ";
    out_smooth << (*ShFaces)(2,IFACE,0) << " \n";
  }

}

//--------------------------------------------------------------------------//

void LaplacianSmoothing3D::setAddress()
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

void LaplacianSmoothing3D::setSize()
{
  XYZSh_New.resize(PhysicsInfo::getnbDim(),
                  PhysicsInfo::getnbShPointsMax());

  ShEdgePoints->resize( 1,
                                         PhysicsInfo::getnbShPointsMax(),
                                         PhysicsInfo::getnbShMax());
}

//--------------------------------------------------------------------------//

void LaplacianSmoothing3D::freeArray()
{
  delete ZRoeShu; delete ZRoeShd;
}

//--------------------------------------------------------------------------//

void LaplacianSmoothing3D::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints = PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  ShEdgePoints = PhysicsData::getInstance().getData <Array3D <unsigned> > ("ShEDGEPoints");
  XYZSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYZSH");
  ShFaces =
   PhysicsData::getInstance().getData <Array3D<unsigned> > ("ShFaces");
  nShockFaces =
    PhysicsData::getInstance().getData <vector<unsigned> > ("nShockFaces");
}

//--------------------------------------------------------------------------//

void LaplacianSmoothing3D::setMeshData()
{
  zroeVect = MeshData::getInstance().getData <vector <double> >("ZROE");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
