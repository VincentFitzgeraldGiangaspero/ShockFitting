// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MeshGeneratorSF/ReSdwInfo3D.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/MeshGenerator.hh"
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
ObjectProvider<ReSdwInfo3D, MeshGenerator> readSdwInfo3DProv("ReSdwInfo3D");

//--------------------------------------------------------------------------//

ReSdwInfo3D::ReSdwInfo3D(const std::string& objectName) :
  MeshGenerator(objectName)
{
}

//--------------------------------------------------------------------------//

ReSdwInfo3D::~ReSdwInfo3D()
{
}

//--------------------------------------------------------------------------//

void ReSdwInfo3D::setup()
{
  LogToScreen(VERBOSE, "ReSdwInfo3D::setup() => start\n");
}

//--------------------------------------------------------------------------//

void ReSdwInfo3D::unsetup()
{
  LogToScreen(VERBOSE, "ReSdwInfo3D::unsetup()\n");
}

//--------------------------------------------------------------------------//

void ReSdwInfo3D::generate()
{
  LogToScreen(INFO, "ReadSdwInfo::generate()\n");

  logfile.Open(getClassName());
  // check if the sh00 file is already opened
  if(file.is_open()) {
   cout << "ReSdwInfo3D::warning => file sh00.dat seems to be already opened\n";}

  file.open(getInputFiles().c_str());



  readShockInfo();

  file.close();

  // de-allocate dynamic arrays
  freeArray();

  logfile.Close();
}

//--------------------------------------------------------------------------//

void ReSdwInfo3D::generate(string processingFile)
{
}

//--------------------------------------------------------------------------//

std::string ReSdwInfo3D::getInputFiles() const
{
  using namespace std;
  assert(m_inputFile.size()==1);
  string name = m_inputFile[0];
  return name;
}

//--------------------------------------------------------------------------//

void ReSdwInfo3D::readShockInfo()
{
  unsigned NSHE, IDUMMY;


  setMeshData();
   setPhysicsData();






  logfile("Open sh00.dat\n");

  file >> (*nShocks);
  logfile("Found n. ",(*nShocks)," shock/discontinuities\n");

  file.precision(18);

   unsigned ISH=0;
   unsigned iShock = ISH+1; // c++ indeces start from 0
   logfile("Shock/Discontinuity n. ",iShock,"\n");
   file >> (*nShockPoints)[ISH] >> (*typeSh)[ISH] ;
   file >> (*nShockFaces)[ISH];
   PhysicsInfo::setnbShPoint(nShockPoints->at(ISH));
   std::cout << "-PhysicsInfo::getnbShPointsMax() " << PhysicsInfo::getnbShPointsMax() <<endl;
   logfile("Kind of discontinuity: ",(*typeSh)[ISH],"\n");
   logfile("n. of points ",(*nShockPoints)[ISH],"\n");


   setAddress();
        setSize();


   for (unsigned K=0; K < (*nShockPoints)[ISH]; K++) {
    for (unsigned I=0; I <3; I++) {
     file >> (*XYZSh)(I,K,ISH);
     logfile((*XYZSh)(I,K,ISH)," ");
    }
    logfile("\n");
   }

   for (unsigned K=0; K < (*nShockPoints)[ISH]; K++) {
    for (unsigned I=0; I < (*ndof); I++) {

     file >>(*ZRoeShd)(I,K,ISH);
     logfile((*ZRoeShd)(I,K,ISH)," ");
    }
    logfile("\n");
   }

   for (unsigned K=0; K < (*nShockPoints)[ISH]; K++) {
    for (unsigned I=0; I < (*ndof); I++) {
     file >> (*ZRoeShu)(I,K,ISH);
     logfile((*ZRoeShu)(I,K,ISH)," ");
    }
    logfile("\n");
    (*NodCodSh)(K,ISH) = 10;
   }

   for (unsigned K=0; K < (*nShockPoints)[ISH]; K++) {
    for (unsigned I=0; I < (*ndof); I++) {
     (*ZRoeShuOld)(I,K,ISH) = (*ZRoeShu)(I,K,ISH);
     (*ZRoeShdOld)(I,K,ISH) = (*ZRoeShd)(I,K,ISH);
    }
   }

	unsigned dummy;

   for (unsigned IFACE=0; IFACE < (*nShockFaces)[ISH]; IFACE++) {
     for (unsigned IV=0; IV < 3; IV++){
       file >> dummy;
       (*ShFaces)(IV,IFACE,ISH)=dummy;
     }
   }



  // initialize NODCODSH which is part of NODCOD
  // If the code -99 is used this means no shock point
  // If the code 10 is used this means shock point
  for (unsigned ISH=0; ISH < PhysicsInfo::getnbShMax(); ISH++) {
    for (unsigned K=0; K < PhysicsInfo::getnbShPointsMax(); K++)
     {(*NodCodSh)(K,ISH) = -99;}
  }


  std::cout << "-Edge detection phase " << endl;


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

    std::cout << "-npoinedge " << npoinedge << " \n";

}

//--------------------------------------------------------------------------//

void ReSdwInfo3D::setSHinSPPs(unsigned NSHE, unsigned ISPPNTS)
{
  for (unsigned K=0; K<NSHE; K++) {
   file >> (*SHinSPPs)(0,K,ISPPNTS) >> (*SHinSPPs)(1,K,ISPPNTS);
   logfile((*SHinSPPs)(0,K,ISPPNTS)," ",(*SHinSPPs)(1,K,ISPPNTS),"\n" );
  }
}

//--------------------------------------------------------------------------//

void ReSdwInfo3D::setSize()
{
  nShockPoints->resize(PhysicsInfo::getnbShMax());
  nShockEdges->resize(PhysicsInfo::getnbShMax());
  typeSpecPoints->resize((*nSpecPoints));
  typeSh->resize(PhysicsInfo::getnbShMax());
  XYZSh->resize(PhysicsInfo::getnbDim(),
               PhysicsInfo::getnbShPointsMax(),
               PhysicsInfo::getnbShMax());
  ZRoeShuOld->resize(PhysicsInfo::getnbDofMax(),
                     PhysicsInfo::getnbShPointsMax(),
                     PhysicsInfo::getnbShMax());
  ZRoeShdOld->resize((*ndof),
                     PhysicsInfo::getnbShPointsMax(),
                     PhysicsInfo::getnbShMax());
  SHinSPPs->resize(2,5,PhysicsInfo::getnbSpecPointsMax());
  ShFaces->resize(3,nShockFaces->at(0),PhysicsInfo::getnbShMax());
  ShEdgePoints->resize(1,nShockPoints->at(0),PhysicsInfo::getnbShMax());
}

//--------------------------------------------------------------------------//

void ReSdwInfo3D::setAddress()
{
  unsigned start;
  start = npoin->at(0);
  NodCodSh = new Array2D <int> (PhysicsInfo::getnbShPointsMax(),
                                PhysicsInfo::getnbShMax(),
                                &nodcod->at(start));
  start = npoin->at(0)*PhysicsInfo::getnbDofMax();
  ZRoeShu =
    new Array3D <double> (PhysicsInfo::getnbDofMax(),
                          PhysicsInfo::getnbShPointsMax(),
                          PhysicsInfo::getnbShMax(),
                          &zroe->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDofMax() +
          PhysicsInfo::getnbShPointsMax() *
          PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDofMax();
  ZRoeShd =
    new Array3D <double> (PhysicsInfo::getnbDofMax(),
                          PhysicsInfo::getnbShPointsMax(),
                          PhysicsInfo::getnbShMax(),
                          &zroe->at(start));

//  ShEdgePoints =  new Array3D<unsigned>(1,
//		  	  	  	  	  	  	  	  	  PhysicsInfo::getnbShPointsMax(),
//		  	  	  	  	  	  	  	  	  PhysicsInfo::getnbShMax());

}

//--------------------------------------------------------------------------//

void ReSdwInfo3D::freeArray()
{
  delete NodCodSh; delete ZRoeShu; delete ZRoeShd;
}

//--------------------------------------------------------------------------//

void ReSdwInfo3D::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nodcod = MeshData::getInstance().getData <vector <int> > ("NODCOD");
  zroe = MeshData::getInstance().getData <vector <double> > ("ZROE");
}

//--------------------------------------------------------------------------//

void ReSdwInfo3D::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  nShockEdges =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  nShockFaces =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockFaces");

  ShEdgePoints = PhysicsData::getInstance().getData <Array3D <unsigned> > ("ShEDGEPoints");

  nSpecPoints = PhysicsData::getInstance().getData <unsigned> ("nSpecPoints");
  typeSpecPoints =
     PhysicsData::getInstance().getData <vector <string> > ("TypeSpecPoints");
  typeSh = PhysicsData::getInstance().getData <vector <string> > ("TYPESH");
  XYZSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYZSH");
  ZRoeShuOld =
       PhysicsData::getInstance().getData <Array3D <double> > ("ZROESHuOLD");
  ZRoeShdOld =
       PhysicsData::getInstance().getData <Array3D <double> > ("ZROESHdOLD");
  ShFaces =
       PhysicsData::getInstance().getData <Array3D <unsigned> > ("ShFaces");
  SHinSPPs =
       PhysicsData::getInstance().getData <Array3D <unsigned> > ("SHinSPPs");
}

//--------------------------------------------------------------------------//

} //namespace ShockFitting
