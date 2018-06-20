// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MeshGeneratorSF/ReadShockMESH3D.hh"
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
ObjectProvider<ReadShockMESH3D, MeshGenerator> ReadShockMESH3DProv("ReadShockMESH3D");

//--------------------------------------------------------------------------//

ReadShockMESH3D::ReadShockMESH3D(const std::string& objectName) :
  MeshGenerator(objectName)
{
}

//--------------------------------------------------------------------------//

ReadShockMESH3D::~ReadShockMESH3D()
{
}

//--------------------------------------------------------------------------//

void ReadShockMESH3D::setup()
{
  LogToScreen(VERBOSE, "ReadShockMESH3D::setup() => start\n");
}

//--------------------------------------------------------------------------//

void ReadShockMESH3D::unsetup()
{
  LogToScreen(VERBOSE, "ReadShockMESH3D::unsetup()\n");
}

//--------------------------------------------------------------------------//

void ReadShockMESH3D::generate()
{
  LogToScreen(INFO, "ReadSdwInfo::generate()\n");

  logfile.Open(getClassName());
  // check if the sh00 file is already opened
  if(file.is_open()) {
   cout << "ReadShockMESH3D::warning => file sh00.dat seems to be already opened\n";}

  file.open(getInputFiles().c_str());


  readShockInfo();

  file.close();

  // de-allocate dynamic arrays
  freeArray();

  logfile.Close();
}

//--------------------------------------------------------------------------//

void ReadShockMESH3D::generate(string processingFile)
{
}

//--------------------------------------------------------------------------//

std::string ReadShockMESH3D::getInputFiles() const
{
  using namespace std;
  assert(m_inputFile.size()==1);
  string name = m_inputFile[0];
  return name;
}

//--------------------------------------------------------------------------//

void ReadShockMESH3D::readShockInfo()
{
  unsigned NSHE, IDUMMY;

  setMeshData();
  setPhysicsData();




  unsigned ISH=0;
  logfile("Open shock3D.o.mesh \n");

  (*nShocks)=1;
  (*typeSh)[0]='S';
  logfile("Found n. ",(*nShocks)," shock/discontinuities\n");





  std::cout << "     => ReadShockMESH3D::readMESHfile()\n";

    std::string dummyfile;
    std::string dummy;

    // reading file
    std::ifstream file;

    dummyfile = "shock3D.o.mesh";

    file.open(dummyfile.c_str());

    while(dummy != "Vertices"){
      file >> dummy;
    }

    file >> nShockPoints->at(0);
    std::cout << "-nShockPoints: " << nShockPoints->at(0) << '\n';
    PhysicsInfo::setnbShPoint(nShockPoints->at(0));


    XYZSh->resize(PhysicsInfo::getnbDim(),
                 PhysicsInfo::getnbShPointsMax(),
                 PhysicsInfo::getnbShMax());

    for (unsigned ISHPOIN=0; ISHPOIN<nShockPoints->at(0); ISHPOIN++) {
      for (unsigned IV=0; IV<3; IV++){
        file >> (*XYZSh)(IV,ISHPOIN,0);            // xyz
      }
      file >> dummy; //1
    }

    file >>dummy;

    file >> nShockFaces->at(0);
    std::cout << "-nShockFaces: " << nShockFaces->at(0) << '\n';

    ShFaces->resize(3,
                  nShockFaces->at(0),
                  PhysicsInfo::getnbShMax());

    ZRoeShu.resize((*ndof),
            nShockPoints->at(0),
            PhysicsInfo::getnbShMax());

    ZRoeShd.resize((*ndof),
                nShockPoints->at(0),
                PhysicsInfo::getnbShMax());


    for (unsigned IFACE=0; IFACE<nShockFaces->at(0); IFACE++) {
      for (unsigned INODE=0; INODE<3; INODE++){
        file >> (*ShFaces)(INODE,IFACE,0);            // node
        (*ShFaces)(INODE,IFACE,0)=(*ShFaces)(INODE,IFACE,0)-1;
      }
      file >> dummy; //1
    }

    file >>dummy;


    setAddress();
          setSize();


    	  // initialize NODCODSH which is part of NODCOD
    	  // If the code -99 is used this means no shock point
    	  // If the code 10 is used this means shock point
    	  for (unsigned ISH=0; ISH < PhysicsInfo::getnbShMax(); ISH++) {
    	    for (unsigned K=0; K < PhysicsInfo::getnbShPointsMax(); K++)
    	     {(*NodCodSh)(K,ISH) = -99;}
    	  }

    	  for (unsigned K=0; K < nShockPoints->at(0); K++) {
    	 	    for (unsigned I=0; I < (*ndof); I++) {

    	 	     ZRoeShd(I,K,ISH)=200.0;
    	 	    ZRoeShu(I,K,ISH)=100.0;
    	 	     logfile(ZRoeShd(I,K,ISH)," ");
    	 	    logfile(ZRoeShu(I,K,ISH)," ");
    	 	    }
    	 	    logfile("\n");
    	 	   }

    	  for (unsigned K=0; K < nShockPoints->at(0); K++) {
    		    for (unsigned I=0; I < (*ndof); I++) {
    		     (*ZRoeShuOld)(I,K,ISH) = ZRoeShu(I,K,ISH);
    		     (*ZRoeShdOld)(I,K,ISH) = ZRoeShd(I,K,ISH);
    		    }
    		   }

}

//--------------------------------------------------------------------------//

void ReadShockMESH3D::setSHinSPPs(unsigned NSHE, unsigned ISPPNTS)
{
  for (unsigned K=0; K<NSHE; K++) {
   file >> (*SHinSPPs)(0,K,ISPPNTS) >> (*SHinSPPs)(1,K,ISPPNTS);
   logfile((*SHinSPPs)(0,K,ISPPNTS)," ",(*SHinSPPs)(1,K,ISPPNTS),"\n" );
  }
}

//--------------------------------------------------------------------------//

void ReadShockMESH3D::setSize()
{
  //nShockPoints->resize(PhysicsInfo::getnbShMax());
  nShockEdges->resize(PhysicsInfo::getnbShMax());
  //typeSpecPoints->resize((*nSpecPoints));
  //typeSh->resize(PhysicsInfo::getnbShMax());
  //XYZSh->resize(PhysicsInfo::getnbDim(),
  //             PhysicsInfo::getnbShPointsMax(),
  //             PhysicsInfo::getnbShMax());
  ZRoeShuOld->resize(PhysicsInfo::getnbDofMax(),
                     PhysicsInfo::getnbShPointsMax(),
                     PhysicsInfo::getnbShMax());
  ZRoeShdOld->resize((*ndof),
                     PhysicsInfo::getnbShPointsMax(),
                     PhysicsInfo::getnbShMax());
  SHinSPPs->resize(2,5,PhysicsInfo::getnbSpecPointsMax());
  //ShFaces->resize(3,nShockFaces->at(0),PhysicsInfo::getnbShMax());
}

//--------------------------------------------------------------------------//

void ReadShockMESH3D::setAddress()
{
  unsigned start;
  start = npoin->at(0);
  NodCodSh = new Array2D <int> (PhysicsInfo::getnbShPointsMax(),
                                PhysicsInfo::getnbShMax(),
                                &nodcod->at(start));

}

//--------------------------------------------------------------------------//

void ReadShockMESH3D::freeArray()
{
  delete NodCodSh;
}

//--------------------------------------------------------------------------//

void ReadShockMESH3D::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nodcod = MeshData::getInstance().getData <vector <int> > ("NODCOD");
  zroe = MeshData::getInstance().getData <vector <double> > ("ZROE");
}

//--------------------------------------------------------------------------//

void ReadShockMESH3D::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  nShockEdges =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
 nShockFaces =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockFaces");
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
