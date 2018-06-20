// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MeshGeneratorSF/ReadShockOFF3D.hh"
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
ObjectProvider<ReadShockOFF3D, MeshGenerator> ReadShockOFF3Dprov("ReadShockOFF3D");

//--------------------------------------------------------------------------//

ReadShockOFF3D::ReadShockOFF3D(const std::string& objectName) :
  MeshGenerator(objectName)
{
}

//--------------------------------------------------------------------------//

ReadShockOFF3D::~ReadShockOFF3D()
{
}

//--------------------------------------------------------------------------//

void ReadShockOFF3D::setup()
{
  LogToScreen(VERBOSE, "ReadShockOFF3D::setup() => start\n");
}

//--------------------------------------------------------------------------//

void ReadShockOFF3D::unsetup()
{
  LogToScreen(VERBOSE, "ReadShockOFF3D::unsetup()\n");
}

//--------------------------------------------------------------------------//

void ReadShockOFF3D::generate()
{
  LogToScreen(INFO, "ReadSdwInfo::generate()\n");

  logfile.Open(getClassName());
  // check if the sh00 file is already opened
  if(file.is_open()) {
   cout << "ReadShockOFF3D::warning => file sh00.dat seems to be already opened\n";}

  file.open(getInputFiles().c_str());


  readShockInfo();

  file.close();

  // de-allocate dynamic arrays
  freeArray();

  logfile.Close();
}

//--------------------------------------------------------------------------//

void ReadShockOFF3D::generate(string processingFile)
{
}

//--------------------------------------------------------------------------//

std::string ReadShockOFF3D::getInputFiles() const
{
  using namespace std;
  assert(m_inputFile.size()==1);
  string name = m_inputFile[0];
  return name;
}

//--------------------------------------------------------------------------//

void ReadShockOFF3D::readShockInfo()
{
	unsigned NSHE, IDUMMY;
	double dummy;
	std::string dummy_s;

	  setMeshData();
	  setPhysicsData();





	  logfile("Open sh00.off\n");

	  file >>dummy_s;
	  (*nShocks)=1;
	  logfile("Found n. ",(*nShocks)," shock/discontinuities\n");

	  file.precision(18);

	   unsigned ISH=0;
	   unsigned iShock = ISH+1; // c++ indeces start from 0
	   logfile("Shock/Discontinuity n. ",iShock,"\n");
	   file >> (*nShockPoints)[ISH];
	   (*typeSh)[ISH]='S';
	   file >> (*nShockFaces)[ISH];
	   PhysicsInfo::setnbShPoint(nShockPoints->at(ISH));
	   logfile("Kind of discontinuity: ",(*typeSh)[ISH],"\n");
	   logfile("n. of points ",(*nShockPoints)[ISH],"\n");

	   setAddress();
	   setSize();

	   std::cout << "-nShockPoints->at(ISH) " << nShockPoints->at(ISH) <<"\n";
	   file >> dummy;
	   for (unsigned K=0; K < (*nShockPoints)[ISH]; K++) {
	    for (unsigned I=0; I <3; I++) {
	     file >> (*XYZSh)(I,K,ISH);
	     logfile((*XYZSh)(I,K,ISH)," ");
	    }
	    logfile("\n");
	   }

	   for (unsigned K=0; K < (*nShockPoints)[ISH]; K++) {
	    for (unsigned I=0; I < (*ndof); I++) {
	    	//file >>(*ZRoeShd)(I,K,ISH);
	    	//file>>dummy;
	     (*ZRoeShd)(I,K,ISH)=100.0;
	     logfile((*ZRoeShd)(I,K,ISH)," ");
	    }
	    logfile("\n");
	   }

	   for (unsigned K=0; K < (*nShockPoints)[ISH]; K++) {
	    for (unsigned I=0; I < (*ndof); I++) {
	    	//file >>(*ZRoeShu)(I,K,ISH);
	    	//file>> dummy;
	     (*ZRoeShu)(I,K,ISH)=100.0;
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

		unsigned dummy_uns;

	   for (unsigned IFACE=0; IFACE < (*nShockFaces)[ISH]; IFACE++) {
		   file >> dummy_uns;
	     for (unsigned IV=0; IV < 3; IV++){
	       file >> dummy_uns;
	       (*ShFaces)(IV,IFACE,ISH)=dummy_uns;
	     }
	   }



	  // initialize NODCODSH which is part of NODCOD
	  // If the code -99 is used this means no shock point
	  // If the code 10 is used this means shock point
	  for (unsigned ISH=0; ISH < PhysicsInfo::getnbShMax(); ISH++) {
	    for (unsigned K=0; K < PhysicsInfo::getnbShPointsMax(); K++)
	     {(*NodCodSh)(K,ISH) = -99;}
	  }
}

//--------------------------------------------------------------------------//

void ReadShockOFF3D::setSHinSPPs(unsigned NSHE, unsigned ISPPNTS)
{
  for (unsigned K=0; K<NSHE; K++) {
   file >> (*SHinSPPs)(0,K,ISPPNTS) >> (*SHinSPPs)(1,K,ISPPNTS);
   logfile((*SHinSPPs)(0,K,ISPPNTS)," ",(*SHinSPPs)(1,K,ISPPNTS),"\n" );
  }
}

//--------------------------------------------------------------------------//

void ReadShockOFF3D::setSize()
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
}

//--------------------------------------------------------------------------//

void ReadShockOFF3D::setAddress()
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

}

//--------------------------------------------------------------------------//

void ReadShockOFF3D::freeArray()
{
  delete NodCodSh; delete ZRoeShu; delete ZRoeShd;
}

//--------------------------------------------------------------------------//

void ReadShockOFF3D::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nodcod = MeshData::getInstance().getData <vector <int> > ("NODCOD");
  zroe = MeshData::getInstance().getData <vector <double> > ("ZROE");
}

//--------------------------------------------------------------------------//

void ReadShockOFF3D::setPhysicsData()
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
