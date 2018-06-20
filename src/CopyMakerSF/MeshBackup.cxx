// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CopyMakerSF/MeshBackup.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsInfo.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<MeshBackup, CopyMaker>
  meshBackupProv("MeshBackup");

//--------------------------------------------------------------------------//

MeshBackup::MeshBackup(const std::string& objectName) :
  CopyMaker(objectName)
{
}

//--------------------------------------------------------------------------//

MeshBackup::~MeshBackup()
{
}

//--------------------------------------------------------------------------//

void MeshBackup::setup()
{
  LogToScreen(VERBOSE,"MeshBackup::setup() => start\n");

  LogToScreen(VERBOSE,"MeshBackup::setup() => end\n");
}

//--------------------------------------------------------------------------//

void MeshBackup::unsetup()
{
  LogToScreen(VERBOSE,"MeshBackup::unsetup()\n");
}

//--------------------------------------------------------------------------//

void MeshBackup::copy()
{
  LogToScreen(INFO,"MeshBackup::copy()\n");

  setMeshData();

  setAddress();

  // set size of the backup bndfac
  bndfacBackup->resize(4,nbfac->at(0) + 2 *		// 3 diventa 4
                       PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShEdgesMax());

  // set size of the backup nodptr
  nodptrBackup->resize(nbpoin->at(0),40);

  // set size of the backup nodcod
  nodcodBackup->resize(npoin->at(0)); 

  // make the backup of bndfac
  for(unsigned IFACE=0; IFACE<nbfac->at(0); IFACE++) {
   for(unsigned I=0; I<4; I++) {
    (*bndfacBackup)(I,IFACE) = (*bndfac)(I,IFACE);
   }
  }

  // make the backup of nodptr
  for(unsigned IBPOIN=0; IBPOIN<nbpoin->at(0); IBPOIN++) {
   for(unsigned I=0; I<10; I++) {
    (*nodptrBackup)(IBPOIN,I) = (*nodptr)(IBPOIN,I);
   }
  }

  // make the backup of nodcod
  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   nodcodBackup->at(IPOIN) = nodcod->at(IPOIN);
  }

  // make the backup of nbpoin
  // Index 2 refers to backup value
  nbpoin->at(2) = nbpoin->at(0);

  // make the backup of nbfac
  // Index 2 refers to backup value
  nbfac->at(2) = nbfac->at(0);

  // de-allocate the dynamic arrays
  freeArray();
}

//--------------------------------------------------------------------------//

void MeshBackup::setAddress()
{
  unsigned totsize = nbfac->at(0) + 2 *
                     PhysicsInfo::getnbShMax()*PhysicsInfo::getnbShEdgesMax();
  bndfac = new Array2D<int> (4,totsize,&bndfacVect->at(0));
  nodptr = new Array2D<int> (nbpoin->at(0),10, &nodptrVect->at(0));
}

//--------------------------------------------------------------------------//

void MeshBackup::freeArray()
{
  delete bndfac; delete nodptr;
}

//--------------------------------------------------------------------------//

void MeshBackup::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nbpoin = MeshData::getInstance().getData <vector<unsigned> > ("NBPOIN");
  nbfac = MeshData::getInstance().getData <vector<unsigned> > ("NBFAC");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  nodptrVect = MeshData::getInstance().getData <vector<int> >("NODPTR");
  bndfacVect = MeshData::getInstance().getData <vector<int> >("BNDFAC");
  nodcodBackup =
     MeshData::getInstance().getData <vector<int> >("NODCODBackup");
  nodptrBackup =
     MeshData::getInstance().getData <Array2D<int> >("NODPTRBackup");
  bndfacBackup = 
     MeshData::getInstance().getData <Array2D<int> >("BNDFACBackup");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
