// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and oc/gpl.txt for the license text.

#include "WritingMeshSF/WriteBackTetgen.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<WriteBackTetgen, WritingMesh> 
  writeBackTetgenProv("WriteBackTetgen");

//--------------------------------------------------------------------------//

WriteBackTetgen::WriteBackTetgen(const std::string& objectName) :
  WritingMesh(objectName)
{
}

//--------------------------------------------------------------------------//

WriteBackTetgen::~WriteBackTetgen()
{
}

//--------------------------------------------------------------------------//

void WriteBackTetgen::setup()
{
  LogToScreen(VERBOSE,"WriteBackTetgen::setup() => start\n");

  LogToScreen(VERBOSE,"WriteBackTetgen::setup() => end\n");
}

//--------------------------------------------------------------------------//

void WriteBackTetgen::unsetup()
{
  LogToScreen(VERBOSE,"WriteBackTetgen::unsetup()\n");
}

//--------------------------------------------------------------------------//

void WriteBackTetgen::writeShockTecplot(unsigned I){
}

//--------------------------------------------------------------------------//

void WriteBackTetgen::write()
{
  LogToScreen(INFO,"WriteBackTetgen::write()\n");

  setMeshData();
  setPhysicsData();

  setAddress();

  *fnameBack = "na99";
  string dummyfile = *fnameBack+".node";

  std::cout << " ---> backup file : " << dummyfile <<" \n";

  file = fopen(dummyfile.c_str(),"w");

  unsigned ilist = npoin->at(0);
  unsigned iPoin;

  fprintf(file, "%u %s %u %s %u %s",ilist,"  ",PhysicsInfo::getnbDim(),"  ",*ndof,"  1\n");
  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   iPoin = IPOIN+1;
   if      (nodcod->at(IPOIN)==-1) {
    fprintf(file,"%u %s",iPoin," ");
    for(unsigned IA=0; IA<PhysicsInfo::getnbDim(); IA++) { 
     fprintf(file,"%.16f %s",(*XYZ)(IA,IPOIN),"  ");}
    for(unsigned IA=0; IA<(*ndof); IA++) { 
     fprintf(file,"%.16f %s",(*Zroe)(IA,IPOIN),"  ");}
    fprintf(file,"%s","0\n");
   }
   else if (nodcod->at(IPOIN)==-2) {
    fprintf(file,"%u %s",iPoin," ");
    for(unsigned IA=0; IA<PhysicsInfo::getnbDim(); IA++) { 
     fprintf(file,"%.16f %s",(*XYZ)(IA,IPOIN),"  ");}
    for(unsigned IA=0; IA<(*ndof); IA++) {  
     fprintf(file,"%.16f %s",(*Zroe)(IA,IPOIN),"  ");}
    fprintf(file,"%s","2\n");
   }
   else {
    fprintf(file,"%u %s",iPoin," ");
    for(unsigned IA=0; IA<PhysicsInfo::getnbDim(); IA++) { 
     fprintf(file,"%.16f %s",(*XYZ)(IA,IPOIN),"  ");}
    for(unsigned IA=0; IA<(*ndof); IA++) {  
     fprintf(file,"%.16f %s",(*Zroe)(IA,IPOIN),"  ");}
    fprintf(file,"%u %s",nodcod->at(IPOIN), "\n");
   }
  }
  
  fclose(file);  

  ///// writing also tecplot file



  std::cout << " ---> backup file 2 : cfbackup.plt " << endl;


    // writing file
    FILE* cfback;




    // write the .plt file
    cfback = fopen("cfbackup.plt", "w");
    std::cout << " here " << endl;

    fprintf(cfback,"%s","TITLE      =  Unstructured grid data\n");
    fprintf(cfback,"%s","VARIABLES  =  \"x0\" \"x1\" \"x2\" ");
    fprintf(cfback,"%s","\"p\" \"u\" \"v\" \"w\" \"T\"\n");

    std::cout << " here " << endl;

    fprintf(cfback,"%s %u","ZONE   T=\"P0 ZONE0 Tetra\", N=",ilist);
    fprintf(cfback,"%s %u %s",", E=",nelem->at(0),"F=FEPOINT, ET=TETRAHEDRON, SOLUTIONTIME=0\n");

    std::cout << " here " << endl;


    // write the grid-points coordinates and state
    for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
    for(unsigned IA=0; IA<PhysicsInfo::getnbDim(); IA++) {
    fprintf(cfback,"%20.16E %s",(*XYZ)(IA,IPOIN), " ");}

    for(unsigned K=0; K<(*ndof); K++) {
    fprintf(cfback,"%20.16E %s",(*Zroe)(K,IPOIN)," ");}
    fprintf(cfback,"%s","\n");
    }

    std::cout << " here " << endl;
    // write the mesh cell nodes
    for(unsigned IELEM=0;IELEM<nelem->at(0);IELEM++) {
    for(unsigned K=0;K<4; K++){
    fprintf(cfback,"%11i %s",(*celnod)(K,IELEM)," ");
    }
    fprintf(cfback,"%s","\n");
    }

    fclose(cfback);
    //



  // de-allocate dynamic arrays
  freeArray();
}

//--------------------------------------------------------------------------//

void WriteBackTetgen::setAddress()
{
  Zroe = new Array2D<double>(PhysicsInfo::getnbDofMax(),
                             npoin->at(0),
                             &zroeVect->at(0));
  XYZ = new Array2D<double>(PhysicsInfo::getnbDim(), npoin->at(0),
                           &coorVect->at(0));

  celnod = new Array2D<int>(4, nelem->at(0),
                            &celnodVect->at(0));
}

//--------------------------------------------------------------------------//

void WriteBackTetgen::freeArray()
{
  delete Zroe; delete XYZ; delete celnod;
}

//--------------------------------------------------------------------------//

void WriteBackTetgen::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nelem = MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  fnameBack = MeshData::getInstance().getData <string>("FNAMEBACK");
}

//--------------------------------------------------------------------------//

void WriteBackTetgen::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
