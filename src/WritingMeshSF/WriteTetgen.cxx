// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and oc/gpl.txt for the license text.

#include "WritingMeshSF/WriteTetgen.hh"
#include "Framework/Log.hh"
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
ObjectProvider<WriteTetgen, WritingMesh> writeTetgenProv("WriteTetgen");

//--------------------------------------------------------------------------//

WriteTetgen::WriteTetgen(const std::string& objectName) :
  WritingMesh(objectName)
{
}

//--------------------------------------------------------------------------//

WriteTetgen::~WriteTetgen()
{
}

//--------------------------------------------------------------------------//

void WriteTetgen::setup()
{
  LogToScreen(VERBOSE, "WriteTetgen::setup() => start\n");

  LogToScreen(VERBOSE, "WriteTetgen::setup() => end\n");
}

//--------------------------------------------------------------------------//

void WriteTetgen::unsetup()
{
  LogToScreen(VERBOSE, "WriteTetgen::unsetup()\n");
}

//--------------------------------------------------------------------------//
void WriteTetgen::writeShockTecplot(unsigned I){

}

//--------------------------------------------------------------------------//

void WriteTetgen::write()
{
  LogToScreen(INFO, "WriteTetgen::write()\n");

  setMeshData();
  setPhysicsData();

  setAddress();

  string dummyfile;

  fname->str(string());

  unsigned nbDig = 0;
  unsigned dummyIstep = MeshData::getInstance().getIstep();

  while(dummyIstep>0) { dummyIstep/=10; nbDig++; }
  *fname << setw(7-nbDig) << setfill('0') << left << string("na").c_str();
  *fname << MeshData::getInstance().getIstep();

  // write node file
  dummyfile = fname->str()+".node";

  file = fopen(dummyfile.c_str(),"w");

  std::cout << "-file name generated  " << dummyfile << "\n";

  ilist = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
                             nShockPoints->at(0);

  // ilist = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
  //                            nShockPoints->at(0);


  // M02M1 and M12M0 are filled with indeces that start
  // from 1 to NPOIN+2*NSHMAX*NPSHMAX+1
  M02M1->resize(ilist+1); // c++ indeces start from 0
  M12M0->resize(ilist+1); // c++ indeces start from 0

  TNPOIN = 0;

  // set map vector for nodcod
  setMapVectorForNodcod();

  icount = npoin->at(0);

  double m_ZROE_UP;
  double m_ZROE_DOWN;
  std::cout << "-Shock Points State check \n";
  bool error =false;

  for (unsigned IV=0; IV<nShockPoints->at(0);IV++) {
    if (((*ZRoeShu)(0,IV,0)<0)||((*ZRoeShd)(0,IV,0)<0)){
      std::cout << "Shock Point n  " << IV  << " with negative pressure \n";
      error =  true;
      if ((*ShEdgePoints)(0,IV,0)==1){
        std::cout << "Edge Point !!!  \n";
      }
    }
    if ((*ZRoeShu)(0,IV,0)>(*ZRoeShd)(0,IV,0)){
      std::cout << "Shock Point n  " << IV  << " with swapped UP and DOWN \n";
      std::cout << "P UP " << (*ZRoeShu)(0,IV,0) << '\n';
      std::cout << "P DOWN " << (*ZRoeShd)(0,IV,0) << '\n';
      error =  true;
      for (unsigned ROE=0; ROE<(*ndof); ROE++){
        m_ZROE_UP=(*ZRoeShu)(ROE,IV,0);
        m_ZROE_DOWN=(*ZRoeShd)(ROE,IV,0);
        (*ZRoeShu)(ROE,IV,0)=m_ZROE_DOWN;
        (*ZRoeShd)(ROE,IV,0)= m_ZROE_UP;
      }
    }
  }

  if (error){exit(1);}

  // set map vector for NodCodSh
  setMapVectorForNodcodSh();

  // set map vector for NodCodSh
  setMapVectorForNodcodSh();

  ilist = TNPOIN;


  fprintf(file,"%u %s %u %s %u %s",ilist," ",PhysicsInfo::getnbDim()," ",*ndof," 1\n");

  // write mesh points coordinates and states on Triangle file
  writeMeshVariables();

  icount = npoin->at(0);

  // write upstream shock points coordinates and states on tetgen file
  writeUpstreamStatus();

  // write downstream shock points coordinates and states on tetgen file
  writeDownstreamStatus();

  fclose(file);


  // write poly file
  dummyfile = fname->str()+".poly";
  file = fopen(dummyfile.c_str(),"w");

  std::cout << "-file name generated  " << dummyfile << "\n";


  fprintf(file,"%s %u %s", "0 ",PhysicsInfo::getnbDim()," 0 1 \n");

  // set map vector for bndfac and write bndfac on poly file
  writeBndfac();


  // compute number of holes
  // computenbHoles();

  // de-allocate dynamic arrays
  freeArray();

  fclose(file);


}

//--------------------------------------------------------------------------//

void WriteTetgen::setMapVectorForNodcod()
{
  for (unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   if (nodcod->at(IPOIN)>=0) { ++TNPOIN;
                               M02M1->at(IPOIN+1) = TNPOIN;
                               M12M0->at(TNPOIN) = IPOIN+1; }
  }

}

//--------------------------------------------------------------------------//

void WriteTetgen::setMapVectorForNodcodSh()
{
  for (unsigned ISH=0; ISH<PhysicsInfo::getnbShMax(); ISH++) {
  //  for(unsigned I=0; I<nShockPoints->at(0); I++) {
   for(unsigned I=0; I<nShockPoints->at(0); I++) {
    ++icount;
    if ((*NodCodSh)(I,ISH)==10) { ++TNPOIN;
                                    M02M1->at(icount) = TNPOIN;
                                    M12M0->at(TNPOIN) = icount;}
   }
  }
}

//--------------------------------------------------------------------------//

void WriteTetgen::writeMeshVariables()
{
  for(unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++) {
   if(nodcod->at(IPOIN)>=0) {
    fprintf(file,"%u %s",M02M1->at(IPOIN+1),"  "); // c++ indeces start from 0
    for(unsigned IA=0; IA<PhysicsInfo::getnbDim() ; IA++) {
     fprintf(file,"%.16f %s",(*XYZ)(IA,IPOIN),"  ");}
    for(unsigned IA=0; IA<(*ndof); IA++) {
     fprintf(file,"%.16f %s",(*Zroe)(IA,IPOIN),"  ");}
    fprintf(file,"%u %s",nodcod->at(IPOIN),"\n");
   }
  }
}

//--------------------------------------------------------------------------//

void WriteTetgen::writeUpstreamStatus()
{
  for(unsigned ISH=0; ISH<PhysicsInfo::getnbShMax(); ISH++) {
//   for(unsigned I=0; I<nShockPoints->at(0); I++) {
   for(unsigned I=0; I<nShockPoints->at(0); I++) {

    ++icount;
    if((*NodCodSh)(I,ISH)==10) {
     fprintf(file,"%u %s",M02M1->at(icount),"  ");
     for(unsigned IA=0; IA<PhysicsInfo::getnbDim() ; IA++) {
     fprintf(file,"%0.16f %s",(*XYZShu)(IA,I,ISH),"  ");}
     for(unsigned IA=0; IA<(*ndof); IA++) {
     fprintf(file,"%0.16f %s",(*ZRoeShu)(IA,I,ISH)," ");}
     fprintf(file,"%u %s",(*NodCodSh)(I,ISH),"\n");
    }
   }
  }
}

//--------------------------------------------------------------------------//

void WriteTetgen::writeDownstreamStatus()
{
  for(unsigned ISH=0; ISH<PhysicsInfo::getnbShMax(); ISH++) {
//   for(unsigned I=0; I<nShockPoints->at(0); I++) {
   for(unsigned I=0; I<nShockPoints->at(0); I++) {

     icount++;
    if((*NodCodSh)(I,ISH)==10) {
     fprintf(file,"%u %s",M02M1->at(icount),"  ");
     for(unsigned IA=0; IA<PhysicsInfo::getnbDim() ; IA++) {
     fprintf(file, "%0.16f %s",(*XYZShd)(IA,I,ISH),"  ");}
     for(unsigned IA=0; IA<(*ndof); IA++) {
     fprintf(file,"%0.16f %s",(*ZRoeShd)(IA,I,ISH),"  ");
   if ((*ZRoeShd)(IA,I,ISH)!=(*ZRoeShd)(IA,I,ISH)){
     std::cout << "here tetgen point " << I   <<" nodcod " << (*NodCodSh)(I,ISH) <<'\n';
     std::cout << "(*XYZShd)(IA,I,ISH) " << (*XYZShd)(0,I,ISH)<< " " << (*XYZShd)(1,I,ISH)<< " ";
     std::cout << (*XYZShd)(2,I,ISH)<< " "<<'\n';
  }
     }
     fprintf(file,"%u %s",(*NodCodSh)(I,ISH),"\n");

    }
   }
  }
}

//--------------------------------------------------------------------------//

void WriteTetgen::writeBndfac()
{
  ///Daje Vincent che sei quasi alla fine
  int IBC;


  ICHECK = 0;
  for(unsigned IFACE=0; IFACE<(*nbfacSh); IFACE++) {
   IBC = (*bndfac)(3,IFACE);    // era 2, ora in 3D é 3 perché ci sta una coordinata in piu
   if(IBC>0) { ++ICHECK; }
  }
  fprintf(file,"%u %s",ICHECK,"  1\n");
  
  ICHECK=0;
  

  unsigned m_nbfac=0;
  for(unsigned IFACE=0; IFACE<(*nbfacSh); IFACE++) {
   IBC = (*bndfac)(3,IFACE);

   if(IBC>0) {
    ++ICHECK;
    fprintf(file,"%5u %5u %3d %s",1,0,IBC," \n");
    fprintf(file,"%5u %s",3, " ");
    for(unsigned IA=0; IA<3; IA++) {
     fprintf(file,"%u %s", M02M1->at((*bndfac)(IA,IFACE)),"  ");
    }
    fprintf(file,"%s","\n");
   }
  }

 fprintf(file,"%s","0\n"); // write number of holes

  // for(unsigned IFACE=0;IFACE<NCLR;IFACE++)  {
  //  for(unsigned IB=m_nbfac;IB<nFacB.at(IFACE)+m_nbfac;IB++) {
  //   NBND.str(string());
  //   NBND << namebnd.at(IFACE); // c++ indeces start from 0
  //   if(NBND.str() == "InnerSup" || NBND.str() == "InnerSub") { NBND << "10"; }
  //   fprintf(Tetgenfile,"%5u",3);   // fprintf(Tetgenfile,"%5u",IB+1);
  //   fprintf(Tetgenfile,"%5i %5i %5i %s",bndfac.at(IB*4+0),bndfac.at(IB*4+1),bndfac.at(IB*4+2)," \n"); // DAAGGIUNGERE
  //
  //  }
  //  m_nbfac=m_nbfac+nFacB.at(IFACE);
  // }
  // fprintf(Tetgenfile,"%s","0\n"); // write number of holes
  // fclose(Tetgenfile);

}

//--------------------------------------------------------------------------//

void WriteTetgen::computenbHoles()
{
  nHoles = 0;
  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   nHoles = nHoles + nShockPoints->at(ISH)-2;
  }

  fprintf(file, "%u %s",nHoles + MeshData::getInstance().getnbAddHoles(),"\n");
  unsigned iHole=0;
  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {
   for(unsigned I=1; I<nShockPoints->at(ISH)-1; I++) {
    ++iHole;
    fprintf(file,"%u %s",iHole,"  ");
    fprintf(file,"%.16f %s %.16f %s",(*XYZSh)(0,I,ISH),"  ",(*XYZSh)(1,I,ISH),"\n");
   }
  }

  for (unsigned I=0; I<MeshData::getInstance().getnbAddHoles(); I++) {
   fprintf(file,"%u %s",iHole+I,"  ");
   for(unsigned j=0; j<2; j++) { fprintf(file,"%.16f %s",caddholes->at(j)," ");}
  }
}

//--------------------------------------------------------------------------//

void WriteTetgen::setAddress()
{
  unsigned start; unsigned totsize;
  start = 0;
  XYZ = new Array2D <double> (PhysicsInfo::getnbDim(), npoin->at(0),
                             &coorVect->at(start));
  Zroe = new Array2D <double> (PhysicsInfo::getnbDofMax(),npoin->at(0),
                               &zroeVect->at(start));
  // totsize = nbfac->at(0) + 2 *
  //                          PhysicsInfo::getnbShMax() *
  //                          PhysicsInfo::getnbShEdgesMax();
  // bndfac = new Array2D<int> (4,totsize,&bndfacVect->at(start));
  totsize = nbfac->at(0) + 2 *
                           PhysicsInfo::getnbShMax() *
                           nShockFaces->at(0);
  bndfac = new Array2D<int> (4,totsize,&bndfacVect->at(start));
  celnod = new Array2D<int> ((*nvt), nelem->at(0), &celnodVect->at(start));
  start = npoin->at(0);
  NodCodSh = new Array2D <int> (nShockPoints->at(0),
                                PhysicsInfo::getnbShMax(),
                                &nodcod->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDofMax();
  ZRoeShu = new Array3D <double> (PhysicsInfo::getnbDofMax(),
                                  nShockPoints->at(0),
                                  PhysicsInfo::getnbShMax(),
                                  &zroeVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDofMax() +
          PhysicsInfo::getnbShPointsMax() * PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDofMax();
  ZRoeShd = new Array3D <double> (PhysicsInfo::getnbDofMax(),
                                  nShockPoints->at(0),
                                  PhysicsInfo::getnbShMax(),
                                  &zroeVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDim() ;
  XYZShu = new Array3D <double> (PhysicsInfo::getnbDim() ,
                                nShockPoints->at(0),
                                PhysicsInfo::getnbShMax(),
                                &coorVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDim()  +
          nShockPoints->at(0) * PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDim() ;
  XYZShd = new Array3D <double> (PhysicsInfo::getnbDim() ,
                                nShockPoints->at(0),
                                PhysicsInfo::getnbShMax(),
                                &coorVect->at(start));
}

//--------------------------------------------------------------------------//

void WriteTetgen::freeArray()
{
  delete XYZ; delete Zroe;
  delete bndfac; delete celnod;
  delete NodCodSh;
  delete ZRoeShu; delete ZRoeShd; delete XYZShu; delete XYZShd;
}

//--------------------------------------------------------------------------//

void WriteTetgen::setMeshData()
{
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  nelem =MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nbfac = MeshData::getInstance().getData <vector<unsigned> > ("NBFAC");
  nbfacSh = MeshData::getInstance().getData<unsigned>("NBFACSH");
  caddholes =
    MeshData::getInstance().getData <vector<double> > ("CADDholes");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  bndfacVect = MeshData::getInstance().getData <vector<int> >("BNDFAC");
  celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
  fname = MeshData::getInstance().getData <stringstream> ("FNAME");
  M12M0 = MeshData::getInstance().getData <vector <int> > ("M12M0");
  M02M1 = MeshData::getInstance().getData <vector <unsigned> > ("M02M1");
}

//--------------------------------------------------------------------------//

void WriteTetgen::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  nShockEdges =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  XYZSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYZSh");
  ShFaces =
     PhysicsData::getInstance().getData <Array3D<unsigned> > ("ShFaces");
  nShockFaces =
     PhysicsData::getInstance().getData <vector<unsigned> > ("nShockFaces");
  ShEdgePoints = PhysicsData::getInstance().getData <Array3D <unsigned> > ("ShEDGEPoints");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
