// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ShockDetectorSF/ShockSurfaceRemeshingMMG.hh"
#include "Framework/ChemicalInfo.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"



//--------------------------------------------------------------------------//


namespace ShockFitting {

//--------------------------------------------------------------------------//

/// Constructor
ShockSurfaceRemeshingMMG::ShockSurfaceRemeshingMMG(){
}

//--------------------------------------------------------------------------//

/// Destructor
ShockSurfaceRemeshingMMG::~ShockSurfaceRemeshingMMG(){
}

//--------------------------------------------------------------------------//

void ShockSurfaceRemeshingMMG::Remeshing(){
  
  LogToScreen(INFO, "ShockSurfaceRemeshingMMG::Remeshing()\n");

  setPhysicsData();

  // createMMGfile();

  command3 = "../../ShockDetectorSF/MMG/meshconv shock_smooth.off -c stl > log/MeshconvExe.log";
  command2 = "../../ShockDetectorSF/MMG/gmsh -3 shock_smooth.stl -format mesh >log/GmshExe.log";
  command = "../../ShockDetectorSF/MMG/mmgs_O3 shock_smooth.mesh -hmax 0.08 > log/MMGExe.log";

  system(command3.c_str());
  if(system(command3.c_str())!=0) {
    std::cout << "meshconv::error => surface remeshing execution failed\n";
    exit(1); 
  }

  system(command2.c_str());
  if(system(command2.c_str())!=0) {
   std::cout << "Gmsh::error => surface remeshing execution failed\n";
   exit(1); }


  system(command.c_str());
  if(system(command.c_str())!=0) {
    std::cout << "MMG::error => surface remeshing execution failed\n";
    exit(1); 
  }

  readMMGfile();

  detectEdges();

}

//--------------------------------------------------------------------------//

void ShockSurfaceRemeshingMMG::createMMGfile()
{

  std::cout << "     => ShockSurfaceRemeshingMMG::createMMGfile()\n";


  FILE* meshfile;
  meshfile = fopen("shock_smooth.mesh","w");

  fprintf(meshfile,"%s","MeshVersionFormatted 2\n");
  fprintf(meshfile,"%s","Dimension\n");
  fprintf(meshfile,"%s" ,"3\n");
  fprintf(meshfile,"%s" ,"Vertices\n");
  fprintf(meshfile,"%u %s" ,nShockPoints->at(0),"\n");

  for(unsigned ISH=0; ISH<(*nShocks); ISH++){
    for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
      fprintf(meshfile,"%f %s" ,(*XYZSh)(0,ISHPOIN,ISH)," ");
      fprintf(meshfile,"%f %s" ,(*XYZSh)(1,ISHPOIN,ISH)," ");
      fprintf(meshfile,"%f %s" ,(*XYZSh)(2,ISHPOIN,ISH)," ");
      fprintf(meshfile,"%s", "1\n");
    }
    fprintf(meshfile,"%s","Triangles\n");
    fprintf(meshfile,"%u %s",nShockFaces->at(ISH),"\n");
    for(unsigned IFACE=0;IFACE<nShockFaces->at(ISH);IFACE++) {
      fprintf(meshfile,"%u %s" ,(*ShFaces)(0,IFACE,ISH)+1," ");
      fprintf(meshfile,"%u %s" ,(*ShFaces)(1,IFACE,ISH)+1," ");
      fprintf(meshfile,"%u %s" ,(*ShFaces)(2,IFACE,ISH)+1," ");
      fprintf(meshfile,"%s", "1\n");
    }
  }

  fprintf(meshfile,"%s", " END\n");

  std::cout << "     => ShockSurfaceRemeshingMMG::ENDcreateMMGfile()\n";

}

//--------------------------------------------------------------------------//

void ShockSurfaceRemeshingMMG::readMMGfile()
{
  std::cout << "     => ShockSurfaceRemeshingMMG::readMMGfile()\n";

  std::string dummyfile;
  std::string dummy;

  // reading file
  std::ifstream file;

  dummyfile = "shock_smooth.o.mesh";

  file.open(dummyfile.c_str()); 

  while(dummy != "Vertices"){
    file >> dummy;
  }

  file >> nShockPoints->at(0);
  std::cout << "-nShockPoints: " << nShockPoints->at(0) << '\n';


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


  for (unsigned IFACE=0; IFACE<nShockFaces->at(0); IFACE++) {
    for (unsigned INODE=0; INODE<3; INODE++){
      file >> (*ShFaces)(INODE,IFACE,0);            // node
      (*ShFaces)(INODE,IFACE,0)=(*ShFaces)(INODE,IFACE,0)-1;
    }
    file >> dummy; //1
  }

  file >>dummy;

}

//--------------------------------------------------------------------------//


void ShockSurfaceRemeshingMMG::detectEdges(){

  std::cout << "     => ShockSurfaceRemeshingMMG::detectEdges() \n";

  unsigned ISH=0;

  ShEdgePoints->resize(1,PhysicsInfo::getnbShPointsMax(),PhysicsInfo::getnbShMax());

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
      (*ShEdgePoints)(0,NPOIN,0)=1;
    }
  }

  std::vector<std::vector<double> > edge_points(3,std::vector <double> (npoinedge));
  unsigned NEDGE=0;

  for (unsigned NPOIN=0; NPOIN < nShockPoints->at(0); NPOIN++){
    if((*ShEdgePoints)(0,NPOIN,0)==1){
      edge_points[0][NEDGE]=(*XYZSh)(0,NPOIN,ISH);
      edge_points[1][NEDGE]=(*XYZSh)(1,NPOIN,ISH);
      edge_points[2][NEDGE]=(*XYZSh)(2,NPOIN,ISH);
      NEDGE++;
    }
  }

  std::cout << "-nShockEdgePoints: " << NEDGE << "\n";

  // plot shock edges
  FILE* shockEdges;
  shockEdges = fopen("ShockEdges.dat","w");

  fprintf(shockEdges,"%s","TITLE = Detected shock edge points\n");
  fprintf(shockEdges,"%s","VARIABLES = \"x0\" \"x1\" \"x2\" \n");
  fprintf(shockEdges,"%s","ZONE T = \"Detected ShockEdges\"\n");
  fprintf(shockEdges,"%s","STRANDID=0, SOLUTIONTIME=0 ");
  fprintf(shockEdges,"%s %u %s","I=",NEDGE,", J=1, K=1, ZONETYPE=Ordered ");
  fprintf(shockEdges,"%s","DATAPACKING=POINT\n");
  fprintf(shockEdges,"%s","DT = (SINGLE, SINGLE, SINGLE)\n");

  for(unsigned ISHPOIN=0;ISHPOIN<NEDGE;ISHPOIN++) {
   fprintf(shockEdges,"%22.14F",edge_points[0][ISHPOIN]);
   fprintf(shockEdges,"%22.14F",edge_points[1][ISHPOIN]);
   fprintf(shockEdges,"%22.14F",edge_points[2][ISHPOIN]);		
   fprintf(shockEdges,"%s","\n");
  }  

  fclose(shockEdges);


}

//--------------------------------------------------------------------------//

void ShockSurfaceRemeshingMMG::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints = PhysicsData::getInstance().getData <std::vector <unsigned> > ("nShockPoints");
  nShockFaces = PhysicsData::getInstance().getData <std::vector <unsigned> > ("nShockFaces");		
  ShEdgePoints = PhysicsData::getInstance().getData <Array3D <unsigned> > ("ShEDGEPoints");
  XYZSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYZSH");
  ShFaces = PhysicsData::getInstance().getData <Array3D<unsigned> > ("ShFaces");
}

//--------------------------------------------------------------------------//


}//namespace ShockFitting
