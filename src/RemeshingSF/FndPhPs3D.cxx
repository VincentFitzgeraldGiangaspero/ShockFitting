// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RemeshingSF/FndPhPs3D.hh"
#include "SConfig/ObjectProvider.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/Remeshing.hh"
#include "MathTools/Array3D.hh"
#include "MathTools/Array2D.hh"
#include "MathTools/Ishel.hh"
#include "RemeshingSF/Rdshp.hh"


#include <CGAL/intersections.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/squared_distance_3.h>
#include <algorithm>
#include <vector>
#include <fstream>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Origin.h>
#include <CGAL/squared_distance_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <omp.h>
#include <sys/time.h>


typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Tetrahedron_3<K> Tetrahedron_3;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;
typedef Polyhedron_3::Vertex_iterator        Vertex_iterator;
typedef CGAL::Triangle_3<K> Triangle_3;

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron_3> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

// vincent 11/17

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<FndPhPs3D, Remeshing> findPhPs3DProv("FndPhPs3D");

//--------------------------------------------------------------------------//

FndPhPs3D::FndPhPs3D(const std::string& objectName) :
  Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

FndPhPs3D::~FndPhPs3D()
{
}

//--------------------------------------------------------------------------//

void FndPhPs3D::setup()
{
  LogToScreen(VERBOSE, "FndPhPs3D::setup() => start\n");

  LogToScreen(VERBOSE, "FndPhPs3D::setup() => end\n");
}

//--------------------------------------------------------------------------//

void FndPhPs3D::unsetup()
{
  LogToScreen(VERBOSE, "FndPhPs3D::unsetup()\n");
}

//--------------------------------------------------------------------------//

void FndPhPs3D::remesh()
{
  LogToScreen(INFO, "FndPhPs3D::remesh()\n");

  setMeshData();
  setPhysicsData();
  logfile.Open(getClassName());

  setAddress();

  int counter_inter=0;
  int counter_close=0;

  xc.resize(3);
  yc.resize(3);
  n.resize(3);
  d.resize(3);

  std::vector<Tetrahedron_3> T;
  std::vector<Polyhedron_3> P;
  std::vector<Triangle_3> A;
  std::vector<double> dist_vect;
  std::vector<unsigned> node_index_vect;
  double dist,dist_max,dist_min;
  dist_min=1;
  dist_max=0;
  unsigned node_index;
  
  std::cout << "-Vectors creation phase: start \n"; 

  // mesh cells creation
  for (unsigned iElem=0; iElem < nelem->at(0); iElem++) {
    Polyhedron_3 Poly;
    Poly.make_tetrahedron(Point_3((*XYZ)(0,(*celnod)(0,iElem)-1),
                        (*XYZ)(1,(*celnod)(0,iElem)-1),
                        (*XYZ)(2,(*celnod)(0,iElem)-1)),
                Point_3((*XYZ)(0,(*celnod)(1,iElem)-1),
                        (*XYZ)(1,(*celnod)(1,iElem)-1),
                        (*XYZ)(2,(*celnod)(1,iElem)-1)),
                Point_3((*XYZ)(0,(*celnod)(2,iElem)-1),
                        (*XYZ)(1,(*celnod)(2,iElem)-1),
                        (*XYZ)(2,(*celnod)(2,iElem)-1)),
                Point_3((*XYZ)(0,(*celnod)(3,iElem)-1),
                        (*XYZ)(1,(*celnod)(3,iElem)-1),
                        (*XYZ)(2,(*celnod)(3,iElem)-1)));
    P.push_back(Poly);

//        T.push_back(Tetrahedron_3(Point_3((*XYZ)(0,(*celnod)(0,iElem)-1),
//                            (*XYZ)(1,(*celnod)(0,iElem)-1),
//                            (*XYZ)(2,(*celnod)(0,iElem)-1)),
//                    Point_3((*XYZ)(0,(*celnod)(1,iElem)-1),
//                            (*XYZ)(1,(*celnod)(1,iElem)-1),
//                            (*XYZ)(2,(*celnod)(1,iElem)-1)),
//                    Point_3((*XYZ)(0,(*celnod)(2,iElem)-1),
//                            (*XYZ)(1,(*celnod)(2,iElem)-1),
//                            (*XYZ)(2,(*celnod)(2,iElem)-1)),
//                    Point_3((*XYZ)(0,(*celnod)(3,iElem)-1),
//                            (*XYZ)(1,(*celnod)(3,iElem)-1),
//                            (*XYZ)(2,(*celnod)(3,iElem)-1))));

  }

  // shock surface creation
  for(unsigned iElemSh=0; iElemSh < nShockFaces->at(0); iElemSh++) {


    A.push_back(Triangle_3(Point_3((*XYZSh)(0,(*ShFaces)(0,iElemSh,0),0),
                                    (*XYZSh)(1,(*ShFaces)(0,iElemSh,0),0),
                                    (*XYZSh)(2,(*ShFaces)(0,iElemSh,0),0)),
                            Point_3((*XYZSh)(0,(*ShFaces)(1,iElemSh,0),0),
                                    (*XYZSh)(1,(*ShFaces)(1,iElemSh,0),0),
                                    (*XYZSh)(2,(*ShFaces)(1,iElemSh,0),0)),
                            Point_3((*XYZSh)(0,(*ShFaces)(2,iElemSh,0),0),
                                    (*XYZSh)(1,(*ShFaces)(2,iElemSh,0),0),
                                    (*XYZSh)(2,(*ShFaces)(2,iElemSh,0),0))));

                                  }


  std::cout << "-Vectors creation phase: end \n"; 
  int start_s=clock();
  unsigned NVT;
  struct timeval start, end;
  gettimeofday(&start, NULL);
  
  


  for (unsigned iSh=0; iSh < (*nShocks); iSh++) {
 
   #pragma omp parallel for 
   for(unsigned iElem=0; iElem < nelem->at(0); iElem++) {
     Tree tree(faces(P[iElem]).first, faces(P[iElem]).second, P[iElem]);
     // #omp 
     for (unsigned iElemSh=0; iElemSh < nShockFaces->at(iSh); iElemSh++) {
      if(tree.do_intersect(A[iElemSh])){
        counter_inter++;
        NVT=0;
        for (Vertex_iterator v = P[iElem].vertices_begin(); v != P[iElem].vertices_end(); ++v){
          dist=CGAL::squared_distance(v->point(),A[iElemSh]);
          // dist=CGAL::squared_distance(T[iElem][NVT],A[iElemSh]);
          if (dist<MeshData::getInstance().getSNDMIN()){
            node_index = (*celnod)(NVT,iElem)-1;
            if(dist>0 && dist<(MeshData::getInstance().getSNDMIN())
              && nodcod->at(node_index)==0)
               nodcod->at(node_index)=-1; 
            // if(dist>0 && dist<(MeshData::getInstance().getSNDMIN()) // -2 is for boundary points
            //    && nodcod->at(node_index)>0)
            //   { nodcod->at(node_index)=-2; }
            counter_close++;
              setPhanPoints();
          } //if distance
          NVT++;
        } // for vertices
     } // if intersect
    } // iElem
   } // iElemSh
  } // iSh
  
  

/*
  for (unsigned iSh=0; iSh < (*nShocks); iSh++) {
  #pragma omp parallel for collapse(2)
  // for(unsigned iElemSh=0; iElemSh < nShockFaces->at(iSh); iElemSh++) {
  //  for (unsigned iElem=0; iElem < nelem->at(0); iElem++) {
    
    
   for(unsigned iElemSh=0; iElemSh < nShockFaces->at(iSh); iElemSh++) {
 		for (unsigned iElem=0; iElem < nelem->at(0); iElem++) {
  
     if(do_intersect(T[iElem],A[iElemSh])) {
  
       // counter_inter++;
        for (unsigned NVT=0; NVT<(*nvt);NVT++){
        	dist=CGAL::squared_distance(T[iElem][NVT],A[iElemSh]);
        	if (dist>dist_max){dist_max=dist;}
        	if (dist<dist_min){dist_min=dist;}


        		if (dist<MeshData::getInstance().getSNDMIN()){
        //std::cout << " distance vertex  " << CGAL::squared_distance(T[iElem][NVT],A[iElemSh]) << " \n";
       // dist_vect.push_back(dist);

        			node_index = (*celnod)(NVT,iElem)-1;
       // node_index_vect.push_back(node_index);
        			if(dist>0 && dist<(MeshData::getInstance().getSNDMIN())
        					&& nodcod->at(node_index)==0)
        			{ nodcod->at(node_index)=-1; }
        // if(dist>0 && dist<(MeshData::getInstance().getSNDMIN())
        //    && nodcod->at(node_index)>0)
        //   { nodcod->at(node_index)=-2; }
      //  counter_close++;
        //  setPhanPoints();
        }
        }
      // if distance is too small the node become a phantom node
     } // if
    } // iElem
   } // iElemSh
  } // iSh
*/


  // execution time 
  int stop_s=clock();
  gettimeofday(&end, NULL);      
  long delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e4;

  //std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << endl;
  std::cout << "-time OMP: " << delta << endl;

  std::cout << "-dist_min = " << dist_min <<endl;
  std::cout << "-dist_max = " << dist_max <<endl;


  countPhanPoints();

  // de-allocate dynamic arrays
  freeArray();

  logfile.Close();
}

//--------------------------------------------------------------------------//

bool FndPhPs3D::cellCrossed(unsigned ISH_index, unsigned ielemsh_index,
                          unsigned ielem_index)
{
  setIndex(ISH_index,ielemsh_index,ielem_index);
  for (unsigned i=0; i<(*nvt); i++) {
   n.at(i) = (*celnod)(i,ielem)-1; //c++ indeces start from 0
   xc.at(i) = (*XYZ)(0,n.at(i));
   yc.at(i) = (*XYZ)(1,n.at(i));
  }


  xs1 = (*XYZSh)(0,ielemsh,ISH);
  ys1 = (*XYZSh)(1,ielemsh,ISH);
  xs2 = (*XYZSh)(0,ielemsh+1,ISH);
  ys2 = (*XYZSh)(1,ielemsh+1,ISH);


  Ishel isCellCrossed (xc,yc,xs1,xs2,ys1,ys2);

  if (isCellCrossed.Ishel1()==0) {
   if (isCellCrossed.Ishel2()==0) {
    return true;}
  }
  return false;
}

//--------------------------------------------------------------------------//

void FndPhPs3D::setPhanPoints()
{
  if(d.at(0)>0 && d.at(0)<(MeshData::getInstance().getSNDMIN())
     && nodcod->at(n.at(0))==0)
    { nodcod->at(n.at(0))=-1; }
  if(d.at(0)>0 && d.at(0)<(MeshData::getInstance().getSNDMIN())
     && nodcod->at(n.at(0))>0)
    { nodcod->at(n.at(0))=-2; }
  if(d.at(1)>0 && d.at(1)<(MeshData::getInstance().getSNDMIN())
     && nodcod->at(n.at(1))==0)
    { nodcod->at(n.at(1))=-1; }
  if(d.at(1)>0 && d.at(1)<(MeshData::getInstance().getSNDMIN())
     && nodcod->at(n.at(1))>0)
    { nodcod->at(n.at(1))=-2; }
  if(d.at(2)>0 && d.at(2)<(MeshData::getInstance().getSNDMIN())
     && nodcod->at(n.at(2))==0)
    { nodcod->at(n.at(2))=-1; }
  if(d.at(2)>0 && d.at(2)<(MeshData::getInstance().getSNDMIN())
     && nodcod->at(n.at(2))>0)
    { nodcod->at(n.at(2))=-2; }
}

//--------------------------------------------------------------------------//

void FndPhPs3D::countPhanPoints()
{

  FILE* BndryPhNodes;
  BndryPhNodes = fopen("log/BndryPhNodes.xyz","w");

  *nPhanPoints = 0; *nBoundPhanPoints = 0;
  for (unsigned K=0; K<npoin->at(0); K++) {
  unsigned inod = K+1; // c++ indeces start from 0
   if (nodcod->at(K)==-1 || nodcod->at(K)==-2) {
     *nPhanPoints = *nPhanPoints + 1; }
   if (nodcod->at(K)==-2) {
     *nBoundPhanPoints = *nBoundPhanPoints + 1; }
   if (nodcod->at(K)==-1) { logfile ("Node ", inod, "has become a phantom\n");}
   if (nodcod->at(K)==-2) {
    logfile ("Node ", inod, "on the boundary has become a phantom\n");
    fprintf(BndryPhNodes,"%10f %s",(*XYZ)(0,K)," ");
    fprintf(BndryPhNodes,"%10f %s",(*XYZ)(1,K)," ");
    fprintf(BndryPhNodes,"%10f",(*XYZ)(2,K));
    fprintf(BndryPhNodes,"%s","\n");
}
  }
  logfile ("Number of Phantom nodes (incl. those on the bndry)",
           *nPhanPoints, "\n");
  std::cout << "-nPhPoints: " << *nPhanPoints << "\n";
  if (*nBoundPhanPoints > 0) {
   logfile("Uh! Oh! there are ", *nBoundPhanPoints,
           "phantom nodes on the boundary");
  }

  fclose(BndryPhNodes);
}

//--------------------------------------------------------------------------//

void FndPhPs3D::setIndex(unsigned ISH_index, unsigned ielemsh_index,
                         unsigned ielem_index)
{
  ISH = ISH_index; ielem = ielem_index; ielemsh = ielemsh_index;
}

//--------------------------------------------------------------------------//

void FndPhPs3D::setAddress()
{
  unsigned start = 0;
  XYZ = new Array2D <double> (PhysicsInfo::getnbDim(),
                             npoin->at(0),
                             &coorVect->at(start));
  celnod = new Array2D <int> ((*nvt), nelem->at(0), &celnodVect->at(start));
}

//--------------------------------------------------------------------------//

void FndPhPs3D::freeArray()
{
  delete XYZ; delete celnod;
}

//--------------------------------------------------------------------------//

void FndPhPs3D::setPhysicsData()
{
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockEdges =
      PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  nShockFaces =
      PhysicsData::getInstance().getData <vector <unsigned> > ("nShockFaces");
  XYZSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYZSH");

  ShFaces = PhysicsData::getInstance().getData <Array3D <unsigned> > ("ShFaces");

}

//--------------------------------------------------------------------------//

void FndPhPs3D::setMeshData()
{
  nedge = MeshData::getInstance().getData <vector<unsigned> > ("NEDGE");
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  nelem = MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  nPhanPoints = MeshData::getInstance().getData <unsigned> ("nPhanPoints");
  nBoundPhanPoints =
             MeshData::getInstance().getData <unsigned> ("nBoundPhanPoints");
  nodcod = MeshData::getInstance().getData <vector<int> >("NODCOD");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  celnodVect = MeshData::getInstance().getData < vector<int> >("CELNOD");
}

//--------------------------------------------------------------------------//

} //namespace ShockFitting
