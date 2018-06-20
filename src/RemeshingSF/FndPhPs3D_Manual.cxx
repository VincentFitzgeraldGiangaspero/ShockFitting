// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RemeshingSF/FndPhPs3D_Manual.hh"
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
#include "MathTools/DotProd.hh"
#include "MathTools/CrossProd.hh"



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

// vincent 03/18

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<FndPhPs3D_Manual, Remeshing> findPhPs3D_manualProv("FndPhPs3D_Manual");

//--------------------------------------------------------------------------//

FndPhPs3D_Manual::FndPhPs3D_Manual(const std::string& objectName) :
  Remeshing(objectName)
{
}

//--------------------------------------------------------------------------//

FndPhPs3D_Manual::~FndPhPs3D_Manual()
{
}

//--------------------------------------------------------------------------//

void FndPhPs3D_Manual::setup()
{
  LogToScreen(VERBOSE, "FndPhPs3D_Manual::setup() => start\n");

  LogToScreen(VERBOSE, "FndPhPs3D_Manual::setup() => end\n");
}

//--------------------------------------------------------------------------//

void FndPhPs3D_Manual::unsetup()
{
  LogToScreen(VERBOSE, "FndPhPs3D_Manual::unsetup()\n");
}

//--------------------------------------------------------------------------//

void FndPhPs3D_Manual::remesh()
{
  LogToScreen(INFO, "FndPhPs3D_Manual::remesh()\n");

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



  double Area;
  double distmin,distap;
  std::vector<double> dst;
  std::vector<double> counter;
  std::vector<double> diff(3);
  std::vector<double> b(3);
  std::vector<double> vertex_a(3);
  std::vector<double> vertex_b(3);
  std::vector<double> vertex_3(3);
  std::vector<double> ab(3),ac(3);


  unsigned verta,vertb,vertc;
  dst.resize(nShockPoints->at(0));
  counter.resize(nShockPoints->at(0));

  CrossProd <double> Cross;
  DotProd <double> Dot;

  double shktol=0.65;
  double pi;
  pi=sqrt(3.0)*shktol;

  int start_s=clock();
  unsigned NVT;
  struct timeval start, end;
  gettimeofday(&start, NULL);

   for(unsigned iElemSh=0; iElemSh < nShockFaces->at(0); iElemSh++) {

	   verta=(*ShFaces)(0,iElemSh,0);
	   vertb=(*ShFaces)(1,iElemSh,0);
	   vertc=(*ShFaces)(2,iElemSh,0);

	   for (unsigned IV=0;IV<3;IV++){
		   ab[IV]=(*XYZSh)(IV,vertb,0)-(*XYZSh)(IV,verta,0);
		   ac[IV]=(*XYZSh)(IV,vertc,0)-(*XYZSh)(IV,verta,0);

	   }
	   Cross.computeCrossProd(ab,ac);
	   diff=Cross.getCrossProd();
	   Dot.computeDotProd(diff,diff);
	   Area=0.5*sqrt(Dot.getDotProd());

	   counter[verta]=counter[verta]+1;
	   counter[vertb]=counter[vertb]+1;
	   counter[vertc]=counter[vertc]+1;

       dst[verta]=dst[verta]+Area;
       dst[vertb]=dst[vertb]+Area;
       dst[vertc]=dst[vertc]+Area;

   }


   for (unsigned IPOIN=0; IPOIN<nShockPoints->at(0); IPOIN++){
	   dst[IPOIN]=dst[IPOIN]*1.0/counter[IPOIN];
	   //std::cout << dst[IPOIN] <<endl;
	   dst[IPOIN]=dst[IPOIN]*pi;
   }

   std::cout << "-fine fase shock " <<endl;

   for (unsigned IPOIN=0; IPOIN<npoin->at(0); IPOIN++){
	   distmin=1.0;

	   for (unsigned ISHPOIN=0; ISHPOIN<nShockPoints->at(0); ISHPOIN++){

		   for (unsigned IV=0; IV<3;IV++){
			   b[IV]=(*XYZSh)(IV,ISHPOIN,0)-(*XYZ)(IV,IPOIN);
		   }

		   Dot.computeDotProd(b,b);
		   distap=Dot.getDotProd();
		  // std::cout << "ciao" << endl;

		   if (distap<distmin){
			distmin=distap;
			vertb=ISHPOIN;
		   }
	   }


	   if (distmin<dst[vertb]&&(nodcod->at(IPOIN)==0)){
		   nodcod->at(IPOIN)=-1;
	   }
   }


   std::cout << "-fine fase mesh point " <<endl;


  // execution time
  int stop_s=clock();
  gettimeofday(&end, NULL);
  long delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e4;

  //std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << endl;
  std::cout << "-time OMP: " << delta << endl;



  countPhanPoints();

  // de-allocate dynamic arrays
  freeArray();

  logfile.Close();
}

//--------------------------------------------------------------------------//

bool FndPhPs3D_Manual::cellCrossed(unsigned ISH_index, unsigned ielemsh_index,
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

void FndPhPs3D_Manual::setPhanPoints()
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

void FndPhPs3D_Manual::countPhanPoints()
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

void FndPhPs3D_Manual::setIndex(unsigned ISH_index, unsigned ielemsh_index,
                         unsigned ielem_index)
{
  ISH = ISH_index; ielem = ielem_index; ielemsh = ielemsh_index;
}

//--------------------------------------------------------------------------//

void FndPhPs3D_Manual::setAddress()
{
  unsigned start = 0;
  XYZ = new Array2D <double> (PhysicsInfo::getnbDim(),
                             npoin->at(0),
                             &coorVect->at(start));
  celnod = new Array2D <int> ((*nvt), nelem->at(0), &celnodVect->at(start));
}

//--------------------------------------------------------------------------//

void FndPhPs3D_Manual::freeArray()
{
  delete XYZ; delete celnod;
}

//--------------------------------------------------------------------------//

void FndPhPs3D_Manual::setPhysicsData()
{
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints = PhysicsData::getInstance().getData <vector<unsigned>> ("nShockPoints");
  nShockEdges =
      PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");
  nShockFaces =
      PhysicsData::getInstance().getData <vector <unsigned> > ("nShockFaces");
  XYZSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYZSH");

  ShFaces = PhysicsData::getInstance().getData <Array3D <unsigned> > ("ShFaces");

}

//--------------------------------------------------------------------------//

void FndPhPs3D_Manual::setMeshData()
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
