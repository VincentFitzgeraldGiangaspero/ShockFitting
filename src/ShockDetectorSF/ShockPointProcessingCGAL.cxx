#include "ShockDetectorSF/ShockPointProcessingCGAL.hh"
#include "MathTools/MinMax.hh"
#include "Framework/ChemicalInfo.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"

#include <iostream>
#include <vector>
#include <fstream>
#include <utility> // defines std::pair
#include <list>

#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/tuple.h>
#include <boost/lexical_cast.hpp>
#include <algorithm>


#include <CGAL/Kernel_traits.h>
#include <CGAL/Origin.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/edge_aware_upsample_point_set.h>
#include <CGAL/property_map.h>
#include <CGAL/bilateral_smooth_point_set.h>
#include <CGAL/tags.h>

#include "SConfig/ObjectProvider.hh"


//--------------------------------------------------------------------------//

//types
// Types Definition
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel_pred; //from normals
typedef CGAL::Simple_cartesian<double> Kernel; //from jet and wlop
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::cpp11::array<std::size_t,3> Facet;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

//--------------------------------------------------------------------------//
namespace std {
  std::ostream&
  operator<<(std::ostream& os, const Facet& f)
  {
    os << "3 " << f[0] << " " << f[1] << " " << f[2];
    return os;
  }
}
//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//


namespace ShockFitting {

//--------------------------------------------------------------------------//

/// Constructor
ShockPointProcessingCGAL::ShockPointProcessingCGAL(){
}

//--------------------------------------------------------------------------//

/// Destructor
ShockPointProcessingCGAL::~ShockPointProcessingCGAL(){
}

//--------------------------------------------------------------------------//

void ShockPointProcessingCGAL::Working(){

  std::cout << "     => SurfaceReconstruction::Working() \n";
}

//--------------------------------------------------------------------------//

void ShockPointProcessingCGAL::PointReconstruction(){

  std::cout << "     => ShockPointProcessingCGAL::PointReconstruction() \n";

  setMeshData();
  setPhysicsData();
  setAddress();

  // check DOMINIUM size
  vector<double> Xdomain(npoin->at(0));
  vector<double> Ydomain(npoin->at(0));
  vector<double> Zdomain(npoin->at(0));
  for(unsigned IPOIN=0;IPOIN<npoin->at(0);IPOIN++) {
    Xdomain.at(IPOIN) = (*XYZ)(0,IPOIN);
    Ydomain.at(IPOIN) = (*XYZ)(1,IPOIN);
    Zdomain.at(IPOIN) = (*XYZ)(2,IPOIN);
  }

  minDomainXCoordinate = *min_element(Xdomain.begin(),Xdomain.end());
  maxDomainXCoordinate = *max_element(Xdomain.begin(),Xdomain.end());
  minDomainYCoordinate = *min_element(Ydomain.begin(),Ydomain.end());
  maxDomainYCoordinate = *max_element(Ydomain.begin(),Ydomain.end());
  minDomainZCoordinate = *min_element(Zdomain.begin(),Zdomain.end());
  maxDomainZCoordinate = *max_element(Zdomain.begin(),Zdomain.end());


  cout << "-Domain boundaries X " << minDomainXCoordinate << " " << maxDomainXCoordinate << "\n";
  cout << "-Domain boundaries Y " << minDomainYCoordinate << " " << maxDomainYCoordinate << "\n";
  cout << "-Domain boundaries Z " << minDomainZCoordinate << " " << maxDomainZCoordinate << "\n";


  unsigned ISH=0;
  unsigned m_nbShPoints=nShockPoints->at(ISH);

  for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
    points_shock.push_back(Point((*XYZSh)(0,ISHPOIN,ISH),(*XYZSh)(1,ISHPOIN,ISH),(*XYZSh)(2,ISHPOIN,ISH)));
  }


  std::ofstream out_gnoffo("log/shock_points.xyz");
  CGAL::write_xyz_points(out_gnoffo, points_shock.begin(), points_shock.end());

  AutoReconstruction();


  FILE* ShockReconPointsCGAL;
  ShockReconPointsCGAL = fopen("log/ShockReconPointsCGAL.dat","w");

  fprintf(ShockReconPointsCGAL,"%s","TITLE = ShockReconPointsCGAL Points\n");
  fprintf(ShockReconPointsCGAL,"%s","VARIABLES = \"x0\" \"x1\" \"x2\"\n");
  fprintf(ShockReconPointsCGAL,"%s","ZONE T = \"ShockReconPointsCGAL points\"\n");
  fprintf(ShockReconPointsCGAL,"%s","STRANDID=0, SOLUTIONTIME=0 ");
  fprintf(ShockReconPointsCGAL,"%s %u %s","I=",m_nbShPoints,", J=1, K=1, ZONETYPE=Ordered ");
  fprintf(ShockReconPointsCGAL,"%s","DATAPACKING=POINT\n");
  fprintf(ShockReconPointsCGAL,"%s","DT = (SINGLE, SINGLE,SINGLE)\n");

  for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
    fprintf(ShockReconPointsCGAL,"%22.14f",(*XYZSh)(0,ISHPOIN,ISH));
    fprintf(ShockReconPointsCGAL,"%22.14f",(*XYZSh)(1,ISHPOIN,ISH));
    fprintf(ShockReconPointsCGAL,"%22.14f",(*XYZSh)(2,ISHPOIN,ISH));
    fprintf(ShockReconPointsCGAL,"%s","\n");
  }

  fclose(ShockReconPointsCGAL);


  freeArray();
}

//--------------------------------------------------------------------------//

void ShockPointProcessingCGAL::AutoReconstruction()
{
  unsigned ISH=0;
  std::cout << "     => ShockPointProcessingCGAL::AutoReconstruction() \n";
  const char* input_filename = "log/shock_points.xyz";

  // Read Input File as a VECTOR
  std::vector<Point> points_shock;
  std::vector<Point> points_wlop;
  std::ifstream stream(input_filename);
  if (!stream || !CGAL::read_xyz_points(stream, std::back_inserter(points_shock)))
  {
    std::cerr << "Error: cannot read file " << input_filename  << std::endl;
  exit(1);
  }


  //WLOP
  //parameters
  const double retain_percentage = 5;   // percentage of points to retain.
  const double neighbor_radius = 0.5;   // neighbors size.
  bool  	   require_uniform_sampling=true; // AGGIUNTO.
  unsigned int 	iter_number_wlop = 200;

  CGAL::wlop_simplify_and_regularize_point_set
  <Concurrency_tag>
  (points_shock.begin(),
    points_shock.end(),
    std::back_inserter(points_wlop),
    retain_percentage,
    neighbor_radius,
    iter_number_wlop,
    require_uniform_sampling
    );

  std::ofstream out_wlop("log/shock_points_w.xyz");
  if (!out_wlop || !CGAL::write_xyz_points(
    out_wlop, points_wlop.begin(), points_wlop.end()))
  {
  exit(1);
  }

  std::cout << "     => ShockPointProcessingCGAL::WLOP() \n";


  //NORMALS
  std::list<PointVectorPair> points;
  std::ifstream input_norm("log/shock_points_w.xyz");

  if (!input_norm ||
    !CGAL::read_xyz_points(input_norm,
      std::back_inserter(points),
      CGAL::First_of_pair_property_map<PointVectorPair>()))
  {
    std::cerr << "Error: cannot read file " << input_norm << std::endl;
  exit(1);
  }
  // Estimates normals direction.
  // Note: pca_estimate_normals() requires an iterator over points
  // as well as property maps to access each point's position and normal.
  const int nb_neighbors = 10; // K-nearest neighbors = 3 rings
  CGAL::pca_estimate_normals<Concurrency_tag>(points.begin(), points.end(),
    CGAL::First_of_pair_property_map<PointVectorPair>(),
    CGAL::Second_of_pair_property_map<PointVectorPair>(),
    nb_neighbors);
  // Orients normals.
  // Note: mst_orient_normals() requires an iterator over points
  // as well as property maps to access each point's position and normal.
  std::list<PointVectorPair>::iterator unoriented_points_begin =
  CGAL::mst_orient_normals(points.begin(), points.end(),
    CGAL::First_of_pair_property_map<PointVectorPair>(),
    CGAL::Second_of_pair_property_map<PointVectorPair>(),
    nb_neighbors);
  // Optional: delete points with an unoriented normal
  // if you plan to call a reconstruction algorithm that expects oriented normals.
  points.erase(unoriented_points_begin, points.end());

  std::ofstream out_normals("log/shock_points_n.xyz");
  if (!out_normals || !CGAL::write_xyz_points_and_normals(
    out_normals, points.begin(), points.end(),
    CGAL::First_of_pair_property_map<PointVectorPair>(),
    CGAL::Second_of_pair_property_map<PointVectorPair>()
    ))
  {
  exit(1);
  }

  std::cout << "     => ShockPointProcessingCGAL::Normals() \n";

  // SMOOTHING
  //Algorithm parameters
  int k = 10;                 // size of neighborhood. The bigger the smoother the result will be.
  // This value should bigger than 1.
  double sharpness_angle = 90; // control sharpness of the result.
  // The bigger the smoother the result will be
  int iter_number = 1;         // number of times the projection is applied

  for (int i = 0; i < iter_number; ++i){
    CGAL::bilateral_smooth_point_set <Concurrency_tag>(
      points.begin(),
      points.end(),
      CGAL::First_of_pair_property_map<PointVectorPair>(),
      CGAL::Second_of_pair_property_map<PointVectorPair>(),
      k,
      sharpness_angle);
  }

  // Save point set.
  std::ofstream out_jet("log/shock_points_w_s.xyz");
  if (!out_jet || !CGAL::write_xyz_points(out_jet, points.begin(), points.end(),
      CGAL::First_of_pair_property_map<PointVectorPair>()))
  {
    exit(1);
  }

  std::cout << "     => ShockPointProcessingCGAL::Smoothing() \n";


  std::ifstream in("log/shock_points_w_s.xyz");
  // 	std::ifstream in((argc>1)?argv[1]:"data/shock_points_u.xyz");
  double per = 0;
  std::vector<Point> points_mesh;


  std::vector<Facet> facets;
  std::copy(std::istream_iterator<Point>(in),
    std::istream_iterator<Point>(),
    std::back_inserter(points_mesh));

  minDomainXCoordinate= minDomainXCoordinate+0.1;
  maxDomainXCoordinate= maxDomainXCoordinate-0.1;

  cout << "-domain boundaries X " << minDomainXCoordinate << " " << maxDomainXCoordinate << "\n";
  cout << "-domain boundaries Y " << minDomainYCoordinate << " " << maxDomainYCoordinate << "\n";
  cout << "-domain boundaries Z " << minDomainZCoordinate << " " << maxDomainZCoordinate << "\n";


  unsigned IPOIN=0;
  for (std::vector<Point>::iterator unwanted_points_begin=points_mesh.begin();
    unwanted_points_begin != points_mesh.end();unwanted_points_begin++){
    if(points_mesh[IPOIN][0]<=minDomainXCoordinate | points_mesh[IPOIN][0]>=maxDomainXCoordinate |
      points_mesh[IPOIN][1]<=minDomainYCoordinate | points_mesh[IPOIN][1]>=maxDomainYCoordinate |
      points_mesh[IPOIN][2]<=minDomainZCoordinate | points_mesh[IPOIN][2]>=maxDomainZCoordinate){
        //cout << "shock point to be removed " << IPOIN << " with coordinates " << points_mesh[IPOIN] << "\n";
        points_mesh.erase(unwanted_points_begin);
    }
    IPOIN++;
  }

  std::cout << "     => ShockPointProcessingCGAL::RemovingPoints : "<< IPOIN << " \n";


  Perimeter perimeter(per);
  CGAL::advancing_front_surface_reconstruction(points_mesh.begin(),
    points_mesh.end(),
    std::back_inserter(facets),
    perimeter);

  std::cout << "     => ShockPointProcessingCGAL::SurfaceReconstruction() \n";


    // .off file of surface meshed
  std::ofstream out_surf("shock.off");
  out_surf << "OFF\n" << points_mesh.size() << " " << facets.size() << " 0\n";
  std::copy(points_mesh.begin(),
    points_mesh.end(),
    std::ostream_iterator<Point>(out_surf, "\n"));
  std::copy(facets.begin(),
    facets.end(),
    std::ostream_iterator<Facet>(out_surf, "\n"));

  nShockPoints->at(0)=points_mesh.size();
  nShockFaces->at(0) = facets.size();

  ShFaces->resize(3,facets.size(),PhysicsInfo::getnbShMax());

  for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
    (*XYZSh)(0,ISHPOIN,ISH)=points_mesh[ISHPOIN][0];
    (*XYZSh)(1,ISHPOIN,ISH)=points_mesh[ISHPOIN][1];
    (*XYZSh)(2,ISHPOIN,ISH)=points_mesh[ISHPOIN][2];
  }


  for( unsigned IFACE=0; IFACE<nShockFaces->at(0); IFACE++){
    (*ShFaces)(0,IFACE,ISH)=facets[IFACE][0];
    (*ShFaces)(1,IFACE,ISH)=facets[IFACE][1];
    (*ShFaces)(2,IFACE,ISH)=facets[IFACE][2];
  }

  cout << "-nShockFaces:  "<< nShockFaces->at(0) << "\n";
  cout << "-nShockPoints: "<< nShockPoints->at(0) << "\n";
}

//--------------------------------------------------------------------------//

void ShockPointProcessingCGAL::Wlop(){
}

//--------------------------------------------------------------------------//

void ShockPointProcessingCGAL::Normals(){
}

//--------------------------------------------------------------------------//

void ShockPointProcessingCGAL::Smoothing(){
}

//--------------------------------------------------------------------------//

void ShockPointProcessingCGAL::setAddress()
{
  unsigned totsize = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
  PhysicsInfo::getnbShPointsMax();
  XYZ = new Array2D<double> (PhysicsInfo::getnbDim(),
    totsize,
    &coorVect->at(0));

}

//--------------------------------------------------------------------------//

void ShockPointProcessingCGAL::freeArray()
{
  delete XYZ;
}

//--------------------------------------------------------------------------//

void ShockPointProcessingCGAL::setMeshData()
{
  nvt = MeshData::getInstance().getData<unsigned>("NVT");
  nelem = MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
            // zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
            // celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
}

//--------------------------------------------------------------------------//

void ShockPointProcessingCGAL::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  nShockFaces =PhysicsData::getInstance().getData <vector <unsigned> > ("nShockFaces");		
  ShFaces =PhysicsData::getInstance().getData <Array3D<unsigned> > ("ShFaces");
  XYZSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYZSH");

}
//--------------------------------------------------------------------------//

}

