// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "StateUpdaterSF/InterpShockPointStateRH.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"

//--------------------------------------------------------------------------//

//vincent 28/04/2017

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

InterpShockPointStateRH::InterpShockPointStateRH()
{
}

//--------------------------------------------------------------------------//

InterpShockPointStateRH::~InterpShockPointStateRH()
{
}

//--------------------------------------------------------------------------//



//--------------------------------------------------------------------------//

void InterpShockPointStateRH::InterpUpState(unsigned POINT)
{

  setMeshData();
  setPhysicsData();
  setAddress();

  unsigned ISH=0;

      // std::vector<Point_3> shock_points_d(nShockPoints->at(ISH));
      std::vector<Point_3> shock_points_u(nShockPoints->at(ISH));



      // leggi gli shock_points downstream e upstream
      for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
        // shock_points_d.at(ISHPOIN)=Point_3(XYZShd(0,ISHPOIN,ISH),XYZShd(1,ISHPOIN,ISH),XYZShd(2,ISHPOIN,ISH));
        shock_points_u.at(ISHPOIN)=Point_3((*XYZSh)(0,ISHPOIN,ISH),(*XYZSh)(1,ISHPOIN,ISH),(*XYZSh)(2,ISHPOIN,ISH));
      }




      // K-neighborous search from CGAL, to get back the ID of the closest node of the mesh

      My_point_property_map ppmap(shock_points_u);

      // Insert number_of_data_points in the tree
      Tree tree(
        boost::counting_iterator<std::size_t>(0),
        boost::counting_iterator<std::size_t>(shock_points_u.size()),
        Splitter(),
        Traits(ppmap));

      Distance tr_dist(ppmap);
      // number of neighbors points
      const unsigned int K=10;

      std::vector<double> new_up_state(K-1);

      std::vector<double> distance(K);
      std::vector<double> point(K);


      unsigned counter_close=0;

      // working variables and vectors


      FILE* RH_interp_up_state;
      RH_interp_up_state = fopen("log/RH_interp_up_state.dat","a");

      fprintf(RH_interp_up_state,"%s","TITLE = RH_interp_up_state containing shock points\n");

      K_Neighbor_search search(tree, shock_points_u[POINT], K,0,true,tr_dist);

      fprintf(RH_interp_up_state,"%s","------------------------------------------\n");
      fprintf(RH_interp_up_state,"%s %u %s %u %s " ,"The ", K , " nearest upstream shock point to the upstream shock point ",
      POINT+1, " with coordinates ");
      fprintf(RH_interp_up_state, "%22.14f %22.14f %22.14f %s" ,(*XYZSh)(0,POINT,ISH) , (*XYZSh)(1,POINT,ISH) , (*XYZSh)(2,POINT,ISH) , "\n");

      for(K_Neighbor_search::iterator it = search.begin(); it != search.end(); it++){

        fprintf(RH_interp_up_state,"%s %u %s %u %s %f %s %f %s %f %s"," -> ", counter_close , " -near-point ",
        it->first, " : ", (*XYZSh)(0,it->first,ISH), " ", (*XYZSh)(1,it->first,ISH), " ", (*XYZSh)(2,it->first,ISH),"\n");
        fprintf(RH_interp_up_state, "%s %f %s"," at distance " , tr_dist.inverse_of_transformed_distance(it->second) , "\n");
        fprintf(RH_interp_up_state, "%s %f %s %f"," with state " , (*ZRoeShu)(0,it->first,ISH) , " ", (*ZRoeShu)(1,it->first,ISH));
        fprintf(RH_interp_up_state, "%s %f %s %f"," ",(*ZRoeShu)(2,it->first,ISH), " " ,(*ZRoeShu)(3,it->first,ISH));
        fprintf(RH_interp_up_state, "%s %f %s"," ",(*ZRoeShu)(4,it->first,ISH), " \n");

        point.at(counter_close)=it->first;
        distance.at(counter_close)=tr_dist.inverse_of_transformed_distance(it->second);
        counter_close ++;
      }


double dist;

for (unsigned I=0; I< distance.size(); I++)
{
  dist=dist+distance.at(I);
}

for (unsigned IV=0; IV<(*ndof); IV++){
for (unsigned I=0; I<distance.size();I++){
new_up_state.at(IV)=new_up_state.at(IV)+(distance.at(I)/dist)*((*ZRoeShu)(IV,point.at(I),ISH));
}
}


std::cout << " new UP state for point " << POINT << " \n";

// for (unsigned IV=0; IV<(*ndof); IV++){ std::cout << " new_up_state.at(IV) " << new_up_state.at(IV) << " \n";}

for (unsigned IV=0; IV<(*ndof); IV++){
(*ZRoeShu)(IV,POINT,ISH)=new_up_state.at(IV);
}
  // de-allocate dynamic arrays
  freeArray();

}

//--------------------------------------------------------------------------//

void InterpShockPointStateRH::InterpDownState(unsigned POINT)
{

  setMeshData();
  setPhysicsData();
  setAddress();
  unsigned IPOIN= POINT;

  unsigned ISH=0;

      // std::vector<Point_3> shock_points_d(nShockPoints->at(ISH));
      std::vector<Point_3> shock_points_d(nShockPoints->at(ISH));



      // leggi gli shock_points downstream e upstream
      for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
        // shock_points_d.at(ISHPOIN)=Point_3(XYZShd(0,ISHPOIN,ISH),XYZShd(1,ISHPOIN,ISH),XYZShd(2,ISHPOIN,ISH));
        shock_points_d.at(ISHPOIN)=Point_3((*XYZSh)(0,ISHPOIN,ISH),(*XYZSh)(1,ISHPOIN,ISH),(*XYZSh)(2,ISHPOIN,ISH));
      }




      // K-neighborous search from CGAL, to get back the ID of the closest node of the mesh

      My_point_property_map ppmap(shock_points_d);

      // Insert number_of_data_points in the tree
      Tree tree(
        boost::counting_iterator<std::size_t>(0),
        boost::counting_iterator<std::size_t>(shock_points_d.size()),
        Splitter(),
        Traits(ppmap));

      Distance tr_dist(ppmap);
      // number of neighbors points
      const unsigned int K=10;

      std::vector<double> new_down_state(K-1);

      std::vector<double> distance(K);
      std::vector<double> point(K);


      unsigned counter_close=0;

      // working variables and vectors

      std::cout << "here 1" << '\n';

      FILE* RH_interp_down_state;
      RH_interp_down_state = fopen("log/RH_interp_down_state.dat","a");

      fprintf(RH_interp_down_state,"%s","TITLE = RH_interp_down_state containing shock points\n");

      K_Neighbor_search search(tree, shock_points_d[IPOIN], K,0,true,tr_dist);

      fprintf(RH_interp_down_state,"%s","------------------------------------------\n");
      fprintf(RH_interp_down_state,"%s %u %s %u %s " ,"The ", K , " nearest downstream shock point to the upstream shock point ",
      POINT+1, " with coordinates ");
      fprintf(RH_interp_down_state, "%22.14f %22.14f %22.14f %s" ,(*XYZSh)(0,IPOIN,ISH) , (*XYZSh)(1,IPOIN,ISH) , (*XYZSh)(2,IPOIN,ISH) , "\n");

      for(K_Neighbor_search::iterator it = search.begin(); it != search.end(); it++){

        fprintf(RH_interp_down_state,"%s %u %s %u %s %f %s %f %s %f %s"," -> ", counter_close , " -near-point ",
        it->first, " : ", (*XYZSh)(0,it->first,ISH), " ", (*XYZSh)(1,it->first,ISH), " ", (*XYZSh)(2,it->first,ISH),"\n");
        fprintf(RH_interp_down_state, "%s %f %s"," at distance " , tr_dist.inverse_of_transformed_distance(it->second) , "\n");
        fprintf(RH_interp_down_state, "%s %f %s %f"," with state " , (*ZRoeShd)(0,it->first,ISH) , " ", (*ZRoeShd)(1,it->first,ISH));
        fprintf(RH_interp_down_state, "%s %f %s %f"," ",(*ZRoeShd)(2,it->first,ISH), " " ,(*ZRoeShd)(3,it->first,ISH));
        fprintf(RH_interp_down_state, "%s %f %s"," ",(*ZRoeShd)(4,it->first,ISH), " \n");

        point.at(counter_close)=it->first;
        distance.at(counter_close)=tr_dist.inverse_of_transformed_distance(it->second);
        counter_close ++;
      }


double dist;

for (unsigned I=0; I< distance.size(); I++)
{
  dist=dist+distance.at(I);
}

for (unsigned IV=0; IV<(*ndof); IV++){
for (unsigned I=0; I<distance.size();I++){
new_down_state.at(IV)=new_down_state.at(IV)+(distance.at(I)/dist)*((*ZRoeShd)(IV,point.at(I),ISH));
}
}


std::cout << " new down state for point " << IPOIN << " \n";

// for (unsigned IV=0; IV<(*ndof); IV++){ std::cout << " new_down_state.at(IV) " << new_down_state.at(IV) << " \n";}

for (unsigned IV=0; IV<(*ndof); IV++){
(*ZRoeShd)(IV,IPOIN,ISH)=new_down_state.at(IV);
}
  // de-allocate dynamic arrays
  freeArray();

}
//--------------------------------------------------------------------------//

void InterpShockPointStateRH::saveNewState(unsigned IPOIN)
{
//   unsigned ISH=0;
// for (unsigned IV=0; IV<(*ndof); IV++){
// (*ZRoeShu)(IV,IPOIN,ISH)=new_up_state.at(IV);
// }

}

//--------------------------------------------------------------------------//

void InterpShockPointStateRH::setAddress()
{
  unsigned start;

  start = PhysicsInfo::getnbDim() *
          (npoin->at(0) + 2 *
           PhysicsInfo::getnbShMax() * nShockPoints->at(0));
  XYZ = new Array2D <double> (PhysicsInfo::getnbDim(),
                             (npoin->at(1) + 2 *
                              PhysicsInfo::getnbShMax() *
                              nShockPoints->at(0)),
                             &coorVect->at(start));
  start = PhysicsInfo::getnbDofMax() *
          (npoin->at(0) + 2 *
           PhysicsInfo::getnbShMax() * nShockPoints->at(0));
  zroe = new Array2D <double>(PhysicsInfo::getnbDofMax(),
                              (npoin->at(1) + 2 *
                              PhysicsInfo::getnbShMax() *
                              nShockPoints->at(0)),
                              &zroeVect->at(start));
  start = (*nvt) * nelem->at(0);
  celnod = new Array2D<int> ((*nvt), nelem->at(1), &celnodVect->at(start));
  // XYZShu and XYZShd have the starting pointers referred to shocked mesh
  start = npoin->at(0) * PhysicsInfo::getnbDim() +
          PhysicsInfo::getnbDim() *
          (npoin->at(0) + 2 *
           PhysicsInfo::getnbShMax() * nShockPoints->at(0));
  XYZShu = new Array3D <double> (PhysicsInfo::getnbDim(),
                                nShockPoints->at(0),
                                PhysicsInfo::getnbShMax(),
                                &coorVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDim() +
          PhysicsInfo::getnbDim() *
          (npoin->at(0) + 2 *
           PhysicsInfo::getnbShMax() * nShockPoints->at(0));

  XYZShd = new Array3D <double> (PhysicsInfo::getnbDim(),
                                nShockPoints->at(0),
                                PhysicsInfo::getnbShMax(),
                                &coorVect->at(start));

  start = npoin->at(0) * PhysicsInfo::getnbDofMax();

  ZRoeShu =
    new Array3D <double> (PhysicsInfo::getnbDofMax(),
                          nShockPoints->at(0),
                          PhysicsInfo::getnbShMax(),
                          &zroeVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDofMax() +
          nShockPoints->at(0) *
          PhysicsInfo::getnbShMax() *
          PhysicsInfo::getnbDofMax();
  ZRoeShd =
    new Array3D <double> (PhysicsInfo::getnbDofMax(),
                          nShockPoints->at(0),
                          PhysicsInfo::getnbShMax(),
                          &zroeVect->at(start));

}

//--------------------------------------------------------------------------//

void InterpShockPointStateRH::freeArray()
{
  delete zroe; delete XYZ;
  delete XYZShu; delete XYZShd;
  delete ZRoeShu; delete ZRoeShd;
  delete celnod;

}

//--------------------------------------------------------------------------//

void InterpShockPointStateRH::setMeshData()
{
  nvt = MeshData::getInstance().getData <unsigned> ("NVT");
  npoin = MeshData::getInstance().getData <std::vector<unsigned> > ("NPOIN");
  nelem = MeshData::getInstance().getData <std::vector<unsigned> > ("NELEM");
  zroeVect = MeshData::getInstance().getData <std::vector<double> > ("ZROE");
  coorVect = MeshData::getInstance().getData <std::vector<double> > ("COOR");
  celnodVect =
    MeshData::getInstance().getData <std::vector<int> > ("CELNOD");
  nodcod =
    MeshData::getInstance().getData <std::vector<int> > ("NODCOD");
}

//--------------------------------------------------------------------------//

void InterpShockPointStateRH::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
   PhysicsData::getInstance().getData <std::vector<unsigned> > ("nShockPoints");
  XYZSh =
   PhysicsData::getInstance().getData <Array3D<double> > ("XYZSH");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
