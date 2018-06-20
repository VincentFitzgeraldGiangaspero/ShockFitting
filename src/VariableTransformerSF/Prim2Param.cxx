// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

/// new implementation of 3D transformation: 09/2017


#include "VariableTransformerSF/Prim2Param.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/MeshData.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

Prim2Param::Prim2Param(const std::string& objectName) :
 VariableTransformer(objectName)
{
}

//--------------------------------------------------------------------------//

Prim2Param::~Prim2Param()
{
}

//--------------------------------------------------------------------------//

void Prim2Param::configure(OptionMap& cmap, const std::string& prefix)
{
  VariableTransformer::configure(cmap, prefix);
}
//--------------------------------------------------------------------------//

void Prim2Param::setup()
{
  LogToScreen(VERBOSE, "Prim2Param::setup() => start \n");

  LogToScreen(VERBOSE, "Prim2Param::setup() => end \n");
}

//--------------------------------------------------------------------------//

void Prim2Param::unsetup()
{
  LogToScreen(VERBOSE, "Prim2Param::unsetup()\n");
}

//--------------------------------------------------------------------------//

void Prim2Param::transform()
{
  LogToScreen(INFO, "Prim2Param::transform()\n");

  setPhysicsData();
  setMeshData();
  setAddress();

  double sqrtr;

  double rhoref = ReferenceInfo::getpref() / ReferenceInfo::getTref() /
                  ReferenceInfo::getRgas();

  ReferenceInfo::setrhoref(rhoref);


  if (PhysicsInfo::getnbDim()==3){ //3D

   for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
     rho = (*zroe)(0,IPOIN) / (*zroe)(4,IPOIN) / ReferenceInfo::getRgas();
     sqrtr = sqrt(rho / ReferenceInfo::getrhoref());
     kinetic = (*zroe)(1,IPOIN) * (*zroe)(1,IPOIN) +
               (*zroe)(2,IPOIN) * (*zroe)(2,IPOIN) +
               (*zroe)(3,IPOIN) * (*zroe)(3,IPOIN);

     kinetic = kinetic*0.5;
     help = ReferenceInfo::getgam()/(ReferenceInfo::getgam()-1);
     h = help * ReferenceInfo::getRgas() * (*zroe)(4,IPOIN) + kinetic;

     (*zroe)(4,IPOIN) = sqrtr * (*zroe)(3,IPOIN) / ReferenceInfo::geturef();
     (*zroe)(3,IPOIN) = sqrtr * (*zroe)(2,IPOIN) / ReferenceInfo::geturef();
     (*zroe)(2,IPOIN) = sqrtr * (*zroe)(1,IPOIN) / ReferenceInfo::geturef();
     (*zroe)(1,IPOIN) = sqrtr * h / (ReferenceInfo::geturef() *
                                     ReferenceInfo::geturef());
     (*zroe)(0,IPOIN) = sqrtr;


    for(unsigned I=0; I<PhysicsInfo::getnbDim(); I++) {
     (*XYZ)(I,IPOIN) = (*XYZ)(I,IPOIN)/ReferenceInfo::getLref(); }
     
    }
  } 
  else { //2D

    for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
     rho = (*zroe)(0,IPOIN) / (*zroe)(3,IPOIN) / ReferenceInfo::getRgas();
     sqrtr = sqrt(rho / ReferenceInfo::getrhoref());
     kinetic = (*zroe)(1,IPOIN) * (*zroe)(1,IPOIN) +
               (*zroe)(2,IPOIN) * (*zroe)(2,IPOIN);
     kinetic = kinetic*0.5;
     help = ReferenceInfo::getgam()/(ReferenceInfo::getgam()-1);
     h = help * ReferenceInfo::getRgas() * (*zroe)(3,IPOIN) + kinetic;

     (*zroe)(3,IPOIN) = sqrtr * (*zroe)(2,IPOIN) / ReferenceInfo::geturef();
     (*zroe)(2,IPOIN) = sqrtr * (*zroe)(1,IPOIN) / ReferenceInfo::geturef();
     (*zroe)(1,IPOIN) = sqrtr * h / (ReferenceInfo::geturef() *
                                     ReferenceInfo::geturef());
     (*zroe)(0,IPOIN) = sqrtr;
    

    for(unsigned I=0; I<PhysicsInfo::getnbDim(); I++) {
     (*XYZ)(I,IPOIN) = (*XYZ)(I,IPOIN)/ReferenceInfo::getLref(); }
     }
  }
    

  

  // de-allocate dynamic arrays
  freeArray();
}

//--------------------------------------------------------------------------//

void Prim2Param::setAddress()
{
  totsize = npoin->at(0) + npoin->at(1) + 4 *
            PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax();
  zroeVect->resize(PhysicsInfo::getnbDofMax() * totsize);
  coorVect->resize(PhysicsInfo::getnbDim() * totsize);
  
  start = PhysicsInfo::getnbDim() * 
          (npoin->at(0) + 2 *
           PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax());
  XYZ = new Array2D <double> (PhysicsInfo::getnbDim(),
                             (npoin->at(1) + 2 *
                              PhysicsInfo::getnbShMax() *
                              PhysicsInfo::getnbShPointsMax()),
                              &coorVect->at(start));
  start = PhysicsInfo::getnbDofMax() * 
          (npoin->at(0) + 2 *
           PhysicsInfo::getnbShMax() * PhysicsInfo::getnbShPointsMax());
  zroe = new Array2D <double> (PhysicsInfo::getnbDofMax(),
                               (npoin->at(1) + 2 *
                                PhysicsInfo::getnbShMax() *
                                PhysicsInfo::getnbShPointsMax()),
                                &zroeVect->at(start));
}

//--------------------------------------------------------------------------//

void Prim2Param::freeArray()
{
  delete XYZ; delete zroe;
}

//--------------------------------------------------------------------------//

void Prim2Param::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> > ("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> > ("COOR");
}

//--------------------------------------------------------------------------//

void Prim2Param::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  IE = PhysicsData::getInstance().getData <unsigned> ("IE");
  IEV = PhysicsData::getInstance().getData <unsigned> ("IEV");
  IX = PhysicsData::getInstance().getData <unsigned> ("IX");
  IY = PhysicsData::getInstance().getData <unsigned> ("IY");
  mm = PhysicsData::getInstance().getData <vector<double> > ("MM");
  hf = PhysicsData::getInstance().getData <vector<double> > ("HF");
  thev = PhysicsData::getInstance().getData <vector<double> > ("THEV");
  gams = PhysicsData::getInstance().getData <vector<double> > ("GAMS");
  typemol =
    PhysicsData::getInstance().getData <vector<string> > ("TYPEMOL");
}

//--------------------------------------------------------------------------//

void Prim2Param::setPhysicsData(string FirstCaptured)
{
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  IE = PhysicsData::getInstance().getData <unsigned> ("IE");
  IEV = PhysicsData::getInstance().getData <unsigned> ("IEV");
  IX = PhysicsData::getInstance().getData <unsigned> ("IX");
  IY = PhysicsData::getInstance().getData <unsigned> ("IY");
  mm = PhysicsData::getInstance().getData <vector<double> > ("MM");
  hf = PhysicsData::getInstance().getData <vector<double> > ("HF");
  thev = PhysicsData::getInstance().getData <vector<double> > ("THEV");
  gams = PhysicsData::getInstance().getData <vector<double> > ("GAMS");
  typemol =
    PhysicsData::getInstance().getData <vector<string> > ("TYPEMOL");
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
