// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cmath>
#include "VariableTransformerSF/Param2Prim.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/MeshData.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

Param2Prim::Param2Prim(const std::string& objectName) :
 VariableTransformer(objectName)
{
}

//--------------------------------------------------------------------------//

Param2Prim::~Param2Prim()
{
}

//--------------------------------------------------------------------------//

void Param2Prim::configure(OptionMap& cmap, const std::string& prefix)
{
  VariableTransformer::configure(cmap, prefix);
}

//--------------------------------------------------------------------------//

void Param2Prim::setup()
{
  LogToScreen(VERBOSE, "Param2Prim::setup() => start \n");

  LogToScreen(VERBOSE, "Param2Prim::setup() => end \n");
}

//--------------------------------------------------------------------------//

void Param2Prim::unsetup()
{
  LogToScreen(VERBOSE, "Param2Prim::unsetup()\n");
}

//--------------------------------------------------------------------------//

void Param2Prim::transform()	// capire in quale formato sono scritte
{
  LogToScreen(INFO, "Param2Prim::transform()\n");

  setPhysicsData();
  setMeshData();
  setAddress();

  double rhoref = ReferenceInfo::getpref() / ReferenceInfo::getTref() /
                  ReferenceInfo::getRgas();

  ReferenceInfo::setrhoref(rhoref);

  for(unsigned IPOIN=0; IPOIN<npoin->at(1); IPOIN++) {
    rho = (*zroe)(0,IPOIN) * (*zroe)(0,IPOIN);
    kinetic =0;
    for (unsigned NDOF=0; NDOF<PhysicsInfo::getnbDim(); NDOF++){    // both 2D and 3D
        kinetic = kinetic + (*zroe)(NDOF+2,IPOIN) * (*zroe)(NDOF+2,IPOIN);
    }
    kinetic = kinetic * 0.5;
    h = (*zroe)(1,IPOIN)/(*zroe)(0,IPOIN);
    help = (ReferenceInfo::getgam()-1)/ReferenceInfo::getgam();
    pres = help * (rho * h - kinetic);
      
    u = (*zroe)(2,IPOIN)/(*zroe)(0,IPOIN);
    v = (*zroe)(3,IPOIN)/(*zroe)(0,IPOIN);
      
    if (PhysicsInfo::getnbDim()==3){
        w = (*zroe)(4,IPOIN)/(*zroe)(0,IPOIN);
    }
      
    (*zroe)(0,IPOIN) = pres * ReferenceInfo::getrhoref() *
                       ReferenceInfo::geturef() * ReferenceInfo::geturef();
    (*zroe)(1,IPOIN) = u * ReferenceInfo::geturef();
    (*zroe)(2,IPOIN) = v * ReferenceInfo::geturef();
      
    if (PhysicsInfo::getnbDim()==3){
    (*zroe)(3,IPOIN) = w * ReferenceInfo::geturef();
    (*zroe)(4,IPOIN) = pres/rho *
       (ReferenceInfo::geturef() * ReferenceInfo::geturef() /
        ReferenceInfo::getRgas());
    }
      else
      {
          (*zroe)(3,IPOIN) = pres/rho *
                       (ReferenceInfo::geturef() * ReferenceInfo::geturef() /
                        ReferenceInfo::getRgas());
      }
    

    for(unsigned I=0; I<PhysicsInfo::getnbDim(); I++) {
     (*XYZ)(I,IPOIN) = (*XYZ)(I,IPOIN)*ReferenceInfo::getLref(); }
  }

  // de-allocate dynamic arrays
  freeArray();
}

//--------------------------------------------------------------------------//

void Param2Prim::setAddress()
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

void Param2Prim::freeArray()
{
  delete XYZ; delete zroe;
}

//--------------------------------------------------------------------------//

void Param2Prim::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> > ("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> > ("COOR");
}

//--------------------------------------------------------------------------//

void Param2Prim::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  nmol = PhysicsData::getInstance().getData <unsigned> ("NMOL");
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
