// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

//Vincent 18/04/2017

#include "StateUpdaterSF/InterpShockPointStateRH.hh"
#include "StateUpdaterSF/ComputeStateDps4Pg3DAssignUpstream.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "SConfig/ObjectProvider.hh"
#include "StateUpdaterSF/CoDc.hh"
#include "StateUpdaterSF/CoShock.hh"


//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<ComputeStateDps4Pg3DAssignUpstream, StateUpdater>
 computeStateDpsPg3DAssignUpstreamProv("ComputeStateDps4Pg3DAssignUpstream");

//----------------------------------------------------------------------------//

ComputeStateDps4Pg3DAssignUpstream::ComputeStateDps4Pg3DAssignUpstream(const std::string& objectName) :
 StateUpdater(objectName)
{
}

//----------------------------------------------------------------------------//

ComputeStateDps4Pg3DAssignUpstream::~ComputeStateDps4Pg3DAssignUpstream()
{
}

//----------------------------------------------------------------------------//

void ComputeStateDps4Pg3DAssignUpstream::setup()
{
  LogToScreen(VERBOSE,"ComputeStateDps4Pg3DAssignUpstream::setup() => start\n");

  LogToScreen(VERBOSE,"ComputeStateDps4Pg3DAssignUpstream::setup() => end\n");
}

//----------------------------------------------------------------------------//

void ComputeStateDps4Pg3DAssignUpstream::unsetup()
{
  LogToScreen(VERBOSE,"ComputeStateDps4Pg3DAssignUpstream::unsetup()\n");
}

//----------------------------------------------------------------------------//

void ComputeStateDps4Pg3DAssignUpstream::update()
{
  LogToScreen(INFO,"ComputeStateDps4Pg3DAssignUpstream::update()\n");

  logfile.Open(getClassName().c_str());

  setMeshData();
  setPhysicsData();

  setAddress();

  setDiscSpeedSize();

  unsigned I;

  std::vector<unsigned> Need_Interpol;

  unsigned counter_interpol=0;

  unsigned counter_no_conv=0;



  // create object of CoShock class
  CoShock computenewStateForShock;

  // create object of CoDc class
  CoDc computenewStateForDc;    // not used right now

  InterpShockPointStateRH newState;

  xd.resize(4);
  xu.resize(4);


  // xd.resize(5);   // 5 in 3D
  // xu.resize(5);

  for(unsigned ISH=0; ISH<(*nShocks); ISH++) {

   unsigned ivalue = ISH+1;
   logfile("Shock/Disc. n. ",ivalue, "\n");

   for(unsigned IV=0; IV<nShockPoints->at(ISH); IV++) {
    ++TotnbShockPoints;

    I=IV;
    dx = (*vShNor)(0,I,ISH);
    dy = (*vShNor)(1,I,ISH);
    dz = (*vShNor)(2,I,ISH);

    // upload downstream status
    recoverDownState(IV,ISH);

    // upload upstream status
    recoverUpState(IV,ISH);

    // initialize discontinuity speed
    WS = 0.0;


if ((*ZRoeShu)(0,IV,ISH)>(*ZRoeShd)(0,IV,ISH)){
  std::cout << "CONTROLLO STATI IN ComputeShockPointsState3D \n";
  std::cout << "Shock Point n  " << IV  << " with swapped UP and DOWN \n";
  std::cout << "P UP " << (*ZRoeShu)(0,IV,ISH) << '\n';
  std::cout << "P DOWN " << (*ZRoeShd)(0,IV,ISH) << '\n';

}


    //  if (ivalue== 303 || ivalue== 334 || ivalue== 847 ){
    //   std::cout << "---------------------------------------------------------" << '\n';
    //   std::cout << "SH point n " << ivalue << " W " << WS << '\n';
    //   std::cout << "dx " << dx << " dy " << dy << " dz " << dz << '\n';
    //   std::cout << "dxt " << dxt << " dyt " << dyt << " dzt " << dzt << '\n';
    //   std::cout << '\n';
     //
    //   std::cout << "xd.at(0) " << xd.at(0) << " xd.at(1) " << xd.at(1) << " xd.at(2) " << xd.at(2) << " xd.at(3) " << xd.at(3) << '\n';
    //   std::cout << "xu.at(0) " << xu.at(0) << " xu.at(1) " << xu.at(1) << " xu.at(2) " << xu.at(2) << " xu.at(3) " << xu.at(3) << '\n';
    //   std::cout << "R2 " << R2 << '\n';
    //   std::cout << '\n';
     //
    //   std::cout << "NORM module " << sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2)) << '\n';
    //   std::cout << "TANG module " << sqrt(pow(dxt,2)+pow(dyt,2)+pow(dzt,2)) << '\n';
    //   std::cout << '\n';
     //
    //   std::cout << "UVW module down " << sqrt(pow((*ZRoeShd)(1,IV,ISH),2)+pow((*ZRoeShd)(2,IV,ISH),2)+pow((*ZRoeShd)(3,IV,ISH),2)) << '\n';
    //   std::cout << "N-T module down " << sqrt(pow(xd.at(2),2)+pow(xd.at(3),2)) << '\n';
    //   std::cout << "UVW module up " << sqrt(pow((*ZRoeShu)(1,IV,ISH),2)+pow((*ZRoeShu)(2,IV,ISH),2)+pow((*ZRoeShu)(3,IV,ISH),2)) << '\n';
    //   std::cout << "N-T module up " << sqrt(pow(xu.at(2),2)+pow(xu.at(3),2)) << '\n';
    //   std::cout << '\n';
     //
    //   std::cout << "(*ZRoeShdOld)(0,IV,ISH) " << (*ZRoeShd)(0,IV,ISH) << '\n';
    //   std::cout << "(*ZRoeShdOld)(1,IV,ISH) " << (*ZRoeShd)(1,IV,ISH) << '\n';
    //   std::cout << "(*ZRoeShdOld)(2,IV,ISH) " << (*ZRoeShd)(2,IV,ISH) << '\n';
    //   std::cout << "(*ZRoeShdOld)(3,IV,ISH) " << (*ZRoeShd)(3,IV,ISH) << '\n';
    //   std::cout << "(*ZRoeShdOld)(4,IV,ISH) " << (*ZRoeShd)(4,IV,ISH) << '\n';
    //   std::cout << '\n';
     //
    //   std::cout << "(*ZRoeShuOld)(0,IV,ISH) " << (*ZRoeShu)(0,IV,ISH) << '\n';
    //   std::cout << "(*ZRoeShuOld)(1,IV,ISH) " << (*ZRoeShu)(1,IV,ISH) << '\n';
    //   std::cout << "(*ZRoeShuOld)(2,IV,ISH) " << (*ZRoeShu)(2,IV,ISH) << '\n';
    //   std::cout << "(*ZRoeShuOld)(3,IV,ISH) " << (*ZRoeShu)(3,IV,ISH) << '\n';
    //   std::cout << "(*ZRoeShuOld)(4,IV,ISH) " << (*ZRoeShu)(4,IV,ISH) << '\n';
    //   std::cout << '\n';
     //
    //   std::cout << "---------------------------------------------------------" << '\n';
     //
    //   std::cout << '\n';
//}

    if(typeSh->at(ISH)=="S") {
     computenewStateForShock.callCoShock(xd,xu,R2);
     xd = computenewStateForShock.getnewDownValues();
     WS = computenewStateForShock.getnewDiscSpeed();
    }


    if(typeSh->at(ISH)=="D") {
     computenewStateForDc.callCoDc(xd,xu);
     xd = computenewStateForDc.getnewDownValues();
     xu = computenewStateForDc.getnewUpValues();
     WS = computenewStateForDc.getnewDiscSpeed();
    }

    // enforce tangential component equality for the shock case
    if(typeSh->at(ISH)=="S") { xd.at(3) = xu.at(3); }

    // save old downstream status
    saveDownState(IV, ISH);


//     if (WS!=WS){
//       std::cout << "---------------------------------------------------------" << '\n';
//       std::cout << "SH point n " << ivalue << " W " << WS << '\n';
//       std::cout << "dx " << dx << " dy " << dy << " dz " << dz << '\n';
//       std::cout << "dxt " << dxt << " dyt " << dyt << " dzt " << dzt << '\n';
//       std::cout << "SH COORDINATES \n";
//       std::cout << "x= " << (*XYZSh)(0,IV,ISH) << " \n";
//       std::cout << "y= " << (*XYZSh)(1,IV,ISH) << " \n";
//       std::cout << "z= " << (*XYZSh)(2,IV,ISH) << " \n";
//
//
//       std::cout << '\n';
//
//       std::cout << "xd.at(0) " << xd.at(0) << " xd.at(1) " << xd.at(1) << " xd.at(2) " << xd.at(2) << " xd.at(3) " << xd.at(3) << '\n';
//       std::cout << "xu.at(0) " << xu.at(0) << " xu.at(1) " << xu.at(1) << " xu.at(2) " << xu.at(2) << " xu.at(3) " << xu.at(3) << '\n';
//       std::cout << "R2 " << R2 << '\n';
//       std::cout << '\n';
//
//       std::cout << "NORM module " << sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2)) << '\n';
//       std::cout << "TANG module " << sqrt(pow(dxt,2)+pow(dyt,2)+pow(dzt,2)) << '\n';
//       std::cout << '\n';
//
//       std::cout << "UVW module down " << sqrt(pow((*ZRoeShdOld)(1,IV,ISH),2)+pow((*ZRoeShdOld)(2,IV,ISH),2)+pow((*ZRoeShdOld)(3,IV,ISH),2)) << '\n';
//       std::cout << "N-T module down " << sqrt(pow(xd.at(2),2)+pow(xd.at(3),2)) << '\n';
//       std::cout << "UVW module up " << sqrt(pow((*ZRoeShuOld)(1,IV,ISH),2)+pow((*ZRoeShuOld)(2,IV,ISH),2)+pow((*ZRoeShuOld)(3,IV,ISH),2)) << '\n';
//       std::cout << "N-T module up " << sqrt(pow(xu.at(2),2)+pow(xu.at(3),2)) << '\n';
//       std::cout << '\n';
//
//       std::cout << "(*ZRoeShdOld)(0,IV,ISH) " << (*ZRoeShdOld)(0,IV,ISH) << '\n';
//       std::cout << "(*ZRoeShdOld)(1,IV,ISH) " << (*ZRoeShdOld)(1,IV,ISH) << '\n';
//       std::cout << "(*ZRoeShdOld)(2,IV,ISH) " << (*ZRoeShdOld)(2,IV,ISH) << '\n';
//       std::cout << "(*ZRoeShdOld)(3,IV,ISH) " << (*ZRoeShdOld)(3,IV,ISH) << '\n';
//       std::cout << "(*ZRoeShdOld)(4,IV,ISH) " << (*ZRoeShdOld)(4,IV,ISH) << '\n';
//       std::cout << '\n';
//       // std::cout << "(*ZRoeShuOld)(0,IV,ISH) " << (*ZRoeShuOld)(0,IV,ISH) << '\n';
//       // std::cout << "(*ZRoeShuOld)(1,IV,ISH) " << (*ZRoeShuOld)(1,IV,ISH) << '\n';
//       // std::cout << "(*ZRoeShuOld)(2,IV,ISH) " << (*ZRoeShuOld)(2,IV,ISH) << '\n';
//       // std::cout << "(*ZRoeShuOld)(3,IV,ISH) " << (*ZRoeShuOld)(3,IV,ISH) << '\n';
//       // std::cout << "(*ZRoeShuOld)(4,IV,ISH) " << (*ZRoeShuOld)(4,IV,ISH) << '\n';
//       std::cout << "(*ZRoeShu)(0,IV,ISH) " << (*ZRoeShu)(0,IV,ISH) << '\n';
//       std::cout << "(*ZRoeShu)(1,IV,ISH) " << (*ZRoeShu)(1,IV,ISH) << '\n';
//       std::cout << "(*ZRoeShu)(2,IV,ISH) " << (*ZRoeShu)(2,IV,ISH) << '\n';
//       std::cout << "(*ZRoeShu)(3,IV,ISH) " << (*ZRoeShu)(3,IV,ISH) << '\n';
//       std::cout << "(*ZRoeShu)(4,IV,ISH) " << (*ZRoeShu)(4,IV,ISH) << '\n';
//
//       std::cout << "---------------------------------------------------------" << '\n';
//       std::cout << '\n';
// }

    // compute downstream variables and assing the new values
    // to ZRoe array (for the shock case)
    computeDownState(IV, ISH);

    // compute upstream variables and assing the new values
    // to ZRoe array (for the shock case)
    computeUpState(IV, ISH);

    // set the new discontinuity speed
    (*WSh)(0,IV,ISH) = WS * dx;
    (*WSh)(1,IV,ISH) = WS * dy;
    (*WSh)(2,IV,ISH) = WS * dz;


    if (WS!=WS){

      // std::cout << "---------------------------------------------------------" << '\n';
      // std::cout << "SH point n " << ivalue << " W " << WS << '\n';
      // std::cout << "dx " << dx << " dy " << dy << " dz " << dz << '\n';
      // std::cout << "dxt " << dxt << " dyt " << dyt << " dzt " << dzt << '\n';
      // std::cout << "SH COORDINATES \n";
      // std::cout << "x= " << (*XYZSh)(0,IV,ISH) << " \n";
      // std::cout << "y= " << (*XYZSh)(1,IV,ISH) << " \n";
      // std::cout << "z= " << (*XYZSh)(2,IV,ISH) << " \n";
      //
      //
      // std::cout << '\n';
      //
      // std::cout << "xd.at(0) " << xd.at(0) << " xd.at(1) " << xd.at(1) << " xd.at(2) " << xd.at(2) << " xd.at(3) " << xd.at(3) << '\n';
      // std::cout << "xu.at(0) " << xu.at(0) << " xu.at(1) " << xu.at(1) << " xu.at(2) " << xu.at(2) << " xu.at(3) " << xu.at(3) << '\n';
      // std::cout << "R2 " << R2 << '\n';
      // std::cout << '\n';
      //
      // std::cout << "NORM module " << sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2)) << '\n';
      // std::cout << "TANG module " << sqrt(pow(dxt,2)+pow(dyt,2)+pow(dzt,2)) << '\n';
      // std::cout << '\n';
      //
      // std::cout << "UVW module down " << sqrt(pow((*ZRoeShdOld)(1,IV,ISH),2)+pow((*ZRoeShdOld)(2,IV,ISH),2)+pow((*ZRoeShdOld)(3,IV,ISH),2)) << '\n';
      // std::cout << "N-T module down " << sqrt(pow(xd.at(2),2)+pow(xd.at(3),2)) << '\n';
      // std::cout << "UVW module up " << sqrt(pow((*ZRoeShuOld)(1,IV,ISH),2)+pow((*ZRoeShuOld)(2,IV,ISH),2)+pow((*ZRoeShuOld)(3,IV,ISH),2)) << '\n';
      // std::cout << "N-T module up " << sqrt(pow(xu.at(2),2)+pow(xu.at(3),2)) << '\n';
      // std::cout << '\n';
      //
      // std::cout << "(*ZRoeShdOld)(0,IV,ISH) " << (*ZRoeShdOld)(0,IV,ISH) << '\n';
      // std::cout << "(*ZRoeShdOld)(1,IV,ISH) " << (*ZRoeShdOld)(1,IV,ISH) << '\n';
      // std::cout << "(*ZRoeShdOld)(2,IV,ISH) " << (*ZRoeShdOld)(2,IV,ISH) << '\n';
      // std::cout << "(*ZRoeShdOld)(3,IV,ISH) " << (*ZRoeShdOld)(3,IV,ISH) << '\n';
      // std::cout << "(*ZRoeShdOld)(4,IV,ISH) " << (*ZRoeShdOld)(4,IV,ISH) << '\n';
      // std::cout << '\n';
      // // std::cout << "(*ZRoeShuOld)(0,IV,ISH) " << (*ZRoeShuOld)(0,IV,ISH) << '\n';
      // // std::cout << "(*ZRoeShuOld)(1,IV,ISH) " << (*ZRoeShuOld)(1,IV,ISH) << '\n';
      // // std::cout << "(*ZRoeShuOld)(2,IV,ISH) " << (*ZRoeShuOld)(2,IV,ISH) << '\n';
      // // std::cout << "(*ZRoeShuOld)(3,IV,ISH) " << (*ZRoeShuOld)(3,IV,ISH) << '\n';
      // // std::cout << "(*ZRoeShuOld)(4,IV,ISH) " << (*ZRoeShuOld)(4,IV,ISH) << '\n';
      // std::cout << "(*ZRoeShu)(0,IV,ISH) " << (*ZRoeShu)(0,IV,ISH) << '\n';
      // std::cout << "(*ZRoeShu)(1,IV,ISH) " << (*ZRoeShu)(1,IV,ISH) << '\n';
      // std::cout << "(*ZRoeShu)(2,IV,ISH) " << (*ZRoeShu)(2,IV,ISH) << '\n';
      // std::cout << "(*ZRoeShu)(3,IV,ISH) " << (*ZRoeShu)(3,IV,ISH) << '\n';
      // std::cout << "(*ZRoeShu)(4,IV,ISH) " << (*ZRoeShu)(4,IV,ISH) << '\n';
      //
      // std::cout << "---------------------------------------------------------" << '\n';
      // std::cout << '\n';
}


if (WS!=WS){

  logfile(" INTERPOLATION UPSTREAM FOR SHOCK POINT ", ivalue, " \n");

  counter_interpol++;
  newState.InterpUpState(IV);

  // upload upstream status
  recoverUpState(IV,ISH);

  // upload downstream status
  recoverDownOlDState(IV,ISH);


  computenewStateForShock.callCoShock(xd,xu,R2);
  xd = computenewStateForShock.getnewDownValues();
  WS = computenewStateForShock.getnewDiscSpeed();


  // compute downstream variables and assing the new values
  // to ZRoe array (for the shock case)
  computeDownState(IV, ISH);

  // compute upstream variables and assing the new values
  // to ZRoe array (for the shock case)
  computeUpState(IV, ISH);

  // set the new discontinuity speed
  (*WSh)(0,IV,ISH) = WS * dx;
  (*WSh)(1,IV,ISH) = WS * dy;
  (*WSh)(2,IV,ISH) = WS * dz;

  // std::cout << "---------------------------------------------------------" << '\n';
  // std::cout << " NEW STATE INTERPOLATED " << '\n';
  //
  // std::cout << "---------------------------------------------------------" << '\n';
  // std::cout << "SH point n " << ivalue << " W " << WS << '\n';
  // std::cout << "dx " << dx << " dy " << dy << " dz " << dz << '\n';
  // std::cout << "dxt " << dxt << " dyt " << dyt << " dzt " << dzt << '\n';
  // std::cout << "SH COORDINATES \n";
  // std::cout << "x= " << (*XYZSh)(0,IV,ISH) << " \n";
  // std::cout << "y= " << (*XYZSh)(1,IV,ISH) << " \n";
  // std::cout << "z= " << (*XYZSh)(2,IV,ISH) << " \n";
  //
  //
  // std::cout << '\n';
  //
  // std::cout << "xd.at(0) " << xd.at(0) << " xd.at(1) " << xd.at(1) << " xd.at(2) " << xd.at(2) << " xd.at(3) " << xd.at(3) << '\n';
  // std::cout << "xu.at(0) " << xu.at(0) << " xu.at(1) " << xu.at(1) << " xu.at(2) " << xu.at(2) << " xu.at(3) " << xu.at(3) << '\n';
  // std::cout << "R2 " << R2 << '\n';
  // std::cout << '\n';
  //
  // std::cout << "NORM module " << sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2)) << '\n';
  // std::cout << "TANG module " << sqrt(pow(dxt,2)+pow(dyt,2)+pow(dzt,2)) << '\n';
  // std::cout << '\n';
  //
  // std::cout << "UVW module down " << sqrt(pow((*ZRoeShdOld)(1,IV,ISH),2)+pow((*ZRoeShdOld)(2,IV,ISH),2)+pow((*ZRoeShdOld)(3,IV,ISH),2)) << '\n';
  // std::cout << "N-T module down " << sqrt(pow(xd.at(2),2)+pow(xd.at(3),2)) << '\n';
  // std::cout << "UVW module up " << sqrt(pow((*ZRoeShuOld)(1,IV,ISH),2)+pow((*ZRoeShuOld)(2,IV,ISH),2)+pow((*ZRoeShuOld)(3,IV,ISH),2)) << '\n';
  // std::cout << "N-T module up " << sqrt(pow(xu.at(2),2)+pow(xu.at(3),2)) << '\n';
  // std::cout << '\n';
  //
  // std::cout << "(*ZRoeShd)(0,IV,ISH) " << (*ZRoeShd)(0,IV,ISH) << '\n';
  // std::cout << "(*ZRoeShd)(1,IV,ISH) " << (*ZRoeShd)(1,IV,ISH) << '\n';
  // std::cout << "(*ZRoeShd)(2,IV,ISH) " << (*ZRoeShd)(2,IV,ISH) << '\n';
  // std::cout << "(*ZRoeShd)(3,IV,ISH) " << (*ZRoeShd)(3,IV,ISH) << '\n';
  // std::cout << "(*ZRoeShd)(4,IV,ISH) " << (*ZRoeShd)(4,IV,ISH) << '\n';
  // std::cout << '\n';
  // std::cout << "(*ZRoeShu)(0,IV,ISH) " << (*ZRoeShu)(0,IV,ISH) << '\n';
  // std::cout << "(*ZRoeShu)(1,IV,ISH) " << (*ZRoeShu)(1,IV,ISH) << '\n';
  // std::cout << "(*ZRoeShu)(2,IV,ISH) " << (*ZRoeShu)(2,IV,ISH) << '\n';
  // std::cout << "(*ZRoeShu)(3,IV,ISH) " << (*ZRoeShu)(3,IV,ISH) << '\n';
  // std::cout << "(*ZRoeShu)(4,IV,ISH) " << (*ZRoeShu)(4,IV,ISH) << '\n';
  //
  // std::cout << "---------------------------------------------------------" << '\n';
  // std::cout << "WSh_x " << (*WSh)(0,IV,ISH)  << "\n";
  // std::cout << "WSh_y " << (*WSh)(1,IV,ISH)  << "\n";
  // std::cout << "WSh_z " << (*WSh)(2,IV,ISH)  << "\n";
  // std::cout << "WSh module " << sqrt(pow((*WSh)(0,IV,ISH),2)+pow((*WSh)(1,IV,ISH),2)+pow((*WSh)(2,IV,ISH),2))  << "\n";
  //
  // std::cout << '\n';

}


    ivalue = IV+1;
    logfile("S/D point nr. ", ivalue, " Speed: ");
    logfile((*WSh)(0,IV,ISH), ", ", (*WSh)(1,IV,ISH),", ", (*WSh)(2,IV,ISH), "\n");
    logfile("WS ", WS, " \n");


    if(WS!=WS){
      Need_Interpol.push_back(IV);
      counter_no_conv++;

      (*WSh)(0,IV,ISH) = 0;
      (*WSh)(1,IV,ISH) = 0;
      (*WSh)(2,IV,ISH) = 0;
    }

   }
  }

  std::cout << "Need_Interpol.size() " <<Need_Interpol.size() <<'\n';
  std::cout << " counter_no_conv " << counter_no_conv<<'\n';




  for(unsigned IV=0; IV<Need_Interpol.size(); IV++) {

    std::cout << "INTERPOLATION UP and DOWN for Shock Point " << Need_Interpol[IV] << "\n";

    unsigned ISH=0;
    logfile(" INTERPOLATION UPSTREAM AND DOWNSTREAM FOR SHOCK POINT ", Need_Interpol[IV], " \n");

    InterpDownState(Need_Interpol[IV]);
    InterpUpState(Need_Interpol[IV]);
    InterpShockVelocity(Need_Interpol[IV],ISH);

  //  std::cout << "(*ZRoeShd)(0,IV,ISH) " << (*ZRoeShd)(0,IV,ISH) << '\n';
  //  std::cout << "(*ZRoeShd)(1,IV,ISH) " << (*ZRoeShd)(1,IV,ISH) << '\n';
  //  std::cout << "(*ZRoeShd)(2,IV,ISH) " << (*ZRoeShd)(2,IV,ISH) << '\n';
  //  std::cout << "(*ZRoeShd)(3,IV,ISH) " << (*ZRoeShd)(3,IV,ISH) << '\n';
  //  std::cout << "(*ZRoeShd)(4,IV,ISH) " << (*ZRoeShd)(4,IV,ISH) << '\n';
  //  std::cout << '\n';
  //  std::cout << "(*ZRoeShu)(0,IV,ISH) " << (*ZRoeShu)(0,IV,ISH) << '\n';
  //  std::cout << "(*ZRoeShu)(1,IV,ISH) " << (*ZRoeShu)(1,IV,ISH) << '\n';
  //  std::cout << "(*ZRoeShu)(2,IV,ISH) " << (*ZRoeShu)(2,IV,ISH) << '\n';
  //  std::cout << "(*ZRoeShu)(3,IV,ISH) " << (*ZRoeShu)(3,IV,ISH) << '\n';
  //  std::cout << "(*ZRoeShu)(4,IV,ISH) " << (*ZRoeShu)(4,IV,ISH) << '\n';

  //  std::cout << "---------------------------------------------------------" << '\n';
  //  std::cout << "WSh_x " << (*WSh)(0,IV,ISH)  << "\n";
  //  std::cout << "WSh_y " << (*WSh)(1,IV,ISH)  << "\n";
  //  std::cout << "WSh_z " << (*WSh)(2,IV,ISH)  << "\n";
  //  std::cout << "WSh module " << sqrt(pow((*WSh)(0,IV,ISH),2)+pow((*WSh)(1,IV,ISH),2)+pow((*WSh)(2,IV,ISH),2))  << "\n";

  //  std::cout << '\n';

  }

  std::cout << " -----> interpolation for " << counter_interpol <<" shock points \n";

  std::cout << " -----> no convergence for " << counter_no_conv <<" shock points \n";


//for(unsigned IV=0; IV<Need_Interpol.size(); IV++) {
//    for (unsigned IA=0; IA<(*ndof); IA++){
//      if ((*ZRoeShd)(IA,IV,0)!=(*ZRoeShd)(IA,IV,0)){
//        std::cout << " NO CONVERGENCE FOR " << IV << " POINT \n";
//      }
//    }
//  }

  // de-allocate dynamic array
  freeArray();

  logfile.Close();
}

//----------------------------------------------------------------------------//

void ComputeStateDps4Pg3DAssignUpstream::InterpUpState(unsigned POINT)
{



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

 for (unsigned IV=0; IV<(*ndof); IV++){ std::cout << " new_up_state.at(IV) " << new_up_state.at(IV) << " \n";}

for (unsigned IV=0; IV<(*ndof); IV++){
(*ZRoeShu)(IV,POINT,ISH)=new_up_state.at(IV);
}
  
}

//----------------------------------------------------------------------------//

void ComputeStateDps4Pg3DAssignUpstream::InterpDownState(unsigned POINT)
{

  
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
for (unsigned I=1; I<distance.size();I++){
if ((*ZRoeShd)(IV,point.at(I),ISH)==(*ZRoeShd)(IV,point.at(I),ISH)){
new_down_state.at(IV)=new_down_state.at(IV)+(distance.at(I)/dist)*((*ZRoeShd)(IV,point.at(I),ISH));
}
}
}


std::cout << " new down state for point " << IPOIN << " \n";

 for (unsigned IV=0; IV<(*ndof); IV++){ std::cout << " new_down_state.at(IV) " << new_down_state.at(IV) << " \n";}

for (unsigned IV=0; IV<(*ndof); IV++){
(*ZRoeShd)(IV,IPOIN,ISH)=new_down_state.at(IV);
}

}
//--------------------------------------------------------------------------//

void ComputeStateDps4Pg3DAssignUpstream::InterpShockVelocity(unsigned POINT, unsigned ISH){

  // std::vector<Point_3> shock_points_d(nShockPoints->at(ISH));
  std::vector<Point_3> shock_points(nShockPoints->at(ISH));

  // leggi gli shock_points downstream e upstream
  for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
    shock_points.at(ISHPOIN)=Point_3((*XYZSh)(0,ISHPOIN,ISH),(*XYZSh)(1,ISHPOIN,ISH),(*XYZSh)(2,ISHPOIN,ISH));
  }

std::cout << "here0 " << '\n';
  // K-neighborous search from CGAL, to get back the ID of the closest node of the mesh
  My_point_property_map ppmap(shock_points);

  // Insert number_of_data_points in the tree
  Tree tree(
    boost::counting_iterator<std::size_t>(0),
    boost::counting_iterator<std::size_t>(shock_points.size()),
    Splitter(),
    Traits(ppmap));

    Distance tr_dist(ppmap);
    // number of neighbors points
    const unsigned int K=10;

    std::vector<double> new_WSh(3);

    std::vector<double> distance(K);
    std::vector<double> point(K);

    unsigned counter_close=0;

    std::cout << "here 1" << '\n';

    // working variables and vectors

    FILE* RH_interp_ShockVelocity;
    RH_interp_ShockVelocity = fopen("log/RH_interp_ShockVelocity.dat","a");

    fprintf(RH_interp_ShockVelocity,"%s","TITLE = RH_interp_ShockVelocity n");

    K_Neighbor_search search(tree, shock_points[POINT], K,0,true,tr_dist);

    fprintf(RH_interp_ShockVelocity,"%s","------------------------------------------\n");
    fprintf(RH_interp_ShockVelocity,"%s %u %s %u %s " ,"The ", K , " nearest upstream shock point to the upstream shock point ",
    POINT+1, " with coordinates ");
    fprintf(RH_interp_ShockVelocity, "%22.14f %22.14f %22.14f %s" ,(*XYZSh)(0,POINT,ISH) , (*XYZSh)(1,POINT,ISH) , (*XYZSh)(2,POINT,ISH) , "\n");

    for(K_Neighbor_search::iterator it = search.begin(); it != search.end(); it++){

      fprintf(RH_interp_ShockVelocity,"%s %u %s %u %s %f %s %f %s %f %s"," -> ", counter_close , " -near-point ",
      it->first, " : ", (*XYZSh)(0,it->first,ISH), " ", (*XYZSh)(1,it->first,ISH), " ", (*XYZSh)(2,it->first,ISH),"\n");
      fprintf(RH_interp_ShockVelocity, "%s %f %s"," at distance " , tr_dist.inverse_of_transformed_distance(it->second) , "\n");
      fprintf(RH_interp_ShockVelocity, "%s %f %s %f"," with state " , (*ZRoeShu)(0,it->first,ISH) , " ", (*ZRoeShu)(1,it->first,ISH));
      fprintf(RH_interp_ShockVelocity, "%s %f %s %f"," ",(*ZRoeShu)(2,it->first,ISH), " " ,(*ZRoeShu)(3,it->first,ISH));
      fprintf(RH_interp_ShockVelocity, "%s %f %s"," ",(*ZRoeShu)(4,it->first,ISH), " \n");
      fprintf(RH_interp_ShockVelocity, "%s %f %s","Shock Velocity WSh =",(*WSh)(0,it->first,ISH), " ");
      fprintf(RH_interp_ShockVelocity, "%f %s %f %s",(*WSh)(1,it->first,ISH), " ", (*WSh)(2,it->first,ISH), "\n");

      point.at(counter_close)=it->first;
      distance.at(counter_close)=tr_dist.inverse_of_transformed_distance(it->second);
      counter_close ++;
    }

std::cout << "here 2" << '\n';

    double dist;

    for (unsigned I=0; I< distance.size(); I++)
    {
      dist=dist+distance.at(I);
    }

    for (unsigned IV=0; IV<3; IV++){
      for (unsigned I=0; I<distance.size();I++){
        new_WSh.at(IV)=new_WSh.at(IV)+(distance.at(I)/dist)*((*WSh)(IV,point.at(I),ISH));
      }
    }


    std::cout << " new Shock velocity for point " << POINT << " WSh "<<new_WSh.at(0) <<" \n";


    for (unsigned IV=0; IV<3; IV++){
      (*WSh)(IV,POINT,ISH)=new_WSh.at(IV);
    }

std::cout << "here end" << '\n';
    fprintf(RH_interp_ShockVelocity, "%s %f %s","NEW Shock Velocity WSh =",(*WSh)(0,POINT,ISH), " ");
    fprintf(RH_interp_ShockVelocity, "%f %s %f %s",(*WSh)(1,POINT,ISH), " ", (*WSh)(2,POINT,ISH), "\n");

  }

  //----------------------------------------------------------------------------//

void ComputeStateDps4Pg3DAssignUpstream::recoverDownOlDState(unsigned IV, unsigned ISH){

  // density
    xd.at(0) = (*ZRoeShdOld)(0,IV,ISH)/((*ZRoeShdOld)(4,IV,ISH)*ReferenceInfo::getRgas()); // perfect gas law

  // pressure
    xd.at(1)=(*ZRoeShdOld)(0,IV,ISH);

  // normal velocity
    xd.at(2) = (*ZRoeShdOld)(1,IV,ISH) * dx + (*ZRoeShdOld)(2,IV,ISH) * dy + (*ZRoeShdOld)(3,IV,ISH) * dz;

  // tangential velocity
    std::vector<double> u_tan(3);
    std::vector<double> u_norm(3);

    u_norm[0]=(*ZRoeShdOld)(1,IV,ISH)* dx;
    u_norm[1]=(*ZRoeShdOld)(2,IV,ISH)* dy;
    u_norm[2]=(*ZRoeShdOld)(3,IV,ISH)* dz;


    u_tan[0]=(*ZRoeShdOld)(1,IV,ISH)-xd.at(2)*dx;
    u_tan[1]=(*ZRoeShdOld)(2,IV,ISH)-xd.at(2)*dy;
    u_tan[2]=(*ZRoeShdOld)(3,IV,ISH)-xd.at(2)*dz;

  // std::cout << "xd.at(2)= u_norm_mod " << xd.at(2) << '\n';
  // std::cout << u_tan[0]<< " " <<u_tan[1] << " " <<u_tan[2] << "\n";


    xd.at(3)=sqrt(u_tan[0]*u_tan[0]+u_tan[1]*u_tan[1]+u_tan[2]*u_tan[2]);

    // std::cout << "xd.at(3)= u_tan_mod " << xd.at(3) << '\n';


    dxt=u_tan[0]/xd.at(3);
    dyt=u_tan[1]/xd.at(3);
    dzt=u_tan[2]/xd.at(3);

    R2 = sqrt(PhysicsInfo::getGam()*xd.at(1)/xd.at(0)) +
         0.5 * (PhysicsInfo::getGam()-1) * xd.at(2);






}



void ComputeStateDps4Pg3DAssignUpstream::recoverDownState(unsigned IV, unsigned ISH)
{

  // at the moment, Zrhoe is still in p,u,v,w,T format

  // xd.at(3) = -(*ZRoeShd)(2,IV,ISH) * dy + (*ZRoeShd)(3,IV,ISH) * dx;
  // xd.at(3) = xd.at(3)/(*ZRoeShd)(0,IV,ISH); // tangential
  // xd.at(2) = (*ZRoeShd)(2,IV,ISH) * dx + (*ZRoeShd)(3,IV,ISH) * dy;
  // xd.at(2) = xd.at(2)/(*ZRoeShd)(0,IV,ISH); // normal
  // xd.at(0) = (*ZRoeShd)(0,IV,ISH)*(*ZRoeShd)(0,IV,ISH); // density
  // help = pow((*ZRoeShd)(2,IV,ISH),2) + pow((*ZRoeShd)(3,IV,ISH),2);
  // xd.at(1) = (PhysicsInfo::getGam()-1)/PhysicsInfo::getGam() *
  //    ((*ZRoeShd)(0,IV,ISH) * (*ZRoeShd)(1,IV,ISH) - 0.5 * help); // pressure
  // R2 = sqrt(PhysicsInfo::getGam()*xd.at(1)/xd.at(0)) +
  //      0.5 * (PhysicsInfo::getGam()-1) * xd.at(2);

// density
  xd.at(0) = (*ZRoeShd)(0,IV,ISH)/((*ZRoeShd)(4,IV,ISH)*ReferenceInfo::getRgas()); // perfect gas law

// pressure
  xd.at(1)=(*ZRoeShd)(0,IV,ISH);

// normal velocity
  xd.at(2) = (*ZRoeShd)(1,IV,ISH) * dx + (*ZRoeShd)(2,IV,ISH) * dy + (*ZRoeShd)(3,IV,ISH) * dz;

// tangential velocity
  std::vector<double> u_tan(3);
  std::vector<double> u_norm(3);

  u_norm[0]=(*ZRoeShd)(1,IV,ISH)* dx;
  u_norm[1]=(*ZRoeShd)(2,IV,ISH)* dy;
  u_norm[2]=(*ZRoeShd)(3,IV,ISH)* dz;


  u_tan[0]=(*ZRoeShd)(1,IV,ISH)-xd.at(2)*dx;
  u_tan[1]=(*ZRoeShd)(2,IV,ISH)-xd.at(2)*dy;
  u_tan[2]=(*ZRoeShd)(3,IV,ISH)-xd.at(2)*dz;

// std::cout << "xd.at(2)= u_norm_mod " << xd.at(2) << '\n';
// std::cout << u_tan[0]<< " " <<u_tan[1] << " " <<u_tan[2] << "\n";


  xd.at(3)=sqrt(u_tan[0]*u_tan[0]+u_tan[1]*u_tan[1]+u_tan[2]*u_tan[2]);

  // std::cout << "xd.at(3)= u_tan_mod " << xd.at(3) << '\n';


  dxt=u_tan[0]/xd.at(3);
  dyt=u_tan[1]/xd.at(3);
  dzt=u_tan[2]/xd.at(3);

  R2 = sqrt(PhysicsInfo::getGam()*xd.at(1)/xd.at(0)) +
       0.5 * (PhysicsInfo::getGam()-1) * xd.at(2);

  logfile("Zd(1) ",(*ZRoeShd)(0,IV,ISH),"\n");
  logfile("Zd(2) ",(*ZRoeShd)(1,IV,ISH),"\n");
  logfile("Zd(3) ",(*ZRoeShd)(2,IV,ISH),"\n");
  logfile("Zd(4) ",(*ZRoeShd)(3,IV,ISH),"\n");
  logfile("Zd(5) ",(*ZRoeShd)(4,IV,ISH),"\n");

  logfile("xd.at(0)= rho ",xd.at(0),"\n");
  logfile("xd.at(1)= p ",xd.at(1),"\n");
  logfile("xd.at(2)= u_n ",xd.at(2),"\n");
  logfile("xd.at(3)= u_t ",xd.at(3),"\n");




  logfile("dx ",dx,"\n");
  logfile("dy ",dy,"\n");
  logfile("dz ",dz,"\n");
  logfile("dxt ",dxt,"\n");
  logfile("dyt ",dyt,"\n");
  logfile("dzt ",dzt,"\n");






}


//----------------------------------------------------------------------------//

void ComputeStateDps4Pg3DAssignUpstream::recoverUpState(unsigned IV, unsigned ISH)
{
  // xu.at(3) = -(*ZRoeShu)(2,IV,ISH) * dy + (*ZRoeShu)(3,IV,ISH) * dx;
  // xu.at(3) = xu.at(3)/(*ZRoeShu)(0,IV,ISH); // tangential
  // xu.at(2) = (*ZRoeShu)(2,IV,ISH) * dx + (*ZRoeShu)(3,IV,ISH) * dy;
  // xu.at(2) = xu.at(2)/(*ZRoeShu)(0,IV,ISH); // normal
  // xu.at(0) = (*ZRoeShu)(0,IV,ISH)*(*ZRoeShu)(0,IV,ISH); // density
  // help = pow((*ZRoeShu)(2,IV,ISH),2) + pow((*ZRoeShu)(3,IV,ISH),2);
  // xu.at(1) = (PhysicsInfo::getGam()-1)/PhysicsInfo::getGam() *
  //    ((*ZRoeShu)(0,IV,ISH) * (*ZRoeShu)(1,IV,ISH) - 0.5 * help); // pressure



    (*ZRoeShu)(1,IV,ISH)=ReferenceInfo::geturef();
    (*ZRoeShu)(2,IV,ISH)=0;
    (*ZRoeShu)(3,IV,ISH)=0;
  // density
    xu.at(0) = (ReferenceInfo::getpref())/(ReferenceInfo::getTref()*ReferenceInfo::getRgas()); // perfect gas law

  // pressure
    xu.at(1)=ReferenceInfo::getpref();

  // normal velocity
    xu.at(2) = (*ZRoeShu)(1,IV,ISH) * dx + (*ZRoeShu)(2,IV,ISH) * dy + (*ZRoeShu)(3,IV,ISH) * dz;

  // tangential velocity
    std::vector<double> u_tan(3);
    std::vector<double> u_norm(3);

    u_norm[0]=(*ZRoeShu)(1,IV,ISH)* dx;
    u_norm[1]=(*ZRoeShu)(2,IV,ISH)* dy;
    u_norm[2]=(*ZRoeShu)(3,IV,ISH)* dz;


    u_tan[0]=(*ZRoeShu)(1,IV,ISH)-xu.at(2)*dx;
    u_tan[1]=(*ZRoeShu)(2,IV,ISH)-xu.at(2)*dy;
    u_tan[2]=(*ZRoeShu)(3,IV,ISH)-xu.at(2)*dz;

    xu.at(3)=sqrt(u_tan[0]*u_tan[0]+u_tan[1]*u_tan[1]+u_tan[2]*u_tan[2]);

    dxt=u_tan[0]/xu.at(3);
    dyt=u_tan[1]/xu.at(3);
    dzt=u_tan[2]/xu.at(3);


  logfile("Zu(1) ",(*ZRoeShu)(0,IV,ISH),"\n");
  logfile("Zu(2) ",(*ZRoeShu)(1,IV,ISH),"\n");
  logfile("Zu(3) ",(*ZRoeShu)(2,IV,ISH),"\n");
  logfile("Zu(4) ",(*ZRoeShu)(3,IV,ISH),"\n");
  logfile("Zu(5) ",(*ZRoeShu)(4,IV,ISH),"\n");

  logfile("xu.at(0)= rho ",xu.at(0),"\n");
  logfile("xu.at(1)= p ",xu.at(1),"\n");
  logfile("xu.at(2)= u_n ",xu.at(2),"\n");
  logfile("xu.at(3)= u_t ",xu.at(3),"\n");

  logfile("dx ",dx,"\n");
  logfile("dy ",dy,"\n");
  logfile("dz ",dz,"\n");
  logfile("dxt ",dxt,"\n");
  logfile("dyt ",dyt,"\n");
  logfile("dzt ",dzt,"\n");


}

//----------------------------------------------------------------------------//

void ComputeStateDps4Pg3DAssignUpstream::computeDownState(unsigned IV, unsigned ISH)
{

  double z1, z2, z3, z4, z5;

  // pressure
  z1 = xd.at(1);

  // u component
  z2=xd.at(2)*dx+xd.at(3)*dxt;

  // v component
  z3=xd.at(2)*dy+xd.at(3)*dyt;

  // w component
  z4=xd.at(2)*dz+xd.at(3)*dzt;

  // T
  z5=xd.at(1)/(xd.at(0)*ReferenceInfo::getRgas()); //perfect gas law


  (*ZRoeShd)(0,IV,ISH) = z1;    // pressure
  (*ZRoeShd)(1,IV,ISH) = z2;    // u
  (*ZRoeShd)(2,IV,ISH) = z3;    // v
  (*ZRoeShd)(3,IV,ISH) = z4;    // w
  (*ZRoeShd)(4,IV,ISH) = z5;    // T

  // for (unsigned IA=0;IA<5;IA++){
  //   if((*ZRoeShd)(IA,IV,ISH)!=(*ZRoeShd)(IA,IV,ISH)){
  //     std::cout << "COMPUTE shock_point n. " << IV <<'\n';
  //   }
  // }
}

//----------------------------------------------------------------------------//

void ComputeStateDps4Pg3DAssignUpstream::computeUpState(unsigned IV, unsigned ISH)
{
  double z1, z2, z3, z4, z5;

  // pressure
  z1 = xu.at(1);

  // u component
  z2=xu.at(2)*dx+xu.at(3)*dxt;

  // v component
  z3=xu.at(2)*dy+xu.at(3)*dyt;

  // w component
  z4=xu.at(2)*dz+xu.at(3)*dzt;

  // T
  z5=xu.at(1)/(xu.at(0)*ReferenceInfo::getRgas()); //perfect gas law


  (*ZRoeShu)(0,IV,ISH) = z1;
  (*ZRoeShu)(1,IV,ISH) = z2;
  (*ZRoeShu)(2,IV,ISH) = z3;
  (*ZRoeShu)(3,IV,ISH) = z4;
  (*ZRoeShu)(4,IV,ISH) = z5;

}

//----------------------------------------------------------------------------//

void ComputeStateDps4Pg3DAssignUpstream::saveDownState(unsigned IV, unsigned ISH)
{
  for(unsigned I=0; I<(*ndof); I++) {
   (*ZRoeShdOld)(I,IV,ISH) = (*ZRoeShd)(I,IV,ISH);
   if ((*ZRoeShd)(I,IV,ISH)!=(*ZRoeShd)(I,IV,ISH)){std::cout << "saveDownState" << '\n';}
  }
}

//----------------------------------------------------------------------------//

void ComputeStateDps4Pg3DAssignUpstream::setAddress()
{
  unsigned start;
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

//----------------------------------------------------------------------------//

void ComputeStateDps4Pg3DAssignUpstream::setDiscSpeedSize()
{
  WSh->resize(PhysicsInfo::getnbDim(),
              nShockPoints->at(0),
              PhysicsInfo::getnbShMax());
}

//----------------------------------------------------------------------------//

void ComputeStateDps4Pg3DAssignUpstream::freeArray()
{
  delete ZRoeShu; delete ZRoeShd;
}

//----------------------------------------------------------------------------//

void ComputeStateDps4Pg3DAssignUpstream::setMeshData()
{
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector <double> > ("ZROE");
  nShockFaces =
        PhysicsData::getInstance().getData <vector<unsigned> > ("nShockFaces");
}

//----------------------------------------------------------------------------//

void ComputeStateDps4Pg3DAssignUpstream::setPhysicsData()
{
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
     PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  typeSh = PhysicsData::getInstance().getData <vector <string> > ("TYPESH");
  vShNor = PhysicsData::getInstance().getData <Array3D<double> > ("VSHNOR");
  WSh = PhysicsData::getInstance().getData <Array3D<double> > ("WSH");
  XYZSh = PhysicsData::getInstance().getData <Array3D<double> > ("XYZSH");
  ZRoeShuOld =
       PhysicsData::getInstance().getData <Array3D <double> > ("ZROESHuOLD");
  ZRoeShdOld =
       PhysicsData::getInstance().getData <Array3D <double> > ("ZROESHdOLD");
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting
