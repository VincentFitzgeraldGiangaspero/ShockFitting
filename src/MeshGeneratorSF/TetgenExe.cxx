// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cstdlib>
#include "MeshGeneratorSF/TetgenExe.hh"
#include "Framework/Log.hh"
#include "Framework/MeshData.hh"
#include "SConfig/ObjectProvider.hh"

//--------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<TetgenExe, MeshGenerator> TetgenExeProv("TetgenExe");

//--------------------------------------------------------------------------//

TetgenExe::TetgenExe(const std::string& objectName) :
  MeshGenerator(objectName)
{
}

//--------------------------------------------------------------------------//

TetgenExe::~TetgenExe()
{
}

//--------------------------------------------------------------------------//

void TetgenExe::setup()
{
  LogToScreen(VERBOSE, "TetgenExe::setup() => start\n");

  LogToScreen(VERBOSE, "TetgenExe::setup() => end\n");
}

//--------------------------------------------------------------------------//

void TetgenExe::unsetup()
{
  LogToScreen(VERBOSE, "TetgenExe::unsetup()\n");
}

//--------------------------------------------------------------------------//

void TetgenExe::generate()
{
  LogToScreen(INFO,"TetgenExe::generate()\n");


  fname = MeshData::getInstance().getData <stringstream>("FNAME");

  command = "../../MeshGeneratorSF/Tetgen/tetgen -pYn "
            + fname->str() +".poly" + " > log/TetgenExe.log";

  std::cout << "-command  => " << command << "\n";

  system(command.c_str());

  if(system(command.c_str())!=0) {
   cout << "TetgenExe::error => Tetgen Mesh Generator execution failed\n";
   exit(1); }
}

//--------------------------------------------------------------------------//

void TetgenExe::generate(string processingFile)
{
  LogToScreen(INFO,"TetgenExe::generate()\n");

  command = "../../MeshGeneratorSF/Tetgen/tetgen -nefp "
            + processingFile +  " > log/TetgenExe.log";

  std::cout << "command " << command << "\n";

  system(command.c_str());

  if(system(command.c_str())!=0) {
   cout << "TetgenExe::error => Tetgen Mesh Generator execution failed\n";
   exit(1); }
}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
