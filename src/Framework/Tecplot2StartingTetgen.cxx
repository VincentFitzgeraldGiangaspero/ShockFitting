// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// modifiche vincent 29/9

#include <stdio.h>
#include <iomanip>
#include "ConverterSF/Tecplot2StartingTetgen.hh"
#include "Framework/Log.hh"
#include "Framework/PhysicsInfo.hh"
#include "MathTools/Jcycl.hh"
#include "SConfig/ObjectProvider.hh"
#include "SConfig/Factory.hh"
#include "SConfig/ConfigFileReader.hh"

//----------------------------------------------------------------------------//

using namespace std;
using namespace SConfig;

//----------------------------------------------------------------------------//

namespace ShockFitting {

//----------------------------------------------------------------------------//

// this variable instantiation activates the self-registration mechanism
ObjectProvider<Tecplot2StartingTetgen, Converter>
Tecplot2StartingTetgenProv("Tecplot2StartingTetgen");

//----------------------------------------------------------------------------//

Tecplot2StartingTetgen::Tecplot2StartingTetgen(const std::string& objectName) :
  Converter(objectName)
{
  // first the *.plt and then the *surf.plt
  m_meshInputfile = vector<string>(); 
  addOption("InputFile",&m_meshInputfile,
            "Tecplot files containing the captured solution");
  // it must be set equal to 0 per Perfect Gas model
  m_nbSpecies = 0;
  addOption("nbSpecies",&m_nbSpecies,
            "Specifies the number of chemical species");
  m_tecplotExtraValues = false;
  addOption("extraValuesPrinted",&m_tecplotExtraValues,
            "Specifies if extra values are printed in the tecplot file");

  m_prim2param.name() = "dummyVariableTransformer";
}

//----------------------------------------------------------------------------//

Tecplot2StartingTetgen::~Tecplot2StartingTetgen()
{
}

//----------------------------------------------------------------------------//

void Tecplot2StartingTetgen::setup()
{
  LogToScreen(VERBOSE, "Tecplot2StartingTetgen::setup() => start\n");

  m_prim2param.ptr()->setup();

  LogToScreen(VERBOSE, "Tecplot2StartingTetgen::setup() => end\n");
}

//----------------------------------------------------------------------------//

void Tecplot2StartingTetgen::unsetup()
{
  LogToScreen(VERBOSE, "Tecplot2StartingTetgen::unsetup()\n");

  m_prim2param.ptr()->unsetup();
}

//----------------------------------------------------------------------------//

void Tecplot2StartingTetgen::configure(OptionMap& cmap, const std::string& prefix)
{
  Converter::configure(cmap, prefix);

  // assign strings on input.case file to variable transformer object
  m_prim2param.name() = m_inFmt+"2"+m_outFmt+m_modelTransf+m_additionalInfo;

  if (ConfigFileReader::isFirstConfig()) {
   m_prim2param.ptr().reset(SConfig::Factory<VariableTransformer>::getInstance().
                            getProvider(m_prim2param.name())
                            ->create(m_prim2param.name()));
  }

  // configure variable transformer object
  configureDeps(cmap, m_prim2param.ptr().get());

}

//----------------------------------------------------------------------------//

void Tecplot2StartingTetgen::convert()
{
  LogToScreen (INFO, "Tecplot2StartingTetgen::convert()\n");

  // read Tecplot format file
  LogToScreen(DEBUG_MIN, "Tecplot2StartingTetgen::reading Tecplot format\n");
  readTecplotFmt();

  // make the transformation from primitive variables to
  // Roe parameter vector variables

// questo é da sistemare

/*  
vector <double> m_prim(ndof);
  vector <double> m_zroe(ndof);
  vector <double> m_XYZ(PhysicsInfo::getnbDim());
  for(unsigned IPOIN=0; IPOIN<npoin; IPOIN++) {
   for(unsigned i=0; i<ndof; i++) 
    { m_prim.at(i) = prim.at(IPOIN*ndof+i); }
   for(unsigned i=0; i<PhysicsInfo::getnbDim(); i++)
    { m_XYZ.at(i) = XYZ.at(IPOIN*PhysicsInfo::getnbDim()+i); }
   m_prim2param.ptr()->transform(m_prim,m_XYZ,m_zroe);
   for(unsigned i=0; i<ndof; i++) { zroe.at(IPOIN*ndof+i) = m_zroe.at(i); }
  }

*/

  // write Tetgen format files
  LogToScreen(DEBUG_MIN, "Tecplot2StartingTetgen::writing Tetgen format\n");
  writeTetgenFmt(); 
}

//----------------------------------------------------------------------------//

void Tecplot2StartingTetgen::readTecplotFmt()
{
  // dummy variables
  string dumstring; double dumvalue;
  string::size_type i;
  vector <unsigned> m_bndfac(3);	//era 2
  unsigned countbnd;

  // number of point for each coloured face
  unsigned npoinFac;

  // map array used to match tecplot id-boundary edges
  // and id-boundary edges
  // it will be resized after
  Array2D <unsigned> array_bnd(1,1); // Array2D <unsigned> array_bnd(1,1);

  // coordinates of the boundary points
  vector <double> XYZ_bnd(3);	//era 2

  // reading file
  ifstream file;

  ndof = m_nbSpecies+5; // ho messo 5 perché é 3D erano 4 prima


  file.open(string(m_meshInputfile.at(0)).c_str());


LogToScreen(INFO, "file plt aperto\n");


  // read TITLE = ...
  getline(file,dumstring);
  // read VARIABLES = x x x ...
  getline(file,dumstring);
 
  // read ZONE   T="ZONE0 TETRA",
  file >> dumstring >> dumstring >> dumstring >> dumstring;

  // read N=....
  file >> dumstring;


  // take the number of points from the string
  i = dumstring.find(string("N="));
  if (i != std::string::npos) { dumstring.erase(i,string("N=").length()); }
  i = dumstring.find(string(","));
  if (i != std::string::npos) { dumstring.erase(i,string(",").length()); }

  npoin = atoi(dumstring.c_str());

  // read E=....
  file >> dumstring;
  
  // take the number of elements from the string
  i = dumstring.find(string("E="));
  if (i != std::string::npos) { dumstring.erase(i,string("E=").length()); }
   i = dumstring.find(string(","));
  if (i != std::string::npos) { dumstring.erase(i,string(",").length()); }

  nelem = atoi(dumstring.c_str());


  // read F=FEPOINT, ET=TETRAHEDRON, SOLUTIONTIME=0
  getline(file,dumstring);

  prim.resize(ndof*npoin);
  zroe.resize(ndof*npoin);
  XYZ.resize(PhysicsInfo::getnbDim()*npoin);
  nodcod.resize(npoin);


  for(unsigned IPOIN=0;IPOIN<npoin;IPOIN++) {
   // read the nodal coordinates
   for(unsigned IV=0; IV<PhysicsInfo::getnbDim(); IV++)
    { file >> XYZ.at(IPOIN*PhysicsInfo::getnbDim()+IV); }
   // read the nodal states
   for(unsigned I=0; I<ndof; I++) { file >> prim.at(IPOIN*ndof+I); }
   if(m_tecplotExtraValues) {
    for(unsigned I=0; I<4; I++) { file >> dumvalue; }
   }
  } 

  // the list of elements will not be read, celnod and celcel will
  // be generated calling Tetgen
 
  file.close(); // close *.plt


LogToScreen(INFO, "file plt chiuso\n");


  file.open(string(m_meshInputfile.at(1)).c_str()); // open *surf.plt 


  // read TITLE = ...
  getline(file,dumstring);
  // read VARIABLES = x x x ...
  getline(file,dumstring);



LogToScreen(INFO, "file surf.plt aperto\n");

  unsigned IFACE=0;
  nbfac = 0;
  while(!file.eof()) {
   
   // read ZONE
   file >> dumstring;
  
  

// i miei file hanno prima la T poi la N quindi devo scambiare le due parti

   if(dumstring=="ZONE" || dumstring=="ZONE " || dumstring=="\nZONE" || dumstring=="ZONE  ") {

// ZONE   T= "ZONE0 2", N=1018,

    // read "T=.."
    file >> dumstring;



    i = dumstring.find(string("T=\""));
    if (i != std::string::npos) { dumstring.erase(i,string("T=\"").length()); }


 //   i = dumstring.find(string(","));
 //   if (i != std::string::npos) { dumstring.erase(i,string(",").length()); }

    namebnd.resize(IFACE+1);


// read ZONE0

    file >> dumstring; 



    i = dumstring.find(string("\"ZONE0"));
    if (i != std::string::npos) { dumstring.erase(i,string("\"ZONE0").length()); }

file >> dumstring; // aggiunto

    i = dumstring.find(string("\","));
    if (i != std::string::npos) { dumstring.erase(i,string("\",").length()); }

    namebnd.at(IFACE)=atoi(dumstring.c_str());



   

    // read "N=.."
    file >> dumstring;

    i = dumstring.find(string("N="));
    if (i != std::string::npos) { dumstring.erase(i,string("N=").length()); }
    i = dumstring.find(string(","));
    if (i != std::string::npos) { dumstring.erase(i,string(",").length()); }

    npoinFac = atoi(dumstring.c_str()); // QUI LEGGE I Numeri di punti di ogni zona
	

    array_bnd.resize(2,npoinFac);	// 


    // read "E=.."
    file >> dumstring;

    i = dumstring.find(string("E="));
    if (i != std::string::npos) { dumstring.erase(i,string("E=").length()); }
    i = dumstring.find(string(","));
    if (i != std::string::npos) { dumstring.erase(i,string(",").length()); }

    nFacB.resize(IFACE+1);
    nFacB.at(IFACE) = atoi(dumstring.c_str());

cout << "nFacB.size " << nFacB.size() << " \n"; 
cout << "npoinFac.size " << npoinFac<< " \n"; 

cout << "nFacB " << nFacB.at(IFACE) << " \n"; 


    bndfac.resize((bndfac.size()+4*(nFacB.at(IFACE)))); // qua ci va il 4 perché ho 3 vertici piu il marker 

cout << "bndfac.size " << bndfac.size() << " \n"; 

    // read F=FEPOINT, ET=Tetgen, SOLUTIONTIME=0
    file >> dumstring  >> dumstring  >> dumstring;

cout << "dumstring " << dumstring << " \n"; 

LogToScreen(INFO, "assegnazione punti\n");

    // read boundary nodes values
    // and assign the correct IDnumber to the boundary points
    countbnd=0;
cout << "npoinFac.size " << npoinFac<< " \n";
    for(unsigned IBPOIN=0;IBPOIN<npoinFac;IBPOIN++) {
//cout << "punto " <<; 
     for(unsigned IV=0;IV<PhysicsInfo::getnbDim();IV++) { 
      file >> XYZ_bnd.at(IV);
     }


     for(unsigned IPOIN=0;IPOIN<npoin;IPOIN++) {
      if(XYZ_bnd.at(0)==XYZ.at(IPOIN*PhysicsInfo::getnbDim()+0) &&
         XYZ_bnd.at(1)==XYZ.at(IPOIN*PhysicsInfo::getnbDim()+1) &&
	 XYZ_bnd.at(2)==XYZ.at(IPOIN*PhysicsInfo::getnbDim()+2)) {
       array_bnd(0,IBPOIN)=IPOIN+1;
       array_bnd(1,IBPOIN)=IBPOIN+1;
       countbnd++;
      }
     }


     for(unsigned IV=0;IV<ndof;IV++) { file >> dumvalue; }
     if(m_tecplotExtraValues) {
      for(unsigned IV=0; IV<4; IV++) { file >> dumvalue; }
     }
    }


cout << " countbound map XYZ " << countbnd << "\n";


    if(countbnd!=npoinFac) {
     cout << "Tecplot2StartingTetgen::error => no match between\n";
     cout << "  coordinates of boundary points\n";
     cout << "countbnd= "<<countbnd<<" # boundary points= "<<npoinFac<<endl;
     exit(1);
    }
 
    // read geometry entities list that is boundary data
    countbnd=0;
    for(unsigned IB=nbfac;IB<nFacB.at(IFACE)+nbfac;IB++) {
     for(int IV=0;IV<3;IV++) { 					// for(int IV=0;IV<2;IV++) {
      file >> m_bndfac.at(IV);
      for(unsigned IBPOIN=0;IBPOIN<npoinFac;IBPOIN++) {
       if(m_bndfac.at(IV)==array_bnd(1,IBPOIN)) {
        bndfac.at(IB*4+IV)=array_bnd(0,IBPOIN); countbnd++; break; }
      } // for IBPOIN
     } // for IV=1:ndim
     bndfac.at(IB*4+3) = namebnd.at(IFACE);
    } 

cout << " countbound map bndfac " << countbnd/3 << "\n";

    if(countbnd!=nFacB.at(IFACE)*3) {   // qui ci va il *3 perche le facce sono fatte da 3 vertici
     cout << "Tecplot2StartingTetgen::error => wrong map for\n";
     cout << "tecplot ID boundary edges and ID boundary edges\n";
     cout << "countbnd= "<<countbnd<<" # edge points= "<<nFacB.at(IFACE)*3<<endl;
     exit(1);
    }

    // update the total number of boundary edges ( FACES )
    nbfac = nFacB.at(IFACE) + nbfac;

    IFACE++;
   } // if dumstring=ZONE
  } // while !eof

  NCLR = IFACE;
/*nbfac=0;
for(unsigned IF=0;IF<NCLR;IF++) {
for(unsigned IB=nbfac;IB<nFacB.at(IF)+nbfac;IB++) {
cout << "IB " << IB << endl;
for(int IV=0;IV<3;IV++) {
cout << bndfac.at(IB*3+IV) << " ";
}
cout << endl;
}
nbfac = nFacB.at(IFACE) + nbfac;
cout << nbfac << endl;
}*/



  file.close(); // close *surf.plt

  // fill nodcod vector
  setNodcod(); // da controllare, mette il 2 in nodcod a tutti i boundary points
}

//----------------------------------------------------------------------------//

void Tecplot2StartingTetgen::setNodcod()
{
  int IPOIN;
/*
cout << " setnodcod " << "\n";
cout << " nodcod size " << nodcod.size() << " \n";
cout << " nbfac = " << nbfac << "\n";
cout << " bndfac.size " << bndfac.size() << " \n";
cout << " bndfac.at(0) " << bndfac.at(0) << " \n";

cout << " bndfac.at(1) " << bndfac.at(1) << " \n";

cout << " bndfac.at(2) " << bndfac.at(2) << " \n";


cout << " bndfac.at(3) " << bndfac.at(3) << " \n";

cout << " bndfac.at(4) " << bndfac.at(4) << " \n";

cout << " bndfac.at(5) " << bndfac.at(5) << " \n";
cout << " bndfac.at(6) " << bndfac.at(6) << " \n";
cout << " bndfac.at(7) " << bndfac.at(7) << " \n";
cout << " bndfac.at(8) " << bndfac.at(8) << " \n";
cout << " bndfac.at(9) " << bndfac.at(9) << " \n";
cout << " bndfac.at(last-1) " << bndfac.at(bndfac.size()-17) << " \n";
cout << " bndfac.at(last-1) " << bndfac.at(bndfac.size()-16) << " \n";
cout << " bndfac.at(last-1) " << bndfac.at(bndfac.size()-15) << " \n";
cout << " bndfac.at(last-1) " << bndfac.at(bndfac.size()-14) << " \n";
cout << " bndfac.at(last-1) " << bndfac.at(bndfac.size()-13) << " \n";
cout << " bndfac.at(last-1) " << bndfac.at(bndfac.size()-12) << " \n";
cout << " bndfac.at(last-1) " << bndfac.at(bndfac.size()-11) << " \n";
cout << " bndfac.at(last-1) " << bndfac.at(bndfac.size()-10) << " \n";
cout << " bndfac.at(last-1) " << bndfac.at(bndfac.size()-9) << " \n";
cout << " bndfac.at(last-1) " << bndfac.at(bndfac.size()-8) << " \n";
cout << " bndfac.at(last-1) " << bndfac.at(bndfac.size()-7) << " \n";
cout << " bndfac.at(last-1) " << bndfac.at(bndfac.size()-6) << " \n";
cout << " bndfac.at(last-1) " << bndfac.at(bndfac.size()-5) << " \n";
cout << " bndfac.at(last-1) " << bndfac.at(bndfac.size()-4) << " \n";
cout << " bndfac.at(last-1) " << bndfac.at(bndfac.size()-3) << " \n";

cout << " bndfac.at(last-1) " << bndfac.at(bndfac.size()-2) << " \n";

cout << " bndfac.at(last) " << bndfac.at(bndfac.size()-1) << " \n";

*/
int g=0;

  for(unsigned IB=0; IB<nbfac; IB++) {
   for(unsigned I=0; I<3; I++) {
    IPOIN = bndfac.at(IB*4+I);
	if (nodcod.at(IPOIN-1)==2){g++;}
    nodcod.at(IPOIN-1) = 2;
   }
  }

int conto=0;

for (unsigned k=0; k<npoin; k++){
if (nodcod.at(k)==2)
{conto++;}
}

cout << " nodcod.at(0) " << nodcod.at(0) << " \n";


cout << " nodcod size = " << nodcod.size() << "\n";

cout << " nodcod gia contati " << g << "\n";
cout << " nodcod = 2 " << conto << " \n" ;


}

//----------------------------------------------------------------------------//

void Tecplot2StartingTetgen::writeTetgenFmt()
{
  // dummystring used to open .node or .poly files
  string dummystring;

  // dummy iface boundary marker
  stringstream NBND;

  // writing file (Tetgen node file to be overwritten)
  FILE* Tetgenfile;

  // write on .node file
  dummystring = "na00.node";

  Tetgenfile = fopen(dummystring.c_str(),"w");

LogToScreen(INFO, "write node\n");


  fprintf(Tetgenfile,"%u %s %u",npoin," ",PhysicsInfo::getnbDim());
  fprintf(Tetgenfile,"%s %u %s"," ",ndof," 1\n");
  for(unsigned IPOIN=0; IPOIN<npoin; IPOIN++) {
   fprintf(Tetgenfile,"%u %s",IPOIN+1," ");
   for(unsigned IA=0; IA<PhysicsInfo::getnbDim(); IA++)
    { fprintf(Tetgenfile,"%20.16F %s",XYZ.at(IPOIN*PhysicsInfo::getnbDim()+IA)," "); }
   for(unsigned IA=0; IA<ndof; IA++)
   { fprintf(Tetgenfile,"%20.16F %s",prim.at(IPOIN*ndof+IA)," "); } // zroe.at era prima
   fprintf(Tetgenfile,"%u %s",nodcod.at(IPOIN),"\n");  	// se viene attivato il setnodcod sono 2 tutti boundary points
  }

int a=0;

for (unsigned IPOIN=0; IPOIN<npoin; IPOIN++) {
   if(nodcod.at(IPOIN)>0) {a++;} 

}



  fclose(Tetgenfile);

cout << " npoin " << npoin << " \n";

cout << " nodcod filled " << a << " \n";
 
  // write on .poly file
  dummystring = "na00.poly";

  Tetgenfile = fopen(dummystring.c_str(),"w");

LogToScreen(INFO, "write poly\n");


  fprintf(Tetgenfile,"%s %u %s", "0 ", PhysicsInfo::getnbDim(), " 0 1\n");
  fprintf(Tetgenfile,"%u %s",nbfac," 1\n");

  unsigned m_nbfac=0;
  for(unsigned IFACE=0;IFACE<NCLR;IFACE++)  {
   for(unsigned IB=m_nbfac;IB<nFacB.at(IFACE)+m_nbfac;IB++) {
    NBND.str(string());
    NBND << namebnd.at(IFACE); // c++ indeces start from 0
    if(NBND.str() == "InnerSup" || NBND.str() == "InnerSub") { NBND << "10"; }
    fprintf(Tetgenfile,"%5u %5u %3s %s",1,0,NBND.str().c_str()," \n");
    fprintf(Tetgenfile,"%5u",3);		// fprintf(Tetgenfile,"%5u",IB+1);
    fprintf(Tetgenfile,"%5i %5i %5i %s",bndfac.at(IB*4+0),bndfac.at(IB*4+1),bndfac.at(IB*4+2)," \n"); // DAAGGIUNGERE
   
   }
   m_nbfac=m_nbfac+nFacB.at(IFACE);
  }

  fprintf(Tetgenfile,"%s","0\n"); // write number of holes   
  fclose(Tetgenfile);
}

//----------------------------------------------------------------------------//

} // namespace ShockFitting

