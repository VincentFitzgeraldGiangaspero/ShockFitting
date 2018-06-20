// Copyright (C) 2014 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// vincent 31/01/2017

#include "ShockDetectorSF/AssignShockPointState3D.hh"

#include "Framework/MeshData.hh"
#include "Framework/PhysicsData.hh"
#include "Framework/PhysicsInfo.hh"
#include "Framework/ReferenceInfo.hh"

#include <cmath>
#include "Framework/FileLogManip.hh"
#include "MathTools/Array3D.hh"


//--------------------------------------------------------------------------//

using namespace std;

//--------------------------------------------------------------------------//

namespace ShockFitting {

//--------------------------------------------------------------------------//

AssignShockPointState3D::AssignShockPointState3D()
{
}

//--------------------------------------------------------------------------//

AssignShockPointState3D::~AssignShockPointState3D()
{
}

//--------------------------------------------------------------------------//

void AssignShockPointState3D::setup(double shLayerThick)
{
  LogToScreen(INFO, "AssignShockPointState3D::setup()\n");


  setMeshData();
  setPhysicsData();

  // assign starting pointers to array
  setAddress();

  // resize vectors and array
  setSize();

  // set the distance used to extract the upstream and downstream coordinates
  setShockLayerThickness(shLayerThick);

  logfile.Open(getClassName());
}

//--------------------------------------------------------------------------//

void AssignShockPointState3D::unsetup()
{

  logfile.Close();

  // de-allocate dynamic array
  freeArray();
}

//--------------------------------------------------------------------------//

void AssignShockPointState3D::extractDownstreamAndUpstreamPoints(Array3D<double> Normals_Array)
{
  cout << "     => AssignShockPointState3D::extractDownstreamAndUpstreamPoints()\n";

  m_shockLayerThickness *= MeshData::getInstance().getDXCELL();  //

  for (unsigned ISHPOIN=0; ISHPOIN<nShockPoints->at(0); ISHPOIN++){
    for (unsigned IV=0; IV<3; IV++){
      normals(IV,ISHPOIN,0)=Normals_Array(IV,ISHPOIN,0);
    }
  }


  cout << "-m_shockLayerThickness = " << m_shockLayerThickness << "\n";
  cout << "-MeshData::getInstance().getDXCELL() = " << MeshData::getInstance().getDXCELL() << "\n";


  vector<double> X(npoin->at(0));       // aggiunte le Z
  vector<double> Y(npoin->at(0));
  vector<double> Z(npoin->at(0));
  for(unsigned IPOIN=0;IPOIN<npoin->at(0);IPOIN++) {
    X.at(IPOIN) = (*XYZ)(0,IPOIN);
    Y.at(IPOIN) = (*XYZ)(1,IPOIN);
    Z.at(IPOIN) = (*XYZ)(2,IPOIN);
  }

  const double minDomainXCoordinate = *min_element(X.begin(),X.end());
  const double maxDomainXCoordinate = *max_element(X.begin(),X.end());
  const double minDomainYCoordinate = *min_element(Y.begin(),Y.end());
  const double maxDomainYCoordinate = *max_element(Y.begin(),Y.end());
  const double minDomainZCoordinate = *min_element(Z.begin(),Z.end());
  const double maxDomainZCoordinate = *max_element(Z.begin(),Z.end());


  double m_Xu, m_Yu;
  double m_Xd, m_Yd;
  double m_Zu, m_Zd;

  double counter_out=0;

  for(unsigned ISH=0;ISH<(*nShocks);ISH++) {
    unsigned m_shPoint = 0;

    for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
      m_Xu = (*XYZSh)(0,ISHPOIN,ISH) +
      normals(0,ISHPOIN,ISH) * m_shockLayerThickness;
      m_Yu = (*XYZSh)(1,ISHPOIN,ISH) +
      normals(1,ISHPOIN,ISH) * m_shockLayerThickness;
      m_Zu = (*XYZSh)(2,ISHPOIN,ISH) +
      normals(2,ISHPOIN,ISH) * m_shockLayerThickness;

      m_Xd = (*XYZSh)(0,ISHPOIN,ISH) -
      normals(0,ISHPOIN,ISH) * m_shockLayerThickness;
      m_Yd = (*XYZSh)(1,ISHPOIN,ISH) -
      normals(1,ISHPOIN,ISH) * m_shockLayerThickness;
      m_Zd = (*XYZSh)(2,ISHPOIN,ISH) -
      normals(2,ISHPOIN,ISH) * m_shockLayerThickness;

      // check if the extracted shock points are out of the domain
      // they are erased
      if(m_Xu>=minDomainXCoordinate && m_Xu<=maxDomainXCoordinate &&
        m_Yu>=minDomainYCoordinate && m_Yu<=maxDomainYCoordinate &&
        m_Zu>=minDomainZCoordinate && m_Zu<=maxDomainZCoordinate &&
        m_Xd>=minDomainXCoordinate && m_Xd<=maxDomainXCoordinate &&
        m_Yd>=minDomainYCoordinate && m_Yd<=maxDomainYCoordinate &&
        m_Zd>=minDomainZCoordinate && m_Zd<=maxDomainZCoordinate  )
        {
          XYZShu(0,m_shPoint,ISH) = m_Xu; XYZShu(1,m_shPoint,ISH) = m_Yu; XYZShu(2,m_shPoint,ISH) = m_Zu;
          XYZShd(0,m_shPoint,ISH) = m_Xd; XYZShd(1,m_shPoint,ISH) = m_Yd; XYZShd(2,m_shPoint,ISH) = m_Zd;
          (*XYZSh)(0,m_shPoint,ISH) = (*XYZSh)(0,ISHPOIN,ISH);
          (*XYZSh)(1,m_shPoint,ISH) = (*XYZSh)(1,ISHPOIN,ISH);
          (*XYZSh)(2,m_shPoint,ISH) = (*XYZSh)(2,ISHPOIN,ISH);
          ++m_shPoint;
        }
        else {
          cout << "     => AssignShockPointState3D:: (!) warning => IDpoint ";
          cout << ISHPOIN+1 << " is out of the domain\n";
          cout << "\n" <<  "         up= " << m_Xu << "," << m_Yu << "," << m_Zu;
          cout << "; down= " << m_Xd << "," << m_Yd << "," << m_Zd;
          cout << " has been erased\n\n";
          counter_out++;
        }


      }
      // update the number of shock points because some of them could have been erased
      nShockPoints->at(ISH) = m_shPoint;
      // nShockEdges->at(ISH) = nShockPoints->at(ISH)-1;
    }


    // plot downstream and upstream coordinates
    // it is useful to check if they are extracted within the shock layer
    FILE* plotState;
    stringstream fileName;

    for(unsigned ISH=0;ISH<(*nShocks);ISH++) {
      fileName.str(string());
      fileName << "./log/UpstreamDownstreamCoor_ShNb_";
      fileName << ISH+1 << ".dat";
      plotState = fopen(fileName.str().c_str(),"w");

      fprintf(plotState,"%s","TITLE = Upstream and downstream coor\n");
      fprintf(plotState,"%s","VARIABLES = \"x0\" \"x1\" \"x2\"\n");
      fprintf(plotState,"%s","ZONE T = \"Up and Down coord\"\n");
      fprintf(plotState,"%s"," STRANDID=0, SOLUTIONTIME=0\n");
      fprintf(plotState,"%s %u %s","I= ",nShockPoints->at(ISH),", J=3, K=1, ZONETYPE=Ordered\n");
      fprintf(plotState,"%s","DATAPACKING=POINT\n");
      fprintf(plotState,"%s","DT = (SINGLE, SINGLE)\n");

      for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
        fprintf(plotState,"%22.14F",XYZShu(0,ISHPOIN,ISH));
        fprintf(plotState,"%22.14F",XYZShu(1,ISHPOIN,ISH));
        fprintf(plotState,"%22.14F",XYZShu(2,ISHPOIN,ISH));
        fprintf(plotState,"%s","\n");

      }

      for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {
        fprintf(plotState,"%22.14F",XYZShd(0,ISHPOIN,ISH));
        fprintf(plotState,"%22.14F",XYZShd(1,ISHPOIN,ISH));
        fprintf(plotState,"%22.14F",XYZShd(2,ISHPOIN,ISH));
        fprintf(plotState,"%s","\n");
      }

      fclose(plotState);
    }
  }

//--------------------------------------------------------------------------//

void AssignShockPointState3D::interpDownstreamAndUpstreamState()
{
  cout << "     => AssignShockPointState3D::interpDownstreamAndUpstreamState()\n";

  unsigned ISH=0;

  std::vector<double> u_tan_up(3);
  std::vector<double> u_norm_up(3);

  std::vector<double> u_tan_d(3);
  std::vector<double> u_norm_d(3);


  std::vector<double> uu(3),ud(3);

  double u_norm_up_mod;


  double Md,Mu;

  double R=ReferenceInfo::getRgas();
  double gamma=PhysicsInfo::getGam();
  double delta=0.2;

  double a,b;

  double u_norm_down_mod;

  double u_down_mod;

  const double pi = 3.1415926535897;

  double sigma;

  double u_tan_up_mod;


  double pd,Td,rhou;
  double pu,Tu,rhod;


  for (unsigned ISHPOIN=0; ISHPOIN < nShockPoints->at(0); ISHPOIN++){

    u_norm_up_mod=0;

    for (unsigned I=0; I<3; I++){
      u_tan_up[I]=0;
      u_norm_up[I]=0;
    }

    (*ZRoeShu)(0,ISHPOIN,ISH) = ReferenceInfo::getpref();    //p
    (*ZRoeShu)(4,ISHPOIN,ISH) = ReferenceInfo::getTref();    //T

    (*ZRoeShu)(1,ISHPOIN,ISH) = ReferenceInfo::geturef();   //u
    (*ZRoeShu)(2,ISHPOIN,ISH) = 0;    //v
    (*ZRoeShu)(3,ISHPOIN,ISH) = 0;    //w

    pu=(*ZRoeShu)(0,ISHPOIN,ISH);
    uu[0]=(*ZRoeShu)(1,ISHPOIN,ISH);
    uu[1]=(*ZRoeShu)(2,ISHPOIN,ISH);
    uu[2]=(*ZRoeShu)(3,ISHPOIN,ISH);
    Tu=(*ZRoeShu)(4,ISHPOIN,ISH);

    rhou=pu/(gamma*R*Tu);


    // Mu=(ReferenceInfo::geturef())/(PhysicsInfo::getGam()*Tu*ReferenceInfo::getRgas());

    Mu=ReferenceInfo::getMachref();

    // std::cout << "here" << '\n';

    for(unsigned IDOF=0;IDOF<(*ndof);IDOF++) { (*ZRoeShd)(IDOF,ISHPOIN,ISH) = 0;}

    //upstream
    u_norm_up_mod=uu[0]* normals(0,ISHPOIN,ISH)+uu[1]* normals(1,ISHPOIN,ISH)+uu[2]* normals(2,ISHPOIN,ISH);

    u_norm_up[0]=u_norm_up_mod* normals(0,ISHPOIN,ISH);
    u_norm_up[1]=u_norm_up_mod* normals(1,ISHPOIN,ISH);
    u_norm_up[2]=u_norm_up_mod* normals(2,ISHPOIN,ISH);


    u_tan_up[0]=uu[0]-u_norm_up_mod*normals(0,ISHPOIN,ISH);
    u_tan_up[1]=uu[1]-u_norm_up_mod*normals(1,ISHPOIN,ISH);
    u_tan_up[2]=uu[2]-u_norm_up_mod*normals(2,ISHPOIN,ISH);

    u_tan_up_mod=sqrt(u_tan_up[0]*u_tan_up[0]+u_tan_up[1]*u_tan_up[1]+u_tan_up[2]*u_tan_up[2]);

 
 
   // if (std::abs(((u_norm_up_mod/ReferenceInfo::geturef())+1))<pow(10.0,(-3))){ // normal shock


      Md=sqrt((0.2*pow(Mu,2)+1)/(gamma*pow(Mu,2)-0.2));
      pd=pu*(1/(gamma+1))*(2*gamma*pow(Mu,2)-gamma+1);
      Td=Tu*(1/pow((gamma+1),2))*(2*gamma*(pow(Mu,2))-gamma+1)*(gamma-1+2/(pow(Mu,2)));
      ud[0]=Md*sqrt(gamma*R*Td);
      ud[1]=0;
      ud[2]=0;
   // }

   // else {


      sigma=atan(u_norm_up_mod/u_tan_up_mod);

      a=(2/tan(sigma))*((pow(Mu,2))*pow(sin(sigma),2))-1;
      b=(pow(Mu,2))*(gamma+cos(2*sigma))+2;
      delta=atan(a/b);
      u_norm_down_mod=u_tan_up_mod*tan(sigma-delta);

      rhod=rhou*(u_norm_up_mod/u_norm_down_mod);

      pd=pu+rhou*u_norm_up_mod*(u_norm_up_mod-u_norm_down_mod);

      ud[0]=u_tan_up[0]+u_norm_down_mod*normals(0,ISHPOIN,ISH);
      ud[1]=u_tan_up[1]+u_norm_down_mod*normals(1,ISHPOIN,ISH);
      ud[2]=u_tan_up[2]+u_norm_down_mod*normals(2,ISHPOIN,ISH);

      Td=pd/(R*rhod);


      if((normals(0,ISHPOIN,ISH)<pow(10,-3))&&(normals(0,ISHPOIN,ISH)>-pow(10,-3))){

        std::cout << "if activated  " << ISHPOIN <<'\n';


        sigma=atan(u_norm_up_mod/u_tan_up_mod)/pi*180;
        delta=abs(sigma)-(pi/180)*abs(atan(u_norm_down_mod/u_tan_up_mod));
        u_down_mod=ud[0]*ud[0]+ud[1]*ud[1]+ud[2]*ud[2];
        Md=u_down_mod/sqrt(gamma*pd/rhod);
      }
   // }

    (*ZRoeShd)(0,ISHPOIN,ISH)=pd;
    (*ZRoeShd)(1,ISHPOIN,ISH)=ud[0];
    (*ZRoeShd)(2,ISHPOIN,ISH)=ud[1];
    (*ZRoeShd)(3,ISHPOIN,ISH)=ud[2];
    (*ZRoeShd)(4,ISHPOIN,ISH)=Td;

    if(ud[0]<0){
      std::cout << "neg velocity " << ISHPOIN <<'\n';
    }
  }
}

//--------------------------------------------------------------------------//

void AssignShockPointState3D::assignDownstreamAndUpstreamState()
{
  cout << "     => AssignShockPointState3D::assignDownstreamAndUpstreamState()\n";

  logfile("Assign downstream and upstream state\n");
  logfile("Number of shocks: ",(*nShocks),"\n");

  // define varID variable used to point the parameter vector variables
  // that assign upstream and downstream state
  // varID is chosen according to the gas model
  unsigned varID;
  if(ChemicalInfo::getModel()=="PG") { varID = 0; }
  else if (ChemicalInfo::getModel()=="TCneq") { varID = (*nsp)+1; }
  else {
    cout << "        AssignShockPointState3D::error => gas model not implemented\n";
    exit(1);
  }

  double dum;

  for(unsigned ISH=0;ISH<(*nShocks);ISH++) {
    logfile("\n------------------------\n");
    logfile("Shock nb => ",ISH+1,"\n");
    logfile("number of shock points is => ",nShockPoints->at(ISH),"\n\n");

    for(unsigned ISHPOIN=0;ISHPOIN<nShockPoints->at(ISH);ISHPOIN++) {

      logfile("Shock point: ",ISHPOIN+1);
      logfile(" -> upstream = ",(*ZRoeShu)(varID,ISHPOIN,ISH));
      logfile(", downstream = ",(*ZRoeShd)(varID,ISHPOIN,ISH),"\n");
      if( (*ZRoeShu)(varID,ISHPOIN,ISH)>(*ZRoeShd)(varID,ISHPOIN,ISH)) {
        logfile("(!) downstream and upstream state are inverted\n");
        for(unsigned IDOF=0;IDOF<(*ndof);IDOF++) {
          dum = (*ZRoeShu)(IDOF,ISHPOIN,ISH);
          (*ZRoeShu)(IDOF,ISHPOIN,ISH) = (*ZRoeShd)(IDOF,ISHPOIN,ISH);
          (*ZRoeShd)(IDOF,ISHPOIN,ISH) = dum;
        }
      }

      logfile("Upstream:: ");
      for(unsigned IDOF=0;IDOF<(*ndof);IDOF++) {
        logfile((*ZRoeShu)(IDOF,ISHPOIN,ISH), " ");
        // make a backup of the zroesh
        (*ZRoeShuOld)(IDOF,ISHPOIN,ISH) = (*ZRoeShu)(IDOF,ISHPOIN,ISH);
      }

      logfile("\nDownstream:: ");
      for(unsigned IDOF=0;IDOF<(*ndof);IDOF++) {
        // make a backup of the zroesh
        logfile((*ZRoeShd)(IDOF,ISHPOIN,ISH), " ");
        (*ZRoeShdOld)(IDOF,ISHPOIN,ISH) = (*ZRoeShd)(IDOF,ISHPOIN,ISH);
      }
      logfile("\n");
    }
  }
}

//--------------------------------------------------------------------------//

void AssignShockPointState3D::setSize()
{
  XYZShu.resize(PhysicsInfo::getnbDim(),
  PhysicsInfo::getnbShPointsMax(),
  (*nShocks));
  XYZShd.resize(PhysicsInfo::getnbDim(),
  PhysicsInfo::getnbShPointsMax(),
  (*nShocks));
  ZRoeShuOld->resize(PhysicsInfo::getnbDofMax(),
  PhysicsInfo::getnbShPointsMax(),
  PhysicsInfo::getnbShMax());
  ZRoeShdOld->resize((*ndof),
  PhysicsInfo::getnbShPointsMax(),
  PhysicsInfo::getnbShMax());
  normals.resize(PhysicsInfo::getnbDim(),
  PhysicsInfo::getnbShPointsMax(),
  PhysicsInfo::getnbShMax());
}

//--------------------------------------------------------------------------//

void AssignShockPointState3D::setAddress()
{
  unsigned totsize = npoin->at(0) + 2 * PhysicsInfo::getnbShMax() *
  PhysicsInfo::getnbShPointsMax();
  XYZ = new Array2D<double> (PhysicsInfo::getnbDim(),
  totsize,
  &coorVect->at(0));
  zroe = new Array2D<double> (PhysicsInfo::getnbDofMax(),
  totsize,
  &zroeVect->at(0));
  celnod = new Array2D<int> ((*nvt),nelem->at(0),
  &celnodVect->at(0));

  unsigned start = npoin->at(0)*PhysicsInfo::getnbDofMax();
  ZRoeShu =
  new Array3D <double> (PhysicsInfo::getnbDofMax(),
  PhysicsInfo::getnbShPointsMax(),
  PhysicsInfo::getnbShMax(),
  &zroeVect->at(start));
  start = npoin->at(0) * PhysicsInfo::getnbDofMax() +
  PhysicsInfo::getnbShPointsMax() *
  PhysicsInfo::getnbShMax() *
  PhysicsInfo::getnbDofMax();
  ZRoeShd =
  new Array3D <double> (PhysicsInfo::getnbDofMax(),
  PhysicsInfo::getnbShPointsMax(),
  PhysicsInfo::getnbShMax(),
  &zroeVect->at(start));
}

//--------------------------------------------------------------------------//

void AssignShockPointState3D::freeArray()
{
  delete XYZ; delete zroe; delete celnod;
  delete ZRoeShu; delete ZRoeShd;
}

//--------------------------------------------------------------------------//

void AssignShockPointState3D::setMeshData()
{
  nvt = MeshData::getInstance().getData<unsigned>("NVT");
  nelem = MeshData::getInstance().getData <vector<unsigned> > ("NELEM");
  npoin = MeshData::getInstance().getData <vector<unsigned> > ("NPOIN");
  zroeVect = MeshData::getInstance().getData <vector<double> >("ZROE");
  coorVect = MeshData::getInstance().getData <vector<double> >("COOR");
  celnodVect = MeshData::getInstance().getData <vector<int> >("CELNOD");
}

//--------------------------------------------------------------------------//

void AssignShockPointState3D::setPhysicsData()
{
  nsp = PhysicsData::getInstance().getData <unsigned> ("NSP");
  ndof = PhysicsData::getInstance().getData <unsigned> ("NDOF");
  nShocks = PhysicsData::getInstance().getData <unsigned> ("nShocks");
  nShockPoints =
  PhysicsData::getInstance().getData <vector <unsigned> > ("nShockPoints");
  nShockEdges =
  PhysicsData::getInstance().getData <vector <unsigned> > ("nShockEdges");    // decidere che fare con questo
  XYZSh = PhysicsData::getInstance().getData <Array3D <double> > ("XYZSH");
  ZRoeShuOld =
  PhysicsData::getInstance().getData <Array3D <double> > ("ZROESHuOLD");
  ZRoeShdOld =
  PhysicsData::getInstance().getData <Array3D <double> > ("ZROESHdOLD");

}

//--------------------------------------------------------------------------//

} // namespace ShockFitting
