#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <iostream>
#include "TF2.h"
#include "TFile.h"
#include "TTree.h"

double me=511e-6;
double Mp=0.938272;
double mpi=0.1345;

double LeptonEnergy;
double HadronEnergy;
double theta_max;
double W2min;

bool collider;

double q2min=1;

TLorentzVector BeamL;
TLorentzVector BeamH;
TLorentzVector Pi0;
TLorentzVector Proton;
TLorentzVector ScattL;

double s;
double q2max;

TF2 *fangle ;TF2 *frapidity ;TF2 *fmomentum ;
TF2 *fangle_pr ;TF2 *frapidity_pr ;TF2 *fmomentum_pr ;

double y, theta_e,Ep, t_min, theta_cm;
double eta_el, p_el, angle_el;
//tmin 
double p_pi0_tmin,eta_pi0_tmin,angle_pi0_tmin;
double p_pr_tmin, pT_pr_tmin,eta_pr_tmin, angle_pr_tmin;
double p_el_tmin,eta_el_tmin;
double mx2_tmin, mE_tmin;
int T_station;
//t=-0.05 GeV^2
double p_pi0_50,eta_pi0_50_0,eta_pi0_50_180,angle_pi0_50;
double p_pr_50,eta_pr_50_0,eta_pr_50_180,pT_pr_50_0,pT_pr_50_180, angle_pr_50_0, angle_pr_50_180 ;
double mx2_50, mE_50;
//t=-0.1 GeV^2
double p_pi0_100,eta_pi0_100_0,eta_pi0_100_180,angle_pi0_100;
double p_pr_100,eta_pr_100_0,eta_pr_100_180,pT_pr_100_0,pT_pr_100_180, angle_pr_100_0, angle_pr_100_180 ;
double mx2_100, mE_100;
//t=-0.2 GeV^2
double p_pi0_200,eta_pi0_200_0,eta_pi0_200_180,angle_pi0_200;
double p_pr_200,eta_pr_200_0,eta_pr_200_180,pT_pr_200_0,pT_pr_200_180, angle_pr_200_0, angle_pr_200_180 ;
double mx2_200, mE_200;
//t=-0.5 GeV^2
double p_pi0_500,eta_pi0_500_0,eta_pi0_500_180,angle_pi0_500;
double p_pr_500,eta_pr_500_0,eta_pr_500_180,pT_pr_500_0,pT_pr_500_180, angle_pr_500_0, angle_pr_500_180 ;
double mx2_500, mE_500;
//t=-1 GeV^2
double p_pi0_1000,eta_pi0_1000_0,eta_pi0_1000_180,angle_pi0_1000;
double p_pr_1000,eta_pr_1000_0,eta_pr_1000_180,pT_pr_1000_0,pT_pr_1000_180, angle_pr_1000_0, angle_pr_1000_180;
double mx2_1000, mE_1000;
double phi_acc_1000;
double phi_acc_500;
double phi_acc_200;
double phi_acc_100;
double phi_acc_50;
//t=tmax
double p_pi0_tmax,eta_pi0_tmax,angle_pi0_tmax;
double p_pr_tmax,eta_pr_tmax, pT_pr_tmax;

double q2, xb;
double q2low, q2high, xblow, xbhigh;
double eps;

//Method to compute the lower and upper bound of the squared momentum transfer to the recoil baryon (Proton, Nucleus, Delta...). It depends on q2, xb, the mass of the recoil baryon mrecoil and the mass of the produced particle mpp 
Double_t tmin(Double_t q, Double_t x, double mrecoil, double mpp){
 Double_t W=TMath::Sqrt(q*(1/x-1)+Mp*Mp);
 Double_t E1cm=(W*W+q+Mp*Mp)/2./W;
 Double_t P1cm=TMath::Sqrt(E1cm*E1cm-Mp*Mp);
 Double_t E3cm=(W*W-mpp*mpp+mrecoil*mrecoil)/2./W;
 Double_t P3cm=TMath::Sqrt(E3cm*E3cm-mrecoil*mrecoil);
 Double_t result=(TMath::Power(q+mpp*mpp+Mp*Mp-mrecoil*mrecoil,2)/4./W/W-TMath::Power(P1cm-P3cm,2));
 return result;
}

Double_t tmax(Double_t q, Double_t x, double mrecoil, double mpp){
 Double_t W=TMath::Sqrt(q*(1/x-1)+Mp*Mp);
 Double_t E1cm=(W*W+q+Mp*Mp)/2./W;
 Double_t P1cm=TMath::Sqrt(E1cm*E1cm-Mp*Mp);
 Double_t E3cm=(W*W-mpp*mpp+mrecoil*mrecoil)/2./W;
 Double_t P3cm=TMath::Sqrt(E3cm*E3cm-mrecoil*mrecoil);
 Double_t result=(TMath::Power(q+mpp*mpp+Mp*Mp-mrecoil*mrecoil,2)/4./W/W-TMath::Power(P1cm+P3cm,2));
 return result;
}

//Setting LorentzVectors for hardon and lepton beam depending on the mode (collider=true, false means fixed-target). Also some angular/momentum bounds are set
void LoadConfig(bool mode, double Ene_L, double Ene_H){
  collider=mode;
  LeptonEnergy=Ene_L;
  HadronEnergy=Ene_H;
  if (collider){
    theta_max=177.9*TMath::DegToRad();
    BeamL.SetPxPyPzE(0,0,-LeptonEnergy,LeptonEnergy);
    BeamH.SetPxPyPzE(0,0,HadronEnergy,TMath::Sqrt(Mp*Mp+HadronEnergy*HadronEnergy));
    W2min=4;//(Mp+mpi)*(Mp+mpi);
    // if (HadronEnergy>90){
      q2low=0; 
      xblow=-5; xbhigh=0;
      /*}
    else{
       q2low=1;
       xblow=0.01; xbhigh=1;   
       }*/
  }
  else{
    theta_max=40*TMath::DegToRad();
    BeamL.SetPxPyPzE(0,0,LeptonEnergy,LeptonEnergy);
    BeamH.SetPxPyPzE(0,0,HadronEnergy,TMath::Sqrt(Mp*Mp+HadronEnergy*HadronEnergy));
    W2min=4;
    q2low=1; q2high=15;
    xblow=0.02; xbhigh=1;   
  }
 
  s=(BeamL+BeamH).Mag2();
}

//THE method which computes all LorentzVector. Q2 and xb has already been set and are kept in memory. Then remains t and phi to define the entire final set
//In truth... it must be validated on the recoil baryon
TLorentzVector Pi0Exclusive(double t, double phi_trento, double mrecoil, double mpp){

  theta_e=2*TMath::ATan(TMath::Sqrt(q2/(1-y)/4./LeptonEnergy/LeptonEnergy));
  Ep=q2/4./LeptonEnergy/TMath::Sin(theta_e/2.)/TMath::Sin(theta_e/2.);
  if (collider) ScattL.SetPxPyPzE(Ep*TMath::Sin(TMath::Pi()+theta_e),0,Ep*TMath::Cos(TMath::Pi()+theta_e),Ep);
  else ScattL.SetPxPyPzE(Ep*TMath::Sin(theta_e),0,Ep*TMath::Cos(theta_e),Ep);
  TLorentzVector q=BeamL-ScattL;
  TLorentzVector cms=q+BeamH;
  //q2=sxy is approximative... need to readjust the value of xb for consistency
  xb=-q.Mag2()/2./(BeamH*q);
  t_min=tmin(q2,xb,mrecoil,mpp);
  if (t==0) t=t_min;

  //The collision is easier to understand in the center-of-mass frame
  ScattL.Boost(-cms.BoostVector());
  q.Boost(-cms.BoostVector());

  //t_min-t=2q_cm q'_cm(1-cos (theta_gamma*/meson))
  double Ecm=(cms.Mag2()+mpp*mpp-mrecoil*mrecoil)/2./cms.Mag();
  double pcm=TMath::Sqrt(Ecm*Ecm-mpp*mpp);
  theta_cm=TMath::ACos(1-(t_min-TMath::Min(t_min,TMath::Max(t,tmax(q2,xb,mrecoil,mpp))))/2./pcm/q.Vect().Mag());
 
  double Ecm_pr=(cms.Mag2()+mrecoil*mrecoil-mpp*mpp)/2./cms.Mag();
  double pcm_pr=TMath::Sqrt(Ecm_pr*Ecm_pr-mrecoil*mrecoil);

  //Phi is a rotation of all four momentum around the axis defined by the virtual photon.
  TVector3 dir_pi0=q.Vect().Unit();
  TVector3 rot_dir;
  if (t==t_min) rot_dir.SetXYZ(dir_pi0.Z(),0,-dir_pi0.X());
  else rot_dir=(q.Vect().Cross(ScattL.Vect())).Unit();
  dir_pi0.Rotate(theta_cm,rot_dir);
  dir_pi0.Rotate(phi_trento,q.Vect().Unit());
 
  Pi0.SetPxPyPzE(pcm*dir_pi0.X(),pcm*dir_pi0.Y(), pcm*dir_pi0.Z(),Ecm);
  Proton.SetPxPyPzE(-pcm_pr*dir_pi0.X(),-pcm_pr*dir_pi0.Y(),-pcm_pr*dir_pi0.Z(),Ecm_pr);

  //Back to lab frame
  Proton.Boost(cms.BoostVector());
  Pi0.Boost(cms.BoostVector());
  q.Boost(cms.BoostVector());
  ScattL.Boost(cms.BoostVector());
  
  return Pi0;
}

//y the inelasticity which is more convenient to use than xb.
double ymax(double q2){
  return TMath::Min(1-q2/4./LeptonEnergy/LeptonEnergy/TMath::Sin(theta_max/2.)/TMath::Sin(theta_max/2.),0.95);
}

double ymin(double q2){
  if (collider) return TMath::Max((W2min-Mp*Mp+q2)/s,0.01);
  else return (W2min-Mp*Mp+q2)/s;
}

double epsilon(TLorentzVector q){
  return 1/(1-2*q.Vect().Mag2()/q.Mag2()*TMath::Tan(theta_e/2.)*TMath::Tan(theta_e/2.));
}

//Compute q2 max
void ComputeQ2Max(){
  double q2temp; double q2step;
  q2temp=1;
  if(collider){
    q2step=1;
  }
  else{
    q2step=0.1;
  }
  while (ymin(q2temp)<=ymax(q2temp)){
    q2max=q2temp;
    q2temp+=q2step;
  }
  if(collider) q2high=TMath::Log10(q2max);
  else q2high=q2max;
}

//Geometrical acceptance funciton for the proton
bool TestProton(){
  bool IsDetected=false;
  double ang=TMath::ACos(Proton.Pz()/Proton.P());
  if (ang>0.035||(ang>0.006&&ang<0.020)||(ang>0.001&&ang<0.005)) IsDetected=true;
  return IsDetected;
}

//This function is the main function... It will compute in a grid of q2/xb fr various t and phi the lorentzvector in lab frame.
void FillTree(bool mode, double EneL, double EneH){
  LoadConfig(mode,EneL,EneH);
  ComputeQ2Max();
 
  TFile *f;
  if (collider) f=new TFile(Form("collider_%d_%d.root",(int) EneL,(int) EneH),"recreate");
  else f=new TFile("clas12.root","recreate");

  TTree* T=new TTree("T","Exclusive pi0");
  T->Branch("q2",&q2,"q2/D");
  T->Branch("xb",&xb,"xb/D");
  T->Branch("y",&y,"y/D");
  T->Branch("t_min",&t_min,"t_min/D");
  //t=tmin
  T->Branch("p_pi0_tmin",&p_pi0_tmin,"p_pi0_tmin/D");
  T->Branch("eta_pi0_tmin",&eta_pi0_tmin,"eta_pi0_tmin/D");
  T->Branch("angle_pi0_tmin",&angle_pi0_tmin,"angle_pi0_tmin/D");
  T->Branch("p_pr_tmin",&p_pr_tmin,"p_pr_tmin/D");
  T->Branch("pT_pr_tmin",&pT_pr_tmin,"pT_pr_tmin/D");
  T->Branch("mx2_tmin",&mx2_tmin,"mx2_tmin/D");
  T->Branch("mE_tmin",&mE_tmin,"mE_tmin/D");
  T->Branch("eta_pr_tmin",&eta_pr_tmin,"eta_pr_tmin/D");
  T->Branch("T_station",&T_station,"T_station/I");
  T->Branch("angle_pr_tmin",&angle_pr_tmin,"angle_pr_tmin/D");
   //t=tmax
  T->Branch("p_pi0_tmax",&p_pi0_tmax,"p_pi0_tmax/D");
  T->Branch("eta_pi0_tmax",&eta_pi0_tmax,"eta_pi0_tmax/D");
  T->Branch("angle_pi0_tmax",&angle_pi0_tmax,"angle_pi0_tmax/D");
  T->Branch("p_pr_tmax",&p_pr_tmax,"p_pr_tmax/D");
   T->Branch("pT_pr_tmax",&pT_pr_tmax,"pT_pr_tmax/D");
  T->Branch("eta_pr_tmax",&eta_pr_tmax,"eta_pr_tmax/D");
   //t=-0.05
  T->Branch("p_pi0_50",&p_pi0_50,"p_pi0_50/D");
  T->Branch("eta_pi0_50_0",&eta_pi0_50_0,"eta_pi0_50_0/D");
  T->Branch("eta_pi0_50_180",&eta_pi0_50_180,"eta_pi0_50_180/D");
  T->Branch("angle_pi0_50",&angle_pi0_50,"angle_pi0_50/D");
  T->Branch("p_pr_50",&p_pr_50,"p_pr_50/D");
  T->Branch("eta_pr_50_0",&eta_pr_50_0,"eta_pr_50_0/D");
  T->Branch("eta_pr_50_180",&eta_pr_50_180,"eta_pr_50_180/D");
   T->Branch("angle_pr_50_0",&angle_pr_50_0,"angle_pr_50_0/D");
     T->Branch("angle_pr_50_180",&angle_pr_50_180,"angle_pr_50_180/D");
     T->Branch("mx2_50",&mx2_50,"mx2_50/D");
      T->Branch("phi_acc_50",&phi_acc_50,"phi_acc_50/D");
      T->Branch("phi_acc_500",&phi_acc_500,"phi_acc_500/D");
  T->Branch("mE_50",&mE_50,"mE_50/D");
   //t=-0.1
  T->Branch("p_pi0_100",&p_pi0_100,"p_pi0_100/D");
  T->Branch("eta_pi0_100_0",&eta_pi0_100_0,"eta_pi0_100_0/D");
  T->Branch("eta_pi0_100_180",&eta_pi0_100_180,"eta_pi0_100_180/D");
  T->Branch("angle_pi0_100",&angle_pi0_100,"angle_pi0_100/D");
  T->Branch("p_pr_100",&p_pr_100,"p_pr_100/D");
  T->Branch("eta_pr_100_0",&eta_pr_100_0,"eta_pr_100_0/D");
  T->Branch("eta_pr_100_180",&eta_pr_100_180,"eta_pr_100_180/D");
   T->Branch("angle_pr_100_0",&angle_pr_100_0,"angle_pr_100_0/D");
     T->Branch("angle_pr_100_180",&angle_pr_100_180,"angle_pr_100_180/D");
     T->Branch("mx2_100",&mx2_100,"mx2_100/D");
      T->Branch("phi_acc_100",&phi_acc_100,"phi_acc_100/D");
      T->Branch("phi_acc_1000",&phi_acc_1000,"phi_acc_1000/D");
  T->Branch("mE_100",&mE_100,"mE_100/D");
  //t=-0.2
  T->Branch("p_pi0_200",&p_pi0_200,"p_pi0_200/D");
  T->Branch("eta_pi0_200_0",&eta_pi0_200_0,"eta_pi0_200_0/D");
  T->Branch("eta_pi0_200_180",&eta_pi0_200_180,"eta_pi0_200_180/D");
  T->Branch("angle_pi0_200",&angle_pi0_200,"angle_pi0_200/D");
  T->Branch("p_pr_200",&p_pr_200,"p_pr_200/D");
  T->Branch("eta_pr_200_0",&eta_pr_200_0,"eta_pr_200_0/D");
  T->Branch("eta_pr_200_180",&eta_pr_200_180,"eta_pr_200_180/D");
   T->Branch("angle_pr_200_0",&angle_pr_200_0,"angle_pr_200_0/D");
     T->Branch("angle_pr_200_180",&angle_pr_200_180,"angle_pr_200_180/D");
     T->Branch("mx2_200",&mx2_200,"mx2_200/D");
      T->Branch("phi_acc_200",&phi_acc_200,"phi_acc_200/D");
      T->Branch("phi_acc_1000",&phi_acc_1000,"phi_acc_1000/D");
  T->Branch("mE_200",&mE_200,"mE_200/D");
  //t=-0.5*q2
  T->Branch("p_pi0_500",&p_pi0_500,"p_pi0_500/D");
  T->Branch("eta_pi0_500_0",&eta_pi0_500_0,"eta_pi0_500_0/D");
  T->Branch("eta_pi0_500_180",&eta_pi0_500_180,"eta_pi0_500_180/D");
  T->Branch("angle_pi0_500",&angle_pi0_500,"angle_pi0_500/D");
  T->Branch("p_pr_500",&p_pr_500,"p_pr_500/D");
  T->Branch("eta_pr_500_0",&eta_pr_500_0,"eta_pr_500_0/D");
  T->Branch("eta_pr_500_180",&eta_pr_500_180,"eta_pr_500_180/D");
   T->Branch("angle_pr_500_0",&angle_pr_500_0,"angle_pr_500_0/D");
     T->Branch("angle_pr_500_180",&angle_pr_500_180,"angle_pr_500_180/D");
     T->Branch("mx2_500",&mx2_500,"mx2_500/D");
      T->Branch("phi_acc_500",&phi_acc_500,"phi_acc_500/D");
      T->Branch("phi_acc_1000",&phi_acc_1000,"phi_acc_1000/D");
  T->Branch("mE_500",&mE_500,"mE_500/D");
   //t=-1*q2
  T->Branch("p_pi0_1000",&p_pi0_1000,"p_pi0_1000/D");
  T->Branch("eta_pi0_1000_0",&eta_pi0_1000_0,"eta_pi0_1000_0/D");
  T->Branch("eta_pi0_1000_180",&eta_pi0_1000_180,"eta_pi0_1000_180/D");
  T->Branch("angle_pi0_1000",&angle_pi0_1000,"angle_pi0_1000/D");
  T->Branch("p_pr_1000",&p_pr_1000,"p_pr_1000/D");
  T->Branch("eta_pr_1000_0",&eta_pr_1000_0,"eta_pr_1000_0/D");
  T->Branch("eta_pr_1000_180",&eta_pr_1000_180,"eta_pr_1000_180/D");
    T->Branch("angle_pr_1000_0",&angle_pr_1000_0,"angle_pr_1000_0/D");
     T->Branch("angle_pr_1000_180",&angle_pr_1000_180,"angle_pr_1000_180/D");
     T->Branch("mx2_1000",&mx2_1000,"mx2_1000/D");
  T->Branch("mE_1000",&mE_1000,"mE_1000/D");
  T->Branch("p_el",&p_el,"p_el/D");
  T->Branch("eta_el",&eta_el,"eta_el/D");
   T->Branch("angle_el",&angle_el,"angle_el/D");
  T->Branch("eps",&eps,"eps/D");

  for (int q=0;q<1000;q++){
    for (int x=0;x<1000;x++){
      if (collider){
	//Grid in log10 of q2 and xb
	q2=TMath::Power(10,q2low+q*(q2high-q2low)/1000.);
	xb=TMath::Power(10,xblow+x*(xbhigh-xblow)/1000.);
      }
      else{
	q2=q2low+q*(q2high-q2low)/1000.;
	xb=xblow+x*(xbhigh-xblow)/1000.;
      }
      	y=q2/s/xb;

      if (q2>=q2min&&q2<q2max&&y<ymax(q2)&&y>ymin(q2)){
	//t=tmin
	T_station=0;
	Pi0Exclusive(0,0);
	angle_el=theta_e*TMath::RadToDeg();
	p_pi0_tmin=Pi0.P(); eta_pi0_tmin=Pi0.Rapidity(); angle_pi0_tmin=2*TMath::ASin(mpi/Pi0.E());
	p_pr_tmin=Proton.P(); eta_pr_tmin=Proton.Rapidity(); pT_pr_tmin=Proton.Pt(); angle_pr_tmin=TMath::ACos(Proton.Pz()/Proton.P());
	mx2_tmin=(BeamL+BeamH-Pi0-Proton-ScattL).Mag2(); mE_tmin=(BeamL+BeamH-Pi0-Proton-ScattL).E();
	if (angle_pr_tmin>0.001&&angle_pr_tmin<0.005) T_station=1;
	if (angle_pr_tmin>0.006&&angle_pr_tmin<0.02) T_station=2;
	if (angle_pr_tmin>0.035) T_station=3;

	//t=-0.05 GeV^2
	Pi0Exclusive(t_min-0.05,0);
	p_pi0_50=Pi0.P(); eta_pi0_50_0=Pi0.Rapidity(); angle_pi0_50=2*TMath::ASin(mpi/Pi0.E());
	p_pr_50=Proton.P(); eta_pr_50_0=Proton.Rapidity();angle_pr_50_0=TMath::ACos(Proton.Pz()/Proton.P());
	Pi0Exclusive(t_min-0.05,TMath::Pi());
	eta_pi0_50_180=Pi0.Rapidity();
	eta_pr_50_180=Proton.Rapidity();
	angle_pr_50_180=TMath::ACos(Proton.Pz()/Proton.P());
	mx2_50=(BeamL+BeamH-Pi0-Proton-ScattL).Mag2(); mE_50=(BeamL+BeamH-Pi0-Proton-ScattL).E();
	phi_acc_50=0;
	for (int i=0; i<100; i++){
	  Pi0Exclusive(t_min-0.05,i*2*TMath::Pi()/100.);
	  if (TestProton()) phi_acc_50++;
	}
	phi_acc_50=phi_acc_50/100.;

	//t=-0.1 GeV^2
	Pi0Exclusive(t_min-0.1,0);
	p_pi0_100=Pi0.P(); eta_pi0_100_0=Pi0.Rapidity(); angle_pi0_100=2*TMath::ASin(mpi/Pi0.E());
	p_pr_100=Proton.P(); eta_pr_100_0=Proton.Rapidity();angle_pr_100_0=TMath::ACos(Proton.Pz()/Proton.P());
	Pi0Exclusive(t_min-0.1,TMath::Pi());
	eta_pi0_100_180=Pi0.Rapidity();
	eta_pr_100_180=Proton.Rapidity();
	angle_pr_100_180=TMath::ACos(Proton.Pz()/Proton.P());
	mx2_100=(BeamL+BeamH-Pi0-Proton-ScattL).Mag2(); mE_100=(BeamL+BeamH-Pi0-Proton-ScattL).E();
	phi_acc_100=0;
	//quick computation of the phi trento (angle between leptonic and hadronic plane)  acceptance for proton detection
	for (int i=0; i<100; i++){
	  Pi0Exclusive(t_min-0.1,i*2*TMath::Pi()/100.);
	  if (TestProton()) phi_acc_100++;
	}
	phi_acc_100=phi_acc_100/100.;
	
	//t=-0.2 GeV^2
	Pi0Exclusive(t_min-0.2,0);
	p_pi0_200=Pi0.P(); eta_pi0_200_0=Pi0.Rapidity(); angle_pi0_200=2*TMath::ASin(mpi/Pi0.E());
	p_pr_200=Proton.P(); eta_pr_200_0=Proton.Rapidity();angle_pr_200_0=TMath::ACos(Proton.Pz()/Proton.P());
	Pi0Exclusive(t_min-0.2,TMath::Pi());
	eta_pi0_200_180=Pi0.Rapidity();
	eta_pr_200_180=Proton.Rapidity();
	angle_pr_200_180=TMath::ACos(Proton.Pz()/Proton.P());
	mx2_200=(BeamL+BeamH-Pi0-Proton-ScattL).Mag2(); mE_200=(BeamL+BeamH-Pi0-Proton-ScattL).E();
	phi_acc_200=0;
	for (int i=0; i<100; i++){
	  Pi0Exclusive(t_min-0.2,i*2*TMath::Pi()/100.);
	  if (TestProton()) phi_acc_200++;
	}
	phi_acc_200=phi_acc_200/100.;

	//t=-0.5 GeV^2
	Pi0Exclusive(t_min-0.5,0);
	p_pi0_500=Pi0.P(); eta_pi0_500_0=Pi0.Rapidity(); angle_pi0_500=2*TMath::ASin(mpi/Pi0.E());
	p_pr_500=Proton.P(); eta_pr_500_0=Proton.Rapidity();angle_pr_500_0=TMath::ACos(Proton.Pz()/Proton.P());
	Pi0Exclusive(t_min-0.5,TMath::Pi());
	eta_pi0_500_180=Pi0.Rapidity();
	eta_pr_500_180=Proton.Rapidity();
	angle_pr_500_180=TMath::ACos(Proton.Pz()/Proton.P());
	mx2_500=(BeamL+BeamH-Pi0-Proton-ScattL).Mag2(); mE_500=(BeamL+BeamH-Pi0-Proton-ScattL).E();
	phi_acc_500=0;
	for (int i=0; i<100; i++){
	  Pi0Exclusive(t_min-0.5,i*2*TMath::Pi()/100.);
	  if (TestProton()) phi_acc_500++;
	}
	phi_acc_500=phi_acc_500/100.;
	//t=-1 GeV^2
	Pi0Exclusive(t_min-1,0);
	p_pi0_1000=Pi0.P(); eta_pi0_1000_0=Pi0.Rapidity(); angle_pi0_1000=2*TMath::ASin(mpi/Pi0.E());
	p_pr_1000=Proton.P(); eta_pr_1000_0=Proton.Rapidity();angle_pr_1000_0=TMath::ACos(Proton.Pz()/Proton.P());

	Pi0Exclusive(t_min-1,TMath::Pi());
	eta_pi0_1000_180=Pi0.Rapidity();
	eta_pr_1000_180=Proton.Rapidity();
	angle_pr_1000_180=TMath::ACos(Proton.Pz()/Proton.P());
	phi_acc_1000=0;
	for (int i=0; i<100; i++){
	  Pi0Exclusive(t_min-1,i*2*TMath::Pi()/100.);
	  if (TestProton()) phi_acc_1000++;
	}
	phi_acc_1000=phi_acc_1000/100.;
	mx2_1000=(BeamL+BeamH-Pi0-Proton-ScattL).Mag2(); mE_1000=(BeamL+BeamH-Pi0-Proton-ScattL).E();
	p_el=ScattL.P(); eta_el=ScattL.Rapidity();
	eps=epsilon(BeamL-ScattL);

	//t=tmax GeV^2
	Pi0Exclusive(tmax(q2, xb, Mp, mpi),0);
	p_pi0_tmax=Pi0.P(); eta_pi0_tmax=Pi0.Rapidity(); angle_pi0_tmax=2*TMath::ASin(mpi/Pi0.E());
	p_pr_tmax=Proton.P(); eta_pr_tmax=Proton.Rapidity(); pT_pr_tmax=Proton.Pt();

	//If everything is alright, save the point (because q2=sxy is approximate, we may have a few kinematic configurations Ok at generation but wrong when xb is accurately computed.)
	if (y<ymax(q2)&&y>ymin(q2)) T->Fill();
      }
    }
  }
  f->cd();
  T->Write("T");
  f->Close();
  return;
}
