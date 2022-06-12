///////////////////////////////////////////////////////
//Macro to select data
//
// This macro is in charge of passing all the filters needed to data files
//
/////////////////////////////////////////////////////
#define data_cxx
#include "data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <sstream>

ostringstream histo;

using namespace std;



void data::Loop()
{
 
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   //-------------------------------SET UP----------------------------------------//

   int icut[50] = {};

   //Set up of histograms

   TH1* hh[1000];

   //Histograms of data

   int ik=1;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Number of leptons (Nlep)",10,-0.5,9.5);

   ik=2;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Type of event (eee, eem, mme, mmm)",15,29.5,41.5);

   ik=3;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Pseudorrapidity (eta) for leptons",15,-3.,3.);

   ik=4;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Transverse momentum of leptons (lepton pt, GeV)",10,20.,200.);

   double pi = 4.*atan(1.);
   ik=5;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Azimuthal angle (phi) for leptons",16,-pi,pi);

   ik=6;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Energy of the lepton",15,0.,300.);

   ik=7;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Missing transverse momentum (Et,miss)",20,0.,160.);

   ik=8;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Transverse mass of W (Mt,W)",20,0.,200.);

   ik=9;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Mass of Z0 (GeV)",20,76.,106.);

   ik=10;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Tranverse momentum of the pair (ptll, GeV)",20,0.,400.);

   ik=11;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Charges of boson W",20,-2.,2.);

   ik=12;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Number of jets",10,0.,10.);

   ik=13;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Difference of phase between boson Z candidates",16,0.,pi);

   ik=14;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Difference of phase between Z candidate and W candidate (same charge)",16,0.,pi);

   ik=15;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Difference of phase between  Z candidate and W candidate (different charges)",16,0.,pi);

   ik=16;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Sum of difference of phase",16,0.,2*pi);

   ik=17;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Longitudinal momentum of the neutrino (GeV)",10,0.,600.);

   ik=18;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Total momentum of the neutrino (GeV)",10,0.,600.);
   
   ik=19;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Total transverse momentum of the system (GeV)",10,0.,400.);

   ik=20;
   histo.str(""); histo <<"hh";histo<<ik;
   cout<<"histo ="<<histo.str().c_str()<<endl;
   hh[ik]=new TH1D(histo.str().c_str(),"Total transverse momentum of the system with jets (GeV)",10,0.,400.);
   
   //---------------------------------------------------LOOP------------------------------------------------------//
   
   for (Long64_t jentry=0; jentry<nentries;jentry++)
     {
      
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      icut[0] = icut[0] + 1;//event counter
  
      //1. Classification of events into 4 categories: eee (33), eemu (35), mumue (37), mumumu (39)
      
      int event_type=0; 
      
      for(int i=0; i<lep_n; i++)
	{
	  event_type = event_type + lep_type->at(i);
	}

       //2. Application of needed triggers

       //electron trigger
      int itrige = 0;
      if(trigE) {itrige = 1;}

      //muon trigger
      int itrigm = 0;
      if(trigM) {itrigm = 1;}

      if(itrige== 0 && itrigm == 0)continue;
       
      if(itrige == 1){icut[1] = icut[1] + 1;}//electron trigger counter
    
      if(itrigm == 1){icut[2] = icut[2] + 1;}//muon trigger counter

      icut[3] = icut[3] + 1; //general trigger counter
      
      //3. Features of the 3 leptons in an array

      int ie = 0;
       
      double vce[20][10] = {}; //array with the information of the leptons
       
      for(int i=0; i<lep_n; i++) //loop over leptons
	{
	  
	  //4. pt>20GeV
	  if(lep_pt->at(i)/1000.<20.)continue;

	  //5. Cut in pseudorrapidity
	  if(lep_type->at(i)==11)
	    {
	      double feta = fabs(lep_eta->at(i));
	      if(feta > 2.47)continue;// feta<2.47 for electrons
	      if(feta > 1.37 && feta < 1.52)continue;// feta<1.37 and feta>1.52
	    }
	  
	  else if(lep_type->at(i)==13)
	    {
	      double feta = fabs(lep_eta->at(i));
	      if(feta > 2.5)continue;// feta<2.5 for muons
	    }

	   //7. Cut in quality
	  int itight=0;
	  if(lep_isTightID->at(i)){itight = 1;}
	  if(itight == 0)continue;

	  ie = ie + 1;//number of leptons involved

	  
	  //Features of involved leptons
	
	  vce[7][ie] = lep_pt->at(i)/1000; //transverse moment (pt)
	  vce[8][ie] = lep_eta->at(i); //pseudorrapidity
	  vce[6][ie] = lep_phi->at(i); // azimuthal angle
	  vce[4][ie] = lep_E->at(i)/1000; //energy
	  vce[9][ie] = lep_ptcone30->at(i)/1000; //trace isolation (ptiso)
	  vce[10][ie] = lep_etcone20->at(i)/1000; // calorimeter isolation (etiso)
	  //
	  vce[11][ie] = double(lep_charge->at(i)); //charge
	  vce[12][ie] = lep_isTightID->at(i); // tight quality
	  vce[14][ie] = lep_type->at(i);//type of lepton
	  //
	  vce[5][ie] = 2.*atan(exp(-vce[8][ie]));//theta
	  vce[13][ie] = vce[7][ie]/sin(vce[5][ie]);//p
	  //
	  vce[1][ie] = vce[4][ie]*sin(vce[5][ie])*cos(vce[6][ie]);//px
	  vce[2][ie] = vce[4][ie]*sin(vce[5][ie])*sin(vce[6][ie]);//py
	  vce[3][ie] = vce[4][ie]*cos(vce[5][ie]);//pz
 
	}
      
      if(ie != 3)continue;// cut in 3 leptons

       icut[7] = icut[7] + 1;
      
      //8. At least one of them has pt>25GeV
      if(vce[7][1]<25. && vce[7][2]<25. && vce[7][3]<25.)continue;

      icut[8] = icut[8] + 1; //counter for this condition
 
      //9. Mass of Z0 and transverse moment of the lepton pair
      
      // The pair of leptons born from the Z boson
      int c1; int c2;  //candidates for the Z
      int cW; //candidate for the W
      double z0[20][5] = {};// cuadrimoment of the z0
      double mll;//mass
      double ptll; //transverse moment of the pair

      // In order to find the pair of leptons, we take the possible pair combinations (12,13,23) and
      // take the pair of the Z as that whose invariant mass is the nearest to the theoretical mass of the Z

      TLorentzVector TZ0; 

      for(int i=1; i<=3; i++)
	{
	  //three possible combinations
	  if(i==1){c1=1; c2=2; cW=3;}
	  else if(i==2){c1=2; c2=3; cW=1;}
	  else if(i==3) {c1=3; c2=1; cW=2;}
	  
	  if(vce[11][c1]*vce[11][c2]>0)continue; //charges of the candidates must be opposite
	  if(vce[14][c1] != vce[14][c2])continue;//flavour of the candidates must be the same

	  //cuadrimoment
	  for(int j=1; j<=4; j++)
	    {
	      z0[j][i] = vce[j][c1]+vce[j][c2];
	    }

	  TZ0.SetPxPyPzE(z0[1][i],z0[2][i],z0[3][i],z0[4][i]);
	  z0[10][i] = TZ0.M();
	  z0[7][i] = TZ0.Pt();
	}

      int iz=0; double ZZ[20]={}; double Mmin=99999.;
      for(int i=1; i<=3; i++)
	{
	  double mZ = fabs(z0[10][i]-91.18);
	  if(mZ<Mmin)
	    {
	      Mmin=mZ;
	      iz=i;
	    }
	}

     if(iz==1){c1=1; c2=2; cW=3;}
     else if(iz==2){c1=2; c2=3; cW=1;}
     else if(iz==3) {c1=3; c2=1; cW=2;}
      
      diff_est = fabs(z0[10][iz]-91.18);
	
      if(diff_est>10.)continue; //those masses that differ more than 10 Gev from the Z0 mass are discarted

      mll = z0[10][iz];
      ptll = z0[7][iz];
      
      icut[9] = icut[9] + 1;

       //10. Miss transverse moment (energy and momentum of the neutrino)
      double miss_pt;
      miss_pt = met_et/1000.;
      if(miss_pt<30.)continue;//only a transverse missing energy greater that 30 GeV
      icut[11] = icut[11] + 1;

      //11. Transverse mass of the W boson
     
      double phi_nu;//phi of the neutrino
      double diff_phi;//phase difference between the neutrino and the lepton cW
      double mtW; //transverse mass of the W boson
      
      phi_nu = met_phi;
      diff_phi = fabs(phi_nu - vce[6][cW]);
      double teg = 180./pi;
      diff_phi = diff_phi*teg;
      
      if(diff_phi>180.){diff_phi = 360. - diff_phi;}

      diff_phi = diff_phi/teg;
      
      mtW = sqrt(2.*vce[7][cW]*miss_pt*(1.-cos(diff_phi)));
      //mtW=sqrt(2*pt_lepton*pt_neutrino*(1-cos(dif_phi)))
      
      if(mtW<30.)continue; //only mtW>30.GeV  

      icut[12] = icut[12] + 1;

      //Good events in histograms
      hh[1]->Fill(double(lep_n));//number of leptons (should be 3)
      hh[2]->Fill(double(event_type));
      
      for(int i=1; i<=ie; i++)
	{
	   //Histograms
	  hh[3]->Fill(vce[8][i]);//eta
	  hh[4]->Fill(vce[7][i]);//pt
	  hh[5]->Fill(vce[6][i]);//phi
	  hh[6]->Fill(vce[4][i]);//energy
	}
      
      hh[7]->Fill(miss_pt);
      hh[8]->Fill(mtW);
      hh[9]->Fill(mll);
      hh[10]->Fill(ptll);

      //Extra
      
      //Separate by charges
      hh[11]->Fill(vce[11][cW]);

      //Number of jets
      hh[12]->Fill(double(jet_n));

      //Angles between leptons
      double diff_tot=0.;

      // c1 & c2
      double phi_c1; double phi_c2;
      double diff_phi_Z;//diferencia de fases 
      
      phi_c1 = vce[6][c1];
      phi_c2 = vce[6][c2];
      diff_phi_Z = fabs(phi_c1 - phi_c2);
      diff_phi_Z = diff_phi_Z*teg;
      
      if(diff_phi_Z>180.){diff_phi_Z = 360. - diff_phi_Z;}

      diff_phi_Z = diff_phi_Z/teg;
      hh[13]->Fill(diff_phi_Z);

      diff_tot = diff_tot + diff_phi_Z;

      // c1 & cW 

      double phi_cW;

      phi_cW = vce[6][cW];
      diff_phi_Z = fabs(phi_c1 - phi_cW);
      diff_phi_Z = diff_phi_Z*teg;
      
      if(diff_phi_Z>180.){diff_phi_Z = 360. - diff_phi_Z;}

      diff_phi_Z = diff_phi_Z/teg; 

      if(vce[11][cW] == vce[11][c1])
	{
	  hh[14]->Fill(diff_phi_Z);
	}
      else
	{
	  hh[15]->Fill(diff_phi_Z);
	}

      diff_tot = diff_tot + diff_phi_Z;
      
      // c2 & cW
 
      diff_phi_Z = fabs(phi_c2 - phi_cW);
      diff_phi_Z = diff_phi_Z*teg;
      
      if(diff_phi_Z>180.){diff_phi_Z = 360. - diff_phi_Z;}

      diff_phi_Z = diff_phi_Z/teg;
      
      if(vce[11][cW] == vce[11][c2])
	{
	  hh[14]->Fill(diff_phi_Z);
	}
      else
	{
	  hh[15]->Fill(diff_phi_Z);
	}

      diff_tot = diff_tot + diff_phi_Z;
      hh[16]->Fill(diff_tot);
      
      //Longitudinal moment of the neutrino
      double p_nu_x=0.; double p_nu_y=0.; double p_nu_z=0.;
      for(int i=1; i<=ie; i++)
	{
	  p_nu_x = p_nu_x - vce[1][i];
	  p_nu_y = p_nu_y - vce[2][i];
	  p_nu_z = p_nu_z - vce[3][i];
	}
      
       double p_nu;
      //p_nu = sqrt(p_nu_x*p_nu_x+p_nu_y*p_nu_y+p_nu_z*p_nu_z);
      
      double pl_nu;
      //pl_nu = sqrt(p_nu*p_nu - miss_pt*miss_pt);
      pl_nu = fabs(p_nu_z);
      p_nu = sqrt(miss_pt*miss_pt+pl_nu*pl_nu);

      hh[17]->Fill(pl_nu);

      hh[18]->Fill(p_nu);

      //Total transverse moment of the system
      double pt_tot_x=0.; double pt_tot_y=0.;
      
      for(int i=1; i<=ie; i++)
	{
	  pt_tot_x = pt_tot_x + vce[7][i]*cos(vce[6][i]);
	  pt_tot_y = pt_tot_y + vce[7][i]*sin(vce[6][i]);
	}
      
      pt_tot_x = pt_tot_x + miss_pt*cos(phi_nu);
      pt_tot_y = pt_tot_y + miss_pt*sin(phi_nu);

      double pt_tot;
      pt_tot = sqrt(pt_tot_x*pt_tot_x+pt_tot_y*pt_tot_y);

      hh[19]->Fill(pt_tot);
      
   }

   //--WRAP UP--//
   cout<<"total number of events = "<<icut[0]<<endl;
   cout<<"total number of events with electron trigger = "<<icut[1]<<endl;
   cout<<"total number of events with muon trigger = "<<icut[2]<<endl;
   cout<<"total number of events with trigger = "<<icut[3]<<endl;
   cout<<"events with 3 leptons  = "<<icut[7]<<endl;
   cout<<"at least one with pt>25 GeV = "<<icut[8]<<endl;
   cout<<"quality cut= "<<icut[6]<<endl;
   cout<<"Z0 cut = "<<icut[9]<<endl;
   cout<<"ptmis cut = "<<icut[11]<<endl;
   cout<<"MT cut = "<<icut[12]<<endl;
   cout<<"Remaining events after all cuts =  "<<icut[10]<<endl;

   TFile myfile("data_histos.root","RECREATE");

    for(int j=1;j<=20;j++)
     {
       hh[j]->Write();
     }
}
