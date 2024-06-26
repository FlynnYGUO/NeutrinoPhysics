#include <iostream>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <vector>
using namespace std;
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"

vector<vector<vector<double>>>* xyz_pos=nullptr;
vector<vector<vector<double>>>* xyz_mom=nullptr;
double geoeff_cut = 0.01;
double numu_e=0;
double e_vis_true=0;
const int NUM_VTX=22, NUM_LAR_DTR=15;
double LAr_position[NUM_LAR_DTR]={-2800.,-2575.,-2400.,-2175.,-2000.,-1775.,-1600.,-1375.,-1200.,-975.,-800.,-575.,-400.,-175.,0.};
double vertex_position[NUM_VTX]={-299.,-292.,-285.,-278.,-271.,-264.,-216.,-168.,-120.,-72.,-24.,24.,72.,120.,168.,216.,264.,271.,278.,285.,292.,299.};
double total_detected[5][NUM_LAR_DTR][NUM_VTX]={};
//float scale[30]={.19,.18,.18,.7,.7,1.05,.04,.034,.034,.07,.07,.07,.23,.21,.21,.73,.65,1.05,1.,1.,1.,1.05,1.05,1.05,.23,.2,.2,.7,.7,1.};
float scale[35]={.14,.12,.12,.16,.14,.14,.14,.36,.24,.24,.38,.27,.27,.25,.5,.35,.35,.54,.4,.38,.38,.6,.56,.56,1.05,1.05,1.05,1.05,.3,.19,.19,.54,.37,.37,.37};

struct Para
{
  //static constexpr const char *const S;
  //constexpr const *char , VTX_X="vtx_x", *VTX_Y="vtx_y", *VTX_Z="vtx_z";
  //const char *LMX="LepMomX", *LMY="LepMomY", *LMZ="LepMomZ";
  char field[20];
  double l;
  double h;
  //vector<vector<vector<double>>>* field_value;
};

Para pr[]= //position is in units of cm, momentum is in units of GeV/c, angle is in units of rad, and energy is in  units of GeV //check if the TTrees have LepNuAngle
{ //match the x-ranges with the PRISM histograms' x-ranges
  {"vtx_x", -300., 300.},
  {"vtx_y", 5.385, 5.39},
  {"vtx_z",  659.995, 660.005},
  {"LepMomX", -0.6, 0.4},
  {"LepMomY", -0.8, 0.4},
  {"LepMomZ", 0., 7.},
  {"LepMomTot", 0., 7.}
  //{"cos_LepNuAngle", 0., 1.}
  //{"LongMom", 0., 7.}
  //{"ND_Gen_numu_E", 0., 10.},
  //{"ND_E_vis_true", 0., 10.}
};

struct sel_type
{
  const char* sel_name;
  const char* eff_name;
  vector<vector<double>>* eff_value=nullptr;
  sel_type() {}
  sel_type(const char* sn, const char* en)
  :sel_name(sn),eff_name(en) {}
};

vector<sel_type> br=
{
  sel_type("muon_contained", "muon_contained_eff"),
  sel_type("muon_tracker", "muon_tracker_eff"),
  sel_type("muon_selected", "muon_selected_eff"),
  sel_type("hadron_selected", "hadron_selected_eff"),
  sel_type("combined", "combined_eff")
};

void populate_histograms(char* eff,char* caf,vector<vector<TH1D*>>& hists1,vector<vector<TH1D*>>& hists2)
{
  TFile eff_file(eff);
  TFile caf_file(caf);
  TTree *event_data=(TTree*)eff_file.Get("event_data");
  TTree *thing=(TTree*)caf_file.Get("throwResults");
  gSystem->Exec("rm -f AutoDict*vector*vector*vector*double*");
  gInterpreter->GenerateDictionary("vector<vector<vector<double>>>", "vector");
  for(auto& sel:br) event_data->SetBranchAddress(sel.eff_name, &(sel.eff_value)); //efficiencies are in 3d array, but energy is in 1d array
  thing->SetBranchAddress("FD_evt_NDLAr_OffAxis_Sim_lep_start_v", &xyz_pos);
  thing->SetBranchAddress("FD_evt_NDLAr_OffAxis_Sim_lep_start_p", &xyz_mom);
  //thing->SetBranchAddress("ND_Gen_numu_E", &numu_e);
  //thing->SetBranchAddress("ND_E_vis_true", &e_vis_true);

  Long64_t nentries1=event_data->GetEntries();
  Long64_t nentries2=thing->GetEntries();
  if (nentries1!=nentries2) {cout<<"The efficiency file "<<eff<<" has "<<nentries2<<" events, and the CAF file "<<caf<<" has "<<nentries1<<" events."<<endl;}
  for (int i=0;i<nentries2;i++) {
    event_data->GetEntry(i);
    thing->GetEntry(i);
    unsigned long lar_pos=14; // on axis pos
    int k=0;
    for (Para item:pr) {
      double var_type=0.0;
      //if (k==7) {var_type=numu_e;}
      //if (k==8) {var_type=e_vis_true;}

      for (unsigned long vtx_pos=0;vtx_pos<NUM_VTX;vtx_pos++)
      {
	      int n=0;
        for (auto& sel:br)
        {
          TH1D* hist1=hists1[n][k];
          TH1D* hist2=hists2[n][k];
          n++;
          if (k<3) var_type=(*xyz_pos)[lar_pos][vtx_pos][k];
          else if (k<6) var_type=(*xyz_mom)[lar_pos][vtx_pos][k-3];
          else if (k==6) {var_type=sqrt(pow((*xyz_mom)[lar_pos][vtx_pos][0],2)+pow((*xyz_mom)[lar_pos][vtx_pos][1],2)+pow((*xyz_mom)[lar_pos][vtx_pos][2],2));cout<<var_type<<endl;}
          vector<vector<double>>* eff_value=sel.eff_value;
          vector<vector<double>>& eff_value2=*eff_value;
          double geo_eff=eff_value2[lar_pos][vtx_pos];
          if (geo_eff>1.) {cout<<"efficiency of event "<<i<<" at position "<<LAr_position[lar_pos]<<", "<<vertex_position[vtx_pos]<<" is "<<geo_eff<<endl;}
          hist1->Fill(var_type);
          if (geo_eff<=geoeff_cut) {
            continue;
          } else {
            hist2->Fill(var_type, geo_eff);
          }
        }
      }// end vtx_pos loop
      k++;
    } // end para item:pr
  }
  eff_file.Close();
  caf_file.Close();
}

void FD_selection_cuts_Bnewversion()
{
  char eff[9999];
  char caf[9999];
  vector<vector<TH1D*>> histograms1;
  vector<vector<TH1D*>> histograms2;
  for(auto& sel:br)
  {
    const char* dt=sel.sel_name;
    vector<TH1D*> histset1, histset2;
    histograms1.push_back(histset1);
    histograms2.push_back(histset2);
    int i=0;
    for (Para item:pr)
    {
      double lowerbound=item.l;
      double upperbound=item.h;
      histograms1.back().push_back(new TH1D(Form("%s_hist_%d",dt,i), Form("raw %s %d", dt, i), 200, lowerbound, upperbound));
      histograms2.back().push_back(new TH1D(Form("%s_hist_%d",dt,i), Form("selected %s %d", dt, i), 200, lowerbound, upperbound));
    i++;
    }
  }

  for (int i=1; i<11; i++)
  {
    for (int j=0; j<5; j++)
    {
      memset(eff, 0, 9999); //clear array each time
      memset(caf, 0, 9999);
      sprintf(eff, "/storage/shared/barwu/FDCAFIntegrationEffFiles/old/%d_p%01d_Eff.root", i, j);
      sprintf(caf, "/storage/shared/barwu/10thTry/FDCAFIntegration4GEC_splitfiles2/%d_p%01d.root", i, j);
      if(access(eff, 0)==0)
      {
        populate_histograms(eff,caf,histograms1,histograms2);
      } else {
        cout<<"Error reading file "<<eff<<endl;
        continue;
      }
    }
  }

  //gStyle->SetOptStat(000000000);
  gStyle->SetOptStat(111111111);
  TCanvas *c=new TCanvas("c","FD-in-ND all graphs",8000, 4000);
  TCanvas *r=new TCanvas("r","Ratio Plots",8000, 4000);
  c->Divide(7,5);
  r->Divide(7,5);
  int i=1;
  int i_select=0;
  for (auto& sel:br)
  {
    const char *dt=sel.sel_name;
    int n=0;
    for(Para item:pr)
    {
      const char *fd=item.field;
      double lowerbound=item.l;
      double upperbound=item.h;
      TVirtualPad *pad=c->cd(i);
      TH1D *hist2=histograms2[i_select][n];
      hist2->SetLineColor(kTeal-3);
      hist2->Draw("hist");
      TH1D *hist1=histograms1[i_select][n];
      hist1->SetLineColor(kPink);
      hist1->Draw("samehist");

      float max1=hist1->GetMaximum();
      hist2->SetAxisRange(lowerbound,upperbound,"X");
      hist2->SetAxisRange(0.,1.16*max1,"Y");
      hist2->SetTitle(Form("%s %s Selection Cut", fd, dt));
      hist2->GetXaxis()->SetTitle(fd);
      hist2->GetYaxis()->SetTitle("# of events");
      TLegend *leg=new TLegend(0.1,0.8,0.35,0.9);
      leg->SetHeader("comparison");
      leg->AddEntry(hist1, "raw distribution");
      leg->AddEntry(hist2, "selection-cut distribution");
      leg->Draw();

      r->cd(i);
      TH1D *rplot=(TH1D*)hist2->Clone();
      rplot->Divide(hist1);
      rplot->SetAxisRange(0.,scale[i-1],"Y");
      rplot->SetLineColor(kBlue);
      rplot->Draw("hist");
      n++;
      i++;
    }
    c->Update();
    r->Update();
    i_select++;
  }
  c->SaveAs(Form("/home/fyguo/testbaroncode/fdhist/%0.3f_eff_veto_cut_FDinND_old_all_hists.png",geoeff_cut));
  r->SaveAs(Form("/home/fyguo/testbaroncode/fdhist/%0.3f_eff_veto_cut_old_all_hists_ratios.png",geoeff_cut));

  TCanvas *cs[5];
  TCanvas *rs[5];
  i=1;
  for (auto& sel:br)
  {
    const char *dt=sel.sel_name;
    cs[i-1]=new TCanvas(Form("c%01d",i),dt,8000, 4000);
    cs[i-1]->Divide(4,2);
    rs[i-1]=new TCanvas(Form("r%01d",i),dt,8000, 4000);
    rs[i-1]->Divide(4,2);
    int n=0;
    for(Para& item:pr)
    {
      const char *fd=item.field;
      double lowerbound=item.l;
      double upperbound=item.h;
      TVirtualPad *pads=cs[i-1]->cd(n+1);
      TH1D *hist2=histograms2[i-1][n];
      hist2->SetLineColor(kTeal-3);
      hist2->Draw("hist");
      TH1D *hist1=histograms1[i-1][n];
      hist1->SetLineColor(kPink);
      hist1->Draw("samehist");

      float max1=hist1->GetMaximum();
      hist2->SetAxisRange(lowerbound,upperbound,"X");
      hist2->SetAxisRange(0.,1.16*max1,"Y");
      hist2->SetTitle(Form("%s %s Selection Cut", fd, dt));
      hist2->GetXaxis()->SetTitle(fd);
      hist2->GetYaxis()->SetTitle("# of events");
      TLegend *leg=new TLegend(0.1,0.8,0.35,0.9);
      leg->SetHeader("comparison");
      leg->AddEntry(hist1, "raw distribution");
      leg->AddEntry(hist2, "selection-cut distribution");
      leg->Draw();

      rs[i-1]->cd(n+1);
      TH1D *rplot=(TH1D*)hist2->Clone();
      rplot->Divide(hist1);
      rplot->SetAxisRange(0.,scale[(i-1)*7+n],"Y");
      rplot->SetLineColor(kBlue);
      rplot->Draw("hist");
      n++;
    }
    cs[i-1]->Update();
    rs[i-1]->Update();
    cs[i-1]->SaveAs(Form("/home/fyguo/testbaroncode/fdhist/%0.3f_eff_veto_cut_FDinND_old_%s_hists.png", geoeff_cut, dt));
    rs[i-1]->SaveAs(Form("/home/fyguo/testbaroncode/fdhist/%0.3f_eff_veto_cut_FDinND_old_%s_hists_ratios.png", geoeff_cut, dt));
    i++;
  }
}
