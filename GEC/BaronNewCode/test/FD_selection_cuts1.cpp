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
double lepnuangle=0.;
double numu_e=0;
double e_vis_true=0;
const int NUM_VTX=22, NUM_LAR_DTR=15;
double LAr_position[NUM_LAR_DTR]={-2800.,-2575.,-2400.,-2175.,-2000.,-1775.,-1600.,-1375.,-1200.,-975.,-800.,-575.,-400.,-175.,0.};
double vertex_position[NUM_VTX]={-299.,-292.,-285.,-278.,-271.,-264.,-216.,-168.,-120.,-72.,-24.,24.,72.,120.,168.,216.,264.,271.,278.,285.,292.,299.};
double total_detected[5][NUM_LAR_DTR][NUM_VTX]={};
//float scale[30]={.19,.18,.18,.7,.7,1.05,.04,.034,.034,.07,.07,.07,.23,.21,.21,.73,.65,1.05,1.,1.,1.,1.05,1.05,1.05,.23,.2,.2,.7,.7,1.};
float scale[45]={.17,.16,.16,.32,.33,1.,1.05,1.,1.,.52,.44,.44,1.05,.9,.88,.88,.88,.88,.7,.6,.6,1.05,.9,1.,1.05,1.,1.,.8,.75,.75,.85,1.05,1.05,1.05,1.05,1.05,
                .55,.45,.45,.75,.9,.95,.9,1.,1.};

struct Para
{
  //static constexpr const char *const S;
  //constexpr const *char , VTX_X="vtx_x", *VTX_Y="vtx_y", *VTX_Z="vtx_z";
  //const char *LMX="LepMomX", *LMY="LepMomY", *LMZ="LepMomZ";
  char field[20];
  const char *units;
  double l;
  double h;
  //vector<vector<vector<double>>>* field_value;
};

Para pr[]= //position is in units of cm, momentum is in units of GeV/c, angle is in units of rad, and energy is in  units of GeV
{ //match the x-ranges with the FD histograms' x-ranges
  {"vtx_x", "cm", -300., 300.},
  //{"vtx_y", "cm", -100., 100.},
  //{"vtx_z", "cm", 50., 350.},
  {"LepMomX", "GeV", -2., 2.},
  {"LepMomY", "GeV", -4.5, 2.,},
  {"LepMomZ", "GeV", -0.5, 7.},
  {"TotMom", "GeV", 0., 7.},
  {"cos_LepNuAngle", "", 0., 1.},
  {"LongMom", "GeV", -1., 7.},
  {"E_vis_true", "GeV", 0., 10.},
  {"Ev", "GeV", 0., 10.}
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

void populate_histograms(char* eff,char* caf,vector<vector<TH1D*>>& hists1,vector<vector<TH1D*>>& hists2,vector<vector<TH1D*>>& hists3,int j)
{
  TFile eff_file(eff);
  TFile caf_file(caf);
  TTree *event_data=(TTree*)eff_file.Get("event_data");
  TTree *thing=(TTree*)caf_file.Get("effTreeND");
  //TTree *cafTree=(TTree*)caf_file.Get("cafTree");
  //TTree *thing=(TTree*)caf_file.Get("effTreeND");
  //gSystem->Exec("rm -f AutoDict*vector*vector*vector*double*");
  gInterpreter->GenerateDictionary("vector<vector<vector<double>>>", "vector");
  for(auto& sel:br) event_data->SetBranchAddress(sel.eff_name, &(sel.eff_value)); //efficiencies are in 3d array, but energy is in 1d array
  thing->SetBranchAddress("ND_OffAxis_Sim_mu_start_v_xyz_LAr", &xyz_pos);
  thing->SetBranchAddress("ND_OffAxis_Sim_mu_start_p_xyz_LAr", &xyz_mom);
  thing->SetBranchAddress("ND_LepNuAngle", &lepnuangle);
  thing->SetBranchAddress("ND_Gen_numu_E", &numu_e);
  thing->SetBranchAddress("ND_E_vis_true", &e_vis_true);

  Long64_t nentries1=event_data->GetEntries();
  Long64_t nentries2=thing->GetEntries();
  if (nentries1!=nentries2) {cout<<"The efficiency file "<<eff<<" has "<<nentries2<<" events, and the CAF file "<<caf<<" has "<<nentries1<<" events."<<endl;}
  for (int i=0;i<nentries2;i++) {
    event_data->GetEntry(i);
    thing->GetEntry(i);
    unsigned long lar_pos=2;
    int k=0;
    for (Para item:pr) {
      double var_type=0.0;
      if (k==7) {var_type=numu_e;}
      if (k==8) {var_type=e_vis_true;}

      for (unsigned long vtx_pos=0;vtx_pos<NUM_VTX;vtx_pos++) {
	      int n=0;
        for (auto& sel:br) {
          TH1D* hist1=hists1[n][k];
          TH1D* hist2=hists2[n][k];
          TH1D* hist3=hists3[n][k];
          n++;
          if (vtx_pos==0||vtx_pos==1||vtx_pos==2||vtx_pos==4||vtx_pos==5||vtx_pos==16||vtx_pos==17||vtx_pos==19||vtx_pos==20||vtx_pos==21) continue;
          if (k==0) var_type=(*xyz_pos)[lar_pos][vtx_pos][0]; //I should implement a switch case block
          else if (k<4) var_type=(*xyz_mom)[lar_pos][vtx_pos][k-1];
          else if (k==4) var_type=sqrt(pow((*xyz_mom)[lar_pos][vtx_pos][0],2)+pow((*xyz_mom)[lar_pos][vtx_pos][1],2)+pow((*xyz_mom)[lar_pos][vtx_pos][2],2));
          else if (k==5) var_type=cos(lepnuangle);
          else if (k==6) var_type=sqrt(pow((*xyz_mom)[lar_pos][vtx_pos][0],2)+pow((*xyz_mom)[lar_pos][vtx_pos][1],2)+pow((*xyz_mom)[lar_pos][vtx_pos][2],2))*cos(lepnuangle);
          vector<vector<double>>* eff_value=sel.eff_value;
          vector<vector<double>>& eff_value2=*eff_value;
          double geo_eff=eff_value2[lar_pos][vtx_pos];
          if (geo_eff>1.) {cout<<"efficiency of event "<<i<<" at position "<<LAr_position[lar_pos]<<", "<<vertex_position[vtx_pos]<<" is "<<geo_eff<<endl;}
          hist1->Fill(var_type);
          if (geo_eff<=0.1) {
            continue;
          } else {
            hist2->Fill(var_type, geo_eff);
            hist3->Fill(var_type);
          }
        }
      }
      k++;
    }
  }
  eff_file.Close();
  caf_file.Close();
}

void FD_selection_cuts1()
{
  char eff[99];
  char caf[99];
  vector<vector<TH1D*>> histograms1;
  vector<vector<TH1D*>> histograms2;
  vector<vector<TH1D*>> histograms3;
  for(auto& sel:br)
  {
    const char* dt=sel.sel_name;
    vector<TH1D*> histset1, histset2, histset3;
    histograms1.push_back(histset1);
    histograms2.push_back(histset2);
    histograms3.push_back(histset3);
    int i=0;
    for (Para item:pr)
    {
      double lowerbound=item.l;
      double upperbound=item.h;
      histograms1.back().push_back(new TH1D(Form("%s_hist_%d",dt,i), Form("raw %s %d", dt, i), 200, lowerbound, upperbound));
      histograms2.back().push_back(new TH1D(Form("%s_hist_%d",dt,i), Form("selected %s %d", dt, i), 200, lowerbound, upperbound));
      histograms3.back().push_back(new TH1D(Form("%s_hist_%d",dt,i), Form("geo corrected %s %d", dt, i), 200, lowerbound, upperbound));
    i++;
    }
  }

  for (int j=1; j<8046; j++)
  {
    memset(eff, 0, 99); //clear array each time
    memset(caf, 0, 99);
    sprintf(eff, "/storage/shared/barwu/FDEff_2811722/FDGeoEff_2811722_%d_Eff.root", j);
    sprintf(caf, "/storage/shared/barwu/10thTry/FDGeoEff_2811722/FDGeoEff_2811722_%d.root", j);
    if(access(eff, 0)==0)
    {
      populate_histograms(eff,caf,histograms1,histograms2,histograms3,j);
    } else {
      cout<<"Warning: missing file:"<<eff<<endl;
      continue;
    }
  }

  gStyle->SetOptStat(000000000);
  //gStyle->SetOptStat(111111111);
  TCanvas *c=new TCanvas("c","Energy Distributions",2000,1000);
  TCanvas *r=new TCanvas("r","Ratio Plots",2000,1000);
  c->Divide(9,5);
  r->Divide(9,5);
  int i=1;
  int i_select=0;
  for (auto& sel:br)
  {
    const char *dt=sel.sel_name;
    int n=0;
    for(Para item:pr)
    {
      const char *fd=item.field;
      const char *var_unit=item.units;
      double lowerbound=item.l;
      double upperbound=item.h;
      TVirtualPad *pad=c->cd(i);
      if (n==5) {pad->SetLogy();} //pad needs to be made logarithmic, not canvas
      TH1D *hist3=histograms3[i_select][n];
      hist3->SetLineColor(kBlue);
      hist3->Draw("hist");
      TH1D *hist2=histograms2[i_select][n];
      hist2->SetLineColor(kTeal-3);
      hist2->Draw("samehist");
      TH1D *hist1=histograms1[i_select][n];
      hist1->SetLineColor(kPink);
      hist1->Draw("samehist");

      float max1=hist1->GetMaximum();
      float max2=hist2->GetMaximum();
      float max3=hist3->GetMaximum();
      float upper_y_bound=max(max(max2,max3), max1)*1.25;
      hist3->SetAxisRange(lowerbound,upperbound,"X");
      if (n!=5) {hist3->SetAxisRange(0.,upper_y_bound,"Y");}
      else {hist3->SetAxisRange(0.1,max1,"Y");}
      hist3->SetTitle(Form("%s: %s",fd,dt));
      hist3->GetXaxis()->SetTitle(Form("%s (%s)",fd,var_unit));
      hist3->GetYaxis()->SetTitle("# of events");
      TLegend *leg=new TLegend(0.1,0.77,0.4,0.9);
      leg->SetHeader("comparison"); 
      leg->AddEntry(hist1, "raw distribution");
      leg->AddEntry(hist2, "selection-cut distribution");
      leg->AddEntry(hist3, "geo corrected distribution");
      leg->Draw();
      gPad->Update();

      r->cd(i);
      TH1D *rplot1=(TH1D*)hist3->Clone();
      rplot1->Divide(hist1);
      rplot1->SetAxisRange(lowerbound,upperbound,"X");
      rplot1->SetAxisRange(0.,1.25,"Y");
      rplot1->SetLineColor(kViolet-3);
      rplot1->Draw("hist");
      TH1D *rplot2=(TH1D*)hist2->Clone();
      rplot2->Divide(hist1);
      rplot2->SetLineColor(kOrange+7);
      rplot2->Draw("samehist");
      TH1D *rplot3=(TH1D*)hist2->Clone();
      rplot3->Divide(hist3);
      rplot3->SetLineColor(kCyan);
      rplot3->Draw("samehist");
      TLegend *rleg=new TLegend(0.1,0.77,0.4,0.9);
      rleg->SetHeader("comparison"); 
      rleg->AddEntry(rplot1, "geo vs raw");
      rleg->AddEntry(rplot2, "sel vs raw");
      rleg->AddEntry(rplot3, "sel vs geo");
      rleg->Draw();
      //TPaveStats *pr;
      //pr=(TPaveStats*)hist2->GetListOfFunctions()->FindObject("stats");
      //pr->SetFillStyle(0);
      n++;
      i++;
    }
    c->Update();
    r->Update();
    i_select++;
  }
  c->SaveAs("/home/barwu/repos/MuonEffNN/images/FD_2811722_less_edge_pos_all_0.1_eff_cut_hists.png");
  r->SaveAs("/home/barwu/repos/MuonEffNN/images/FD_2811722_less_edge_pos_all_0.1_eff_cut_hists_ratios.png");

  TCanvas *cs[5];
  TCanvas *rs[5];
  i=1;
  for (auto& sel:br)
  {
    const char *dt=sel.sel_name;
    cs[i-1]=new TCanvas(Form("c%01d",i),dt,2000,1000);
    cs[i-1]->Divide(3,3);
    rs[i-1]=new TCanvas(Form("r%01d",i),dt,2000,1000);
    rs[i-1]->Divide(3,3);
    int n=0;
    for(Para& item:pr)
    {
      const char *fd=item.field;
      const char *var_unit=item.units;
      double lowerbound=item.l;
      double upperbound=item.h;
      TVirtualPad *pads=cs[i-1]->cd(n+1);
      if (n==5) {pads->SetLogy();}
      TH1D *hist3=histograms3[i-1][n];
      hist3->SetLineColor(kBlue);
      hist3->Draw("hist");
      TH1D *hist2=histograms2[i-1][n];
      hist2->SetLineColor(kTeal-3);
      hist2->Draw("samehist");
      TH1D *hist1=histograms1[i-1][n];
      hist1->SetLineColor(kPink);
      hist1->Draw("samehist");

      float max1=hist1->GetMaximum();
      float max2=hist2->GetMaximum();
      float max3=hist3->GetMaximum();
      float upper_y_bound=max(max(max2,max3), max1)*1.1;
      hist3->SetAxisRange(lowerbound,upperbound,"X");
      if (n!=5) hist3->SetAxisRange(0.,upper_y_bound,"Y");
      else hist3->SetAxisRange(10,max1,"Y");
      hist3->SetTitle(Form("%s: %s",fd,dt));
      hist3->GetXaxis()->SetTitle(Form("%s (%s)",fd,var_unit));
      hist3->GetYaxis()->SetTitle("# of events");
      TLegend *leg=new TLegend(0.1,0.75,0.33,0.9);
      leg->SetHeader("comparison"); 
      leg->AddEntry(hist1, "raw distribution");
      leg->AddEntry(hist2, "selection-cut distribution");
      leg->AddEntry(hist3, "geo corrected distribution");
      leg->Draw();
      gPad->Update();

      rs[i-1]->cd(n+1);
      TH1D *rplot1=(TH1D*)hist3->Clone();
      rplot1->Divide(hist1);
      rplot1->SetAxisRange(lowerbound,upperbound,"X");
      rplot1->SetAxisRange(0.,1.25,"Y");
      rplot1->SetLineColor(kViolet-3);
      rplot1->Draw("hist");
      TH1D *rplot2=(TH1D*)hist2->Clone();
      rplot2->Divide(hist1);
      rplot2->SetLineColor(kOrange+7);
      rplot2->Draw("samehist");
      TH1D *rplot3=(TH1D*)hist2->Clone();
      rplot3->Divide(hist3);
      rplot3->SetLineColor(kCyan);
      rplot3->Draw("samehist");
      TLegend *rleg=new TLegend(0.1,0.77,0.4,0.9);
      rleg->SetHeader("comparison"); 
      rleg->AddEntry(rplot1, "geo vs raw");
      rleg->AddEntry(rplot2, "sel vs raw");
      rleg->AddEntry(rplot3, "sel vs geo");
      rleg->Draw();
      n++;
    }
    cs[i-1]->Update();
    rs[i-1]->Update();
    //cs[i-1]->SaveAs(Form("/home/barwu/repos/MuonEffNN/images/FD_2811722_less_edge_pos_%s_0.1_eff_cut_hists.png", dt));
    //rs[i-1]->SaveAs(Form("/home/barwu/repos/MuonEffNN/images/FD_2811722_less_edge_pos_%s_0.1_eff_cut_hists_ratios.png", dt));
    i++;
  }
}