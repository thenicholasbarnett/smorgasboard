// imports
#include "TROOT.h"
#include "TClass.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include <vector>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TH1.h"

// saving only in root file

void save_h1d(TH1D *h, TString hname){
    h->SetName(hname);
    h->Write();
}

void save_h2d(TH2D *h, TString hname){
    h->SetName(hname);
    h->Write();
}

void normalizeh(TH1D *h){
    double a = h->Integral();
    h->Scale(1.0/a);
}

// fitting functions

Double_t fitf1(Double_t *x, Double_t *par) {
    Double_t arg1 = 0;
    Double_t arg2 = 0;
    if (par[2]!=0) arg1 = (x[0] - par[1])/par[2];
    if (par[4]!=0) arg2 = (x[0] - par[1])/par[4];
    Double_t fitval = (par[0]*TMath::Exp(-0.5*arg1*arg1) + par[3]*TMath::Exp(-0.5*arg2*arg2));
    return fitval;
}

Double_t fitf2(Double_t *x, Double_t *par) {
    Double_t arg = 0;
    if (par[2]!=0) arg = (x[0] - par[1])/par[2];
    Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
    return fitval;
}

// saving pngs and in root file

void save_h1d_1(TH1D *h, TString xtitle, TString ytitle, TString hname){

    // canvas
    TCanvas *c = new TCanvas();
    c->SetTitle("");
    c->SetName(hname);
    c->cd();
    // c->SetLogy();

    // clone of hist
    TH1D *h_c = (TH1D*)h->Clone(hname);

    // // tline
    // int binxmax = h_c->FindLastBinAbove(h_c->GetBinLowEdge(1));
    // TLine *line1 = new TLine(h_c->GetBinLowEdge(1),1,(h_c->GetBinLowEdge(binxmax)+h_c->GetBinWidth(binxmax)),1);
    // line1->SetLineWidth(1);
    // line1->SetLineColor(kBlack);
    // line1->SetLineStyle(2);

    // markers
    h_c->Draw("e1p");
    // line1->Draw("same");
    // h_c->SetMarkerStyle(1);
    h_c->SetMarkerStyle(8);
    h_c->SetMarkerColor(kBlack);
    h_c->SetLineColor(kBlack);

    // titles
    h_c->SetTitle("");
    h_c->SetName(hname);
    h_c->GetYaxis()->SetTitle(ytitle);
    h_c->GetYaxis()->CenterTitle(true);
    h_c->GetXaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        c->SetLogy();
        h_c->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    }
    if(xtitle != "pt"){h_c->GetXaxis()->SetTitle(xtitle);}

    // // legend
    // TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    // l->SetBorderSize(0);
    // l->SetFillStyle(0);
    // l->AddEntry(h_c,hname,"pl");
    // l->Draw("same");

    // // y axis limits
    // h_c->SetMinimum(0.8);
    // h_c->SetMaximum(1.2);

    // saving
    c->Write();
    h->Write();
    c->SaveAs("plots/other/"+hname+".png");

    // deleting
    delete c;
    // delete line1;
    delete h_c;
    // delete l;
}

// saving alpha hists as pngs and in root file

void save_alpha(TH1D *h, TString etalo, TString etahi, TString ptlo, TString pthi, TString hname){

    // canvas
    TCanvas *c = new TCanvas();
    c->SetTitle("");
    c->SetName(hname);
    c->cd();

    // clone of hist
    TH1D *h_c = (TH1D*)h->Clone(hname);

    // tline
    int binxmax = h_c->FindLastBinAbove(h_c->GetBinLowEdge(1));
    TLine *line1 = new TLine(h_c->GetBinLowEdge(1),1,(h_c->GetBinLowEdge(binxmax)+h_c->GetBinWidth(binxmax)),1);
    line1->SetLineWidth(1);
    line1->SetLineColor(kBlack);
    line1->SetLineStyle(2);

    // markers
    h_c->Draw("e1p");
    line1->Draw("same");
    h_c->SetMarkerStyle(8);
    h_c->SetMarkerColor(kBlack);
    h_c->SetLineColor(kBlack);

    // titles
    h_c->SetTitle("");
    h_c->SetName(hname);
    h_c->GetXaxis()->SetTitle("#alpha");
    h_c->GetYaxis()->CenterTitle(true);
    h_c->GetXaxis()->CenterTitle(true);
    h_c->GetYaxis()->SetTitle("R^{MC}/R^{DATA}");

    // legend
    TLegend *l = new TLegend(0.55,0.65,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry((TObject*)0, etalo+" < #eta^{probe} < "+etahi, "");
    l->AddEntry((TObject*)0, ptlo+" < p_{T}^{avg} < "+pthi, "");
    l->Draw("same");

    // // y axis limits
    // h_c->SetMinimum(0.8);
    // h_c->SetMaximum(1.2);

    // saving
    h->Write();
    c->SaveAs("plots/alphas/"+hname+".png");

    // deleting
    delete c;
    delete line1;
    delete h_c;
    delete l;
}

// saving A distributions

void save_adist(TH1D *h, TString etalo, TString etahi, TString ptlo, TString pthi, TString alfa, int tagg, int num, double chi2, TString hname){

    // canvas
    TCanvas *c = new TCanvas();
    c->SetTitle("");
    c->SetName(hname);
    c->cd();

    // clone of hist
    TH1D *h_c = (TH1D*)h->Clone(hname);

    // markers
    h_c->Draw("e1p");
    h_c->SetMarkerStyle(1);
    h_c->SetMarkerColor(kBlack);
    h_c->SetLineColor(kBlack);

    // titles
    h_c->SetTitle("");
    h_c->SetName(hname);
    h_c->GetXaxis()->SetTitle("A value");
    h_c->GetYaxis()->SetTitle("Counts");
    h_c->GetXaxis()->CenterTitle(true);
    h_c->GetYaxis()->CenterTitle(true);

    // legend
    TLegend *l = new TLegend(0.55,0.65,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry((TObject*)0, etalo+" < #eta^{probe} < "+etahi, "");
    l->AddEntry((TObject*)0, ptlo+" < p_{T}^{avg} < "+pthi, "");
    l->AddEntry((TObject*)0, "#alpha < "+alfa, "");
    TString num_ = Form("%.0d",num);
    l->AddEntry((TObject*)0, "Total A count is "+num_, "");
    if(tagg==1){
        l->AddEntry((TObject*)0, "DATA", "");
    }
    if(tagg==0){
        l->AddEntry((TObject*)0, "MC", "");
    }
    TString chi2_ = Form("%0.4f",chi2);
    l->AddEntry((TObject*)0, "#chi^{2}/N_{df} = "+chi2_, "");
    // TString meann_ = Form("%0.4f",meann);
    // l->AddEntry((TObject*)0, "mean is "+meann_, "");
    l->Draw("same");

    // y axis limits
    h_c->SetMinimum(0);

    // saving
    h->Write();
    // if((num>100)&&(checkval_%20==0)){
    if(num>100){
        c->SaveAs("plots/adists/"+hname+".png");
        if(chi2>5){
            c->SaveAs("plots/adists/bad/"+hname+".png");
        }
    }

    // deleting
    delete c;
    delete h_c;
    delete l;
}

// saving final hists

void save_finals(TH1D *h, TString ptlo, TString pthi, TString hname){

    // canvas
    TCanvas *c = new TCanvas();
    c->SetTitle("");
    c->SetName(hname);
    c->cd();

    // clone of hist
    TH1D *h_c = (TH1D*)h->Clone(hname);

    // tline
    int binxmax = h_c->FindLastBinAbove(h_c->GetBinLowEdge(1));
    TLine *line1 = new TLine(h_c->GetBinLowEdge(1),1,(h_c->GetBinLowEdge(binxmax)+h_c->GetBinWidth(binxmax)),1);
    line1->SetLineWidth(1);
    line1->SetLineColor(kBlack);
    line1->SetLineStyle(2);

    // markers
    h_c->Draw("e1p");
    line1->Draw("same");
    h_c->SetMarkerStyle(1);
    // h_c->SetMarkerStyle(8);
    h_c->SetMarkerColor(kBlack);
    h_c->SetLineColor(kBlack);

    // titles
    h_c->SetTitle("");
    h_c->SetName(hname);
    h_c->GetXaxis()->SetTitle("#eta^{probe}");
    h_c->GetYaxis()->SetTitle("(R^{MC}/R^{DATA})_{#alpha-->0}");
    h_c->GetXaxis()->CenterTitle(true);
    h_c->GetYaxis()->CenterTitle(true);

    // legend
    TLegend *l = new TLegend(0.7,0.8,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry((TObject*)0, ptlo+" < p_{T}^{avg} < "+pthi, "");
    l->Draw("same");

    // // y axis limits
    h_c->SetMinimum(0.8);
    // h_c->SetMaximum(1.2);

    // saving
    h->Write();
    c->SaveAs("plots/finals/"+hname+".png");

    // deleting
    delete c;
    delete line1;
    delete h_c;
    delete l;
}

// new subjet code

void save_h1d_2(TH1D *h1, TH1D *h2, TString hname1, TString hname2, TString xtitle, TString ytitle, TString hname){

    // canvas
    TCanvas *c = new TCanvas();
    c->SetTitle("");
    c->SetName(hname);
    c->cd();
    // c->SetLogy();

    // // tline
    // TLine *line1 = new TLine(-5.2,1,5.2,1);
    // line1->SetLineWidth(1);
    // line1->SetLineColor(kBlack);
    // line1->SetLineStyle(2);
    
    // clones of hists
    TH1D *h1_c = (TH1D*)h1->Clone(hname1);
    TH1D *h2_c = (TH1D*)h2->Clone(hname2);

    // markers
    h1_c->Draw("e1p");
    // h1_c->SetMaximum(1.2);
    // h1_c->SetMinimum(0.8);
    h2_c->Draw("same");
    // line1->Draw("same");
    // h1_c->SetMarkerStyle(1);
    h1_c->SetMarkerStyle(1);
    h2_c->SetMarkerStyle(1);
    h1_c->SetMarkerColor(kBlack);
    h2_c->SetMarkerColor(kRed);
    h1_c->SetLineColor(kBlack);
    h2_c->SetLineColor(kRed);

    // titles
    h1_c->SetTitle("");
    h1_c->SetName(hname);
    h1_c->GetYaxis()->SetTitle(ytitle);
    h1_c->GetYaxis()->CenterTitle(true);
    h1_c->GetXaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        c->SetLogy();
        h1_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    }
    if(xtitle != "pt"){h1_c->GetXaxis()->SetTitle(xtitle);}
    
    // legend
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry(h1_c,hname1,"pl");
    l->AddEntry(h2_c,hname2,"pl");
    l->Draw("same");

    // saving
    c->Write();
    c->SaveAs("plots/"+hname+".png");

    // deleting
    delete c;
    // delete line1;
    delete h1_c;
    delete h2_c;
    delete l;
}

void save_h1d_panels(TH1D *h1, TH1D *h2, TH1D *h3, TString hname1, TString hname2, TString hname3, TString xtitle, TString ytitle, TString hname){
    
    // canvas
    TCanvas *c = new TCanvas();
    c->SetTitle("");
    c->SetName(hname);
    c->cd();

    // tpads
    TPad *pad1 = new TPad("pad1","pad1", 0,0.4,1,1);
    TPad *pad2 = new TPad("pad2","pad2", 0,0.0,1,0.4);
    pad1->SetLeftMargin(0.2);
    pad2->SetLeftMargin(0.2);
    pad1->SetBottomMargin(0.0);
    pad2->SetBottomMargin(0.2);
    pad1->SetTopMargin(0.1);
    pad2->SetTopMargin(0.0);
    pad1->Draw();
    pad2->Draw();

    // tline
    TLine *line1 = new TLine(-5.2,1,5.2,1);
    line1->SetLineWidth(1);
    line1->SetLineColor(kBlack);
    line1->SetLineStyle(2);

    // clones of hists
    TH1D *h1_c = (TH1D*)h1->Clone(hname1);
    TH1D *h2_c = (TH1D*)h2->Clone(hname2);
    TH1D *h3_c = (TH1D*)h3->Clone(hname3);

    // tpad1
    pad1->cd();
    // markers
    h1_c->Draw("e1p");
    h2_c->Draw("same");
    line1->Draw("same");
    // y limits
    h1_c->SetMinimum(0.8);
    h1_c->SetMaximum(1.3);
    // x axis title
    h1_c->GetYaxis()->SetTitle("R");
    h1_c->GetYaxis()->SetTitleSize(.06);
    h1_c->GetYaxis()->CenterTitle(true);
    h1_c->SetMarkerStyle(1);
    h2_c->SetMarkerStyle(1);
    h1_c->SetMarkerColor(kBlue);
    h1_c->SetLineColor(kBlue);
    h2_c->SetLineColor(kRed);
    // titles
    h1_c->SetTitle("");
    h1_c->SetName(hname);
    // legend
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry(h1_c,hname1,"pl");
    l->AddEntry(h2_c,hname2,"pl");
    l->AddEntry(h3_c,hname3,"pl");
    l->Draw("same");

    // tpad2
    pad2->cd();
    // markers
    h3_c->Draw("e1p");
    line1->Draw("same");
    h3_c->SetMarkerStyle(1);
    h3_c->SetMarkerColor(kBlack);
    h3_c->SetLineColor(kBlack);
    // y limits
    h3_c->SetMinimum(0.8);
    h3_c->SetMaximum(1.2);
    // titles
    h3_c->SetTitle("");
    h3_c->SetName(hname);
    h3_c->GetYaxis()->SetTitle(ytitle);
    h3_c->GetYaxis()->SetTitleSize(.08);
    h3_c->GetXaxis()->SetTitleSize(.08);
    h3_c->GetYaxis()->CenterTitle(true);
    h3_c->GetXaxis()->CenterTitle(true);
    if(xtitle == "pt"){
        c->SetLogy();
        h3_c->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    }
    if(xtitle != "pt"){h3_c->GetXaxis()->SetTitle(xtitle);}

    // saving
    c->Write();
    c->SaveAs("plots/panels/"+hname+".png");

    // deleting
    delete c;
    delete line1;
    delete h1_c;
    delete h2_c;
    delete h3_c;
    delete l;
}

// the script all runs in this function
void pt_balance_ratio_generator_2023ppRef_11_5()
{
    // taking into account the appropriate errors
    TH1::SetDefaultSumw2();
    
    // getting rid of legends in hists in root file
    gStyle->SetOptStat(0);

    // for generating random numbers
    TRandom2 *rand = new TRandom2(1);

    // ptslices
    // number of pt slices
    const Int_t ptslicenum = 4;
    const Int_t ptslicenum1 = 5;
    // the low and high pt values for each pt slice
    // other high pt (HP)
    double pts[ptslicenum1] = {55,80,120,170,1000};
    double ptlow[ptslicenum] = {55,80,120,170};
    double pthigh[ptslicenum] = {80,120,170,1000};
    // // low pt (ZB)
    // double pts[ptslicenum1] = {15,25,55,80,1000};
    // float ptlow[ptslicenum] = {15,25,55,80};
    // float pthigh[ptslicenum] = {25,55,80,1000};
    
    // eta slices
    // number of eta slices
    const Int_t etaslicenum = 36;
    const Int_t etaslicenum1 = 37;
    // the low and high eta values for each eta slice
    double etas[etaslicenum1] = {-5.191,-3.839,-3.489,-3.139,-2.964,-2.853,-2.65,-2.5,-2.322,-2.172,-1.93,-1.653,-1.479,-1.305,-1.044,-0.783,-0.522,-0.261,0,0.261,0.522,0.783,1.044,1.305,1.479,1.653,1.93,2.172,2.322,2.5,2.65,2.853,2.964,3.139,3.489,3.839,5.191};
    float etahigh[etaslicenum] = {-3.839,-3.489,-3.139,-2.964,-2.853,-2.65,-2.5,-2.322,-2.172,-1.93,-1.653,-1.479,-1.305,-1.044,-0.783,-0.522,-0.261,0,0.261,0.522,0.783,1.044,1.305,1.479,1.653,1.93,2.172,2.322,2.5,2.65,2.853,2.964,3.139,3.489,3.839,5.191};
    float etalow[etaslicenum] = {-5.191,-3.839,-3.489,-3.139,-2.964,-2.853,-2.65,-2.5,-2.322,-2.172,-1.93,-1.653,-1.479,-1.305,-1.044,-0.783,-0.522,-0.261,0,0.261,0.522,0.783,1.044,1.305,1.479,1.653,1.93,2.172,2.322,2.5,2.65,2.853,2.964,3.139,3.489,3.839};

    // alpha bins
    // number of alpha bins
    const Int_t alphabinnum = 9;
    const Int_t alphabinnum1 = 10;
    double alphas[alphabinnum] = {0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45};
    double alphabins[alphabinnum1] = {0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45};

    double ah1d0[3] = {100,-1,1};

    // pt balance study arrays
    // pt slices
    // MC
    // R for each pt slice
    Double_t MCptslicesR[alphabinnum][ptslicenum][etaslicenum];
    // the uncertainty or error in R for each pt slice
    Double_t MCptslicesRerr[alphabinnum][ptslicenum][etaslicenum];
    // DATA
    Double_t DATAptslicesR[alphabinnum][ptslicenum][etaslicenum];
    Double_t DATAptslicesRerr[alphabinnum][ptslicenum][etaslicenum];
    // DATA & MC
    // eta values on x axis
    Double_t ptslicesx[etaslicenum];
    // making tgraphs with information so the error in x is half the bin width
    Double_t ptslicesxerr[etaslicenum];

    // the only plots we are grabbing from the root files
    TH2D *hMCptslicesA[alphabinnum][ptslicenum];
    TH2D *hDATAptslicesA[alphabinnum][ptslicenum];

    // the slices of those plots
    TH1D *hMCetasbins_of_ptslicesA[alphabinnum][ptslicenum][etaslicenum];
    TH1D *hDATAetasbins_of_ptslicesA[alphabinnum][ptslicenum][etaslicenum];

    // alpha extrapolation
    TH1D *halpha[ptslicenum][etaslicenum];
    TH1D *halpha1[ptslicenum][etaslicenum];  

    // for mc R values 1d hists
    TH1D *hMCptslicesR[alphabinnum][ptslicenum];
    TH1D *hMCptslicesR_c[alphabinnum][ptslicenum];
    // for data R values 1d hists
    TH1D *hDATAptslicesR[alphabinnum][ptslicenum];
    TH1D *hDATAptslicesR_c[alphabinnum][ptslicenum];
    // for data/mc values 1d hists
    TH1D *hptslicesR[alphabinnum][ptslicenum];
    TH1D *hptslicesR_c[alphabinnum][ptslicenum];
    TH1D *hptslicesR_final[ptslicenum];

    // pointing fi to the file holding the jet info of interest

    TFile *fimc = TFile::Open("HP_MC_ptbalance_10_15_2024.root","read");
    TFile *fidata = TFile::Open("HP_DATA_ptbalance_10_15_2024.root","read");

    // TFile *fimc = TFile::Open("ZB_MC_ptbalance_10_15_2024.root","read");
    // TFile *fidata = TFile::Open("ZB_DATA_ptbalance_10_15_2024.root","read");

    // making an output root file to store stuff
    TFile *f1 = new TFile("pt_balance_ratios_HP_11_5_2024.root", "recreate");
    // TFile *f1 = new TFile("pt_balance_ratios_ZB_10_22_2024_SF.root", "recreate");
    f1->cd();
    
    // percent difference pt spectrum
    const double pth1d0[3] = {100,15,500};

    // gaussian fits for the jer sf values
    // 0 is MC, 1 is DATA
    TF1 *fitty[2][alphabinnum][ptslicenum][etaslicenum];

    for(unsigned int a=0; a<alphabinnum; a++){
        for(unsigned int p=0; p<ptslicenum; p++){
            // initializing some hists 
            if(a==0){
                TString bhtitle00 = Form("final_R_ptslice_%.0f_%.0f",ptlow[p],pthigh[p]);
                hptslicesR_final[p] = new TH1D(bhtitle00,"",etaslicenum,etas);
            }
            // name of hist being retrieved
            TString bhtitle0 = Form("A_ptslice_%.0f_%.0f__alpha_%.0f",ptlow[p],pthigh[p],alphas[a]*100);
            // getting the 2d hists
            hMCptslicesA[a][p] = (TH2D*) fimc->Get(bhtitle0);
            hDATAptslicesA[a][p] = (TH2D*) fidata->Get(bhtitle0);
            // save_h2d(hMCptslicesA[a][p],"MC_"+bhtitle0);
            // save_h2d(hDATAptslicesA[a][p],"DATA_"+bhtitle0);
            
            // further initializing the histograms
            for(unsigned int q=0; q<etaslicenum; q++){

                // some random numbers to check in order to not save every a distribution
                int checkval1 = rand->Integer(100);
                int checkval2 = rand->Integer(100);
                
                // important strings
                TString alfa = Form("%.2f",alphas[a]);
                TString etahi = Form("%.2f",etahigh[q]);
                TString etalo = Form("%.2f",etalow[q]);
                TString pthi = Form("%.0f",pthigh[p]);
                TString ptlo = Form("%.0f",ptlow[p]);

                // bin
                int q1 = q+1;

                if(q==0){
                    // initializing pt slice hists of R vs eta
                    TString hname01 = Form("ptslice_%.0f_%.0f__alpha_%.0f",ptlow[p],pthigh[p],alphas[a]*100);
                    hMCptslicesR[a][p] = new TH1D("MC"+hname01,"",etaslicenum,etas);
                    hDATAptslicesR[a][p] = new TH1D("DATA"+hname01,"",etaslicenum,etas);
                }

                // intializing hists that are the bins of the slices
                if(etalow[q]<0){

                    // MC
                    
                    // slicing
                    TString dhtitle0_ = Form("MC_A_ptslice_%.0f_%.0f__etabin_%.0f_%.0f_n__alpha_%.0f",ptlow[p],pthigh[p],TMath::Abs(etalow[q]*100),TMath::Abs(etahigh[q]*100),alphas[a]*100);
                    hMCetasbins_of_ptslicesA[a][p][q] = new TH1D(dhtitle0_,dhtitle0_,ah1d0[0],ah1d0[1],ah1d0[2]);
                    hMCetasbins_of_ptslicesA[a][p][q] = hMCptslicesA[a][p]->ProjectionY(dhtitle0_,q1,q1,"");

                    // Defining gaussian fit
                    fitty[0][a][p][q] = new TF1("fitty","gaus",-0.7,0.7);

                    // // Defining double gaussian fit
                    // fitty[0][a][p][q] = new TF1("fitty","([0]/[2]*exp(-0.5*pow((x-[1])/[2],2)))+([3]/[4]*exp(-0.5*pow((x-[1])/[4],2)))",-0.7,0.7);
                    // fitty[0][a][p][q]->SetParLimits(2, 0.0, 1.0);
                    // fitty[0][a][p][q]->SetParLimits(4, 0.0, 1.0);

                    // Fitting the slices
                    fitty[0][a][p][q]->SetLineColor(kRed);
                    fitty[0][a][p][q]->SetLineStyle(0);
                    if(hMCetasbins_of_ptslicesA[a][p][q]->GetEntries()>100){
                        hMCetasbins_of_ptslicesA[a][p][q]->Fit(fitty[0][a][p][q], "QR");
                        save_adist(hMCetasbins_of_ptslicesA[a][p][q], etalo, etahi, ptlo, pthi, alfa, 0, hMCetasbins_of_ptslicesA[a][p][q]->GetEntries(), fitty[0][a][p][q]->GetChisquare()/fitty[0][a][p][q]->GetNDF(), dhtitle0_);
                    }

                    // DATA
                    
                    // slicing
                    TString dhtitle1_ = Form("DATA_A_ptslice_%.0f_%.0f__etabin_%.0f_%.0f_n__alpha_%.0f",ptlow[p],pthigh[p],TMath::Abs(etalow[q]*100),TMath::Abs(etahigh[q]*100),alphas[a]*100);
                    hDATAetasbins_of_ptslicesA[a][p][q] = new TH1D(dhtitle1_,dhtitle1_,ah1d0[0],ah1d0[1],ah1d0[2]);
                    hDATAetasbins_of_ptslicesA[a][p][q] = hDATAptslicesA[a][p]->ProjectionY(dhtitle1_,q1,q1,"");
                    
                    // Defining gaussian fit
                    fitty[1][a][p][q] = new TF1("fitty","gaus",-0.7,0.7);
                    
                    // // Defining double gaussian fit
                    // fitty[1][a][p][q] = new TF1("fitty","([0]/[2]*exp(-0.5*pow((x-[1])/[2],2)))+([3]/[4]*exp(-0.5*pow((x-[1])/[4],2)))",-0.7,0.7);
                    // fitty[1][a][p][q]->SetParLimits(2, 0.0, 1.0);
                    // fitty[1][a][p][q]->SetParLimits(4, 0.0, 1.0);

                    // Fitting the slices
                    fitty[1][a][p][q]->SetLineColor(kRed);
                    fitty[1][a][p][q]->SetLineStyle(0);
                    if(hDATAetasbins_of_ptslicesA[a][p][q]->GetEntries()>100){
                        hDATAetasbins_of_ptslicesA[a][p][q]->Fit(fitty[1][a][p][q], "QR");
                        save_adist(hDATAetasbins_of_ptslicesA[a][p][q], etalo, etahi, ptlo, pthi, alfa, 1, hDATAetasbins_of_ptslicesA[a][p][q]->GetEntries(), fitty[1][a][p][q]->GetChisquare()/fitty[1][a][p][q]->GetNDF(), dhtitle1_);
                        // cout<< hDATAetasbins_of_ptslicesA[a][p][q]->GetEntries()/(hDATAetasbins_of_ptslicesA[a][p][q]->Integral()) <<" is the integral of the A hist"<<endl;
                    }
                    
                    // Initializing Alpha Extrapolation histograms
                    if(a==0){
                        TString dhtitle2_ = Form("alphas__ptslice_%.0f_%.0f__etabin_%.0f_%.0f_n",ptlow[p],pthigh[p],TMath::Abs(etalow[q]*100),TMath::Abs(etahigh[q]*100));
                        halpha[p][q] = new TH1D(dhtitle2_,dhtitle2_,alphabinnum,alphabins);
                        TString dhtitle2_1 = Form("alphas1__ptslice_%.0f_%.0f__etabin_%.0f_%.0f_n",ptlow[p],pthigh[p],TMath::Abs(etalow[q]*100),TMath::Abs(etahigh[q]*100));
                        halpha1[p][q] = new TH1D(dhtitle2_1,dhtitle2_1,alphabinnum,alphabins);
                    }
                }
                if(etalow[q]>0||etalow[q]==0){

                    // MC
                    
                    // slicing
                    TString dhtitle0 = Form("MC_A_ptslice_%.0f_%.0f__etabin_%.0f_%.0f__alpha_%.0f",ptlow[p],pthigh[p],etalow[q]*100,etahigh[q]*100,alphas[a]*100);
                    hMCetasbins_of_ptslicesA[a][p][q] = new TH1D(dhtitle0,dhtitle0,ah1d0[0],ah1d0[1],ah1d0[2]);
                    hMCetasbins_of_ptslicesA[a][p][q] = hMCptslicesA[a][p]->ProjectionY(dhtitle0,q1,q1,"");
                    
                    // Defining gaussian fit
                    fitty[0][a][p][q] = new TF1("fitty","gaus",-0.7,0.7);

                    // // Defining double gaussian fit
                    // fitty[0][a][p][q] = new TF1("fitty","([0]/[2]*exp(-0.5*pow((x-[1])/[2],2)))+([3]/[4]*exp(-0.5*pow((x-[1])/[4],2)))",-0.7,0.7);
                    // fitty[0][a][p][q]->SetParLimits(2, 0.0, 1.0);
                    // fitty[0][a][p][q]->SetParLimits(4, 0.0, 1.0);

                    // Fitting the slices
                    fitty[0][a][p][q]->SetLineColor(kRed);
                    fitty[0][a][p][q]->SetLineStyle(0);
                    if(hMCetasbins_of_ptslicesA[a][p][q]->GetEntries()>100){
                        hMCetasbins_of_ptslicesA[a][p][q]->Fit(fitty[0][a][p][q], "QR");
                        save_adist(hMCetasbins_of_ptslicesA[a][p][q], etalo, etahi, ptlo, pthi, alfa, 0, hMCetasbins_of_ptslicesA[a][p][q]->GetEntries(), fitty[0][a][p][q]->GetChisquare()/fitty[0][a][p][q]->GetNDF(), dhtitle0);
                    }

                    // DATA

                    // Slicing
                    TString dhtitle1 = Form("DATA_A_ptslice_%.0f_%.0f__etabin_%.0f_%.0f__alpha_%.0f",ptlow[p],pthigh[p],etalow[q]*100,etahigh[q]*100,alphas[a]*100);
                    hDATAetasbins_of_ptslicesA[a][p][q] = new TH1D(dhtitle1,dhtitle1,ah1d0[0],ah1d0[1],ah1d0[2]);
                    hDATAetasbins_of_ptslicesA[a][p][q] = hDATAptslicesA[a][p]->ProjectionY(dhtitle1,q1,q1,"");
                    
                    // Defining gaussian fitting
                    fitty[1][a][p][q] = new TF1("fitty","gaus",-0.7,0.7);

                    // // Defining double guassian fitting
                    // fitty[1][a][p][q] = new TF1("fitty","([0]/[2]*exp(-0.5*pow((x-[1])/[2],2)))+([3]/[4]*exp(-0.5*pow((x-[1])/[4],2)))",-0.7,0.7);
                    // fitty[1][a][p][q]->SetParLimits(2, 0.0, 1.0);
                    // fitty[1][a][p][q]->SetParLimits(4, 0.0, 1.0);

                    // Fitting the slices
                    fitty[1][a][p][q]->SetLineColor(kRed);
                    fitty[1][a][p][q]->SetLineStyle(0);
                    if(hDATAetasbins_of_ptslicesA[a][p][q]->GetEntries()>100){
                        hDATAetasbins_of_ptslicesA[a][p][q]->Fit(fitty[1][a][p][q], "QR");
                        save_adist(hDATAetasbins_of_ptslicesA[a][p][q], etalo, etahi, ptlo, pthi, alfa, 1, hDATAetasbins_of_ptslicesA[a][p][q]->GetEntries(), fitty[1][a][p][q]->GetChisquare()/fitty[1][a][p][q]->GetNDF(), dhtitle1);
                    }

                    // Initializing Alpha Extrapolation histograms
                    if(a==0){
                        TString dhtitle2_ = Form("alphas__ptslice_%.0f_%.0f__etabin_%.0f_%.0f",ptlow[p],pthigh[p],etalow[q]*100,etahigh[q]*100);
                        halpha[p][q] = new TH1D(dhtitle2_,dhtitle2_,alphabinnum,alphabins);
                        TString dhtitle2_1 = Form("alphas1__ptslice_%.0f_%.0f__etabin_%.0f_%.0f",ptlow[p],pthigh[p],etalow[q]*100,etahigh[q]*100);
                        halpha1[p][q] = new TH1D(dhtitle2_1,dhtitle2_1,alphabinnum,alphabins);
                    }
                }
                
                // finding stuff for scale factors

                // MC
                
                // // for double gaussian
                // Double_t x_mc = fitty[0][a][p][q]->GetParameter(2);
                // Double_t dx_mc = fitty[0][a][p][q]->GetParError(2);
                // Double_t y_mc = fitty[0][a][p][q]->GetParameter(4);
                // Double_t dy_mc = fitty[0][a][p][q]->GetParError(4);
                // Double_t MCptaavg = TMath::Sqrt(x_mc*x_mc + y_mc*y_mc);
                // Double_t MCptaavgerr = TMath::Sqrt(x_mc*x_mc*dx_mc*dx_mc/(x_mc*x_mc+y_mc*y_mc) + y_mc*y_mc*dy_mc*dy_mc/(x_mc*x_mc+y_mc*y_mc));
                
                // // for gaussian
                // Double_t MCptaavg = fitty[0][a][p][q]->GetParameter(2);
                // Double_t MCptaavgerr = fitty[0][a][p][q]->GetParError(2);
                
                // if(fitty[0][a][p][q]->GetChisquare()/fitty[0][a][p][q]->GetNDF() > 5){
                //     MCptaavg=0;
                //     MCptaavgerr=0;
                // }

                // DATA

                // // for double gaussian
                // Double_t x_ = fitty[1][a][p][q]->GetParameter(2);
                // Double_t dx_ = fitty[1][a][p][q]->GetParError(2);
                // Double_t y_ = fitty[1][a][p][q]->GetParameter(4);
                // Double_t dy_ = fitty[1][a][p][q]->GetParError(4);
                // Double_t DATAptaavg = TMath::Sqrt(x_*x_ + y_*y_);
                // Double_t DATAptaavgerr = TMath::Sqrt(x_*x_*dx_*dx_/(x_*x_+y_*y_) + y_*y_*dy_*dy_/(x_*x_+y_*y_));
                
                // // for gaussian
                // Double_t DATAptaavg = fitty[1][a][p][q]->GetParameter(2);
                // Double_t DATAptaavgerr = fitty[1][a][p][q]->GetParError(2);
                
                // if(fitty[1][a][p][q]->GetChisquare()/fitty[1][a][p][q]->GetNDF() > 5){
                //     DATAptaavg=0;
                //     DATAptaavgerr=0;
                // }
                
                // finding stuff for tgraphs
                // MC
                Double_t MCptaavg = hMCetasbins_of_ptslicesA[a][p][q]->GetMean();
                Double_t MCptaavgerr = hMCetasbins_of_ptslicesA[a][p][q]->GetMeanError();
                MCptslicesR[a][p][q] = ((1+MCptaavg)/(1-MCptaavg)); 
                MCptslicesRerr[a][p][q] = (MCptaavgerr*2/((1-MCptaavg)*(1-MCptaavg)));
                // DATA
                Double_t DATAptaavg = hDATAetasbins_of_ptslicesA[a][p][q]->GetMean();
                Double_t DATAptaavgerr = hDATAetasbins_of_ptslicesA[a][p][q]->GetMeanError();
                
                // setting bin content to these values
                hMCptslicesR[a][p]->SetBinContent(q1, ((1+MCptaavg)/(1-MCptaavg)));
                hMCptslicesR[a][p]->SetBinError(q1, (MCptaavgerr*2/((1-MCptaavg)*(1-MCptaavg))));
                hDATAptslicesR[a][p]->SetBinContent(q1, ((1+DATAptaavg)/(1-DATAptaavg)));
                hDATAptslicesR[a][p]->SetBinError(q1, (DATAptaavgerr*2/((1-DATAptaavg)*(1-DATAptaavg))));
            }
        }
    }

    // making the ratios
    for(unsigned int a=0; a<alphabinnum; a++){
        Int_t a1 = a+1;
        // pt slices
        for(unsigned int p=0; p<ptslicenum; p++){
            TString hnamep = Form("ptslice_%.0f_%.0f__alpha_%.0f",ptlow[p],pthigh[p],alphas[a]*100);
            // making the ratio of the Rval hists (R_MC / R_DATA)
            hMCptslicesR_c[a][p] = (TH1D*)hMCptslicesR[a][p]->Clone();
            hDATAptslicesR_c[a][p] = (TH1D*)hDATAptslicesR[a][p]->Clone();
            hptslicesR[a][p] = (TH1D*)hMCptslicesR_c[a][p]->Clone();
            hptslicesR[a][p]->Divide(hMCptslicesR_c[a][p],hDATAptslicesR_c[a][p],1,1,"");
            for(unsigned int q=0; q<etaslicenum; q++){
                Int_t q1 = q+1;
                Double_t DATApta_entries = hDATAetasbins_of_ptslicesA[a][p][q]->GetEntries();
                Double_t MCpta_entries = hMCetasbins_of_ptslicesA[a][p][q]->GetEntries();
                int minentriesperbin = 100;
                if((DATApta_entries<minentriesperbin)||(MCpta_entries<minentriesperbin)){
                    hptslicesR[a][p]->SetBinContent(q1, 0);
                    hptslicesR[a][p]->SetBinError(q1, 0);
                }
                if(DATApta_entries<minentriesperbin){
                    hDATAptslicesR[a][p]->SetBinContent(q1, 0);
                    hDATAptslicesR[a][p]->SetBinError(q1, 0);
                }
                if(MCpta_entries<minentriesperbin){
                    hMCptslicesR[a][p]->SetBinContent(q1, 0);
                    hMCptslicesR[a][p]->SetBinError(q1, 0);
                }

                // Alpha Extrapolation
                halpha[p][q]->SetBinContent(a1, hptslicesR[a][p]->GetBinContent(q1));
                halpha[p][q]->SetBinError(a1, hptslicesR[a][p]->GetBinError(q1));
            }
            hptslicesR_c[a][p] = (TH1D*)hptslicesR[a][p]->Clone();
            // hMCptslicesR[p]->SetMaximum(1.1);
            // hMCptslicesR[p]->SetMinimum(0.9);
            // hDATAptslicesR[p]->SetMaximum(1.1);
            // hDATAptslicesR[p]->SetMinimum(0.9);
            // hptslicesR[a][p]->SetMaximum(1.1);
            // hptslicesR[a][p]->SetMinimum(0.9);
            // // saving the Rval hists
            save_h1d_1(hMCptslicesR[a][p], "#eta^{Probe}", "R^{MC}", "R_MC_"+hnamep);
            save_h1d_1(hDATAptslicesR[a][p], "#eta^{Probe}", "R^{DATA}", "R_DATA_"+hnamep);
            // saving Rval ratio hist
            // save_h1d_1(hptslicesR[a][p], "#eta^{Probe}", "R^{MC}/R^{DATA}", hnamep);
            save_h1d_panels(hMCptslicesR[a][p], hDATAptslicesR[a][p], hptslicesR[a][p], "MC", "DATA", "MC/DATA", "#eta^{Probe}", "R^{MC}/R^{DATA}", hnamep+"_panel");
            
        }
    }
    
    for(unsigned int a=0; a<alphabinnum; a++){
        Int_t a1 = a+1;
        for(unsigned int p=0; p<ptslicenum; p++){
            for(unsigned int q=0; q<etaslicenum; q++){

                // // not considering alpha bins with only 2 jet events
                // if(((halpha[p][q]->GetBinError(a1))==(halpha[p][q]->GetBinError(a)))&&((halpha[p][q]->GetBinError(a1))==(halpha[p][q]->GetBinError(a)))){continue;}

                // setting bin values
                halpha1[p][q]->SetBinContent(a1, halpha[p][q]->GetBinContent(a1));
                halpha1[p][q]->SetBinError(a1, halpha[p][q]->GetBinError(a1));
            }
        }
    }

    TF1 *fit1[ptslicenum][etaslicenum];
    
    for(unsigned int p=0; p<ptslicenum; p++){
        for(unsigned int q=0; q<etaslicenum; q++){
            Int_t q1 = q+1;
            fit1[p][q] = new TF1("fit1","[0] + x*[1]",0.0,0.45);
            fit1[p][q]->SetLineColor(kRed);
            fit1[p][q]->SetLineStyle(0);
            halpha1[p][q]->Fit(fit1[p][q]);
            hptslicesR_final[p]->SetBinContent(q1, fit1[p][q]->GetParameter(0));
            hptslicesR_final[p]->SetBinError(q1, fit1[p][q]->GetParError(0));
        }
    }
    
    for(unsigned int p=0; p<ptslicenum; p++){
        TString ptlo = Form("%.0f",ptlow[p]);
        TString pthi = Form("%.0f",pthigh[p]);
        TString dhtitle20 = Form("final_R_ptslice_%.0f_%.0f",ptlow[p],pthigh[p]);
        // save_h1d_1(hptslicesR_final[p], "#eta^{probe}", "(R^{MC}/R^{DATA})_{#alpha-->0}", dhtitle20);
        save_finals(hptslicesR_final[p], ptlo, pthi, dhtitle20);
        TString yaxistitle1 = Form("((R^{MC}/R^{DATA}) / (R^{MC}/R^{DATA})_{#alpha<%.2f})_{#alpha-->%.0f}",alphas[p],alphas[3]);
        for(unsigned int q=0; q<etaslicenum; q++){
            // halpha1[p][q]->SetMinimum(0.8);
            // halpha1[p][q]->SetMaximum(1.2);
            TString etalo = Form("%.2f",etalow[q]);
            TString etahi = Form("%.2f",etahigh[q]);
            if(etalow[q]<0){
                TString dhtitle21 = Form("alphas__ptslice_%.0f_%.0f__etabin_%.0f_%.0f_n",ptlow[p],pthigh[p],TMath::Abs(etalow[q]*100),TMath::Abs(etahigh[q]*100));
                // save_h1d_1(halpha1[p][q], "#alpha", "R^{MC}/R^{DATA}", dhtitle21);
                save_alpha(halpha1[p][q], etalo, etahi, ptlo, pthi, dhtitle21);
            }
            if(etalow[q]>0||etalow[q]==0){
                TString dhtitle22 = Form("alphas__ptslice_%.0f_%.0f__etabin_%.0f_%.0f",ptlow[p],pthigh[p],etalow[q]*100,etahigh[q]*100);
                // save_h1d_1(halpha1[p][q], "#alpha", "R^{MC}/R^{DATA}", dhtitle22);
                save_alpha(halpha1[p][q], etalo, etahi, ptlo, pthi, dhtitle22);
                
            }
        }
    }
    f1->Close();
}