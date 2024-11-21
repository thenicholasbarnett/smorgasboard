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

void normalizeh(TH1D *h){
    double a = h->Integral();
    h->Scale(1.0/a);
}

void save_tg_2(TGraphAsymmErrors *h, TString xtitle, TString ytitle, TString hname, Int_t ptcutbin, Int_t etaslicebin, Int_t eta_or_abseta){

    // canvas
    TCanvas *c = new TCanvas();
    c->SetTitle("");
    c->SetName(hname);
    c->cd();

    // trigger 
    TString trignames[5] = {"HLT_AK4PFJet40_v8","HLT_AK4PFJet60_v8","HLT_AK4PFJet80_v8","HLT_AK4PFJet100_v8","HLT_AK4PFJet120_v8"};
    TString trigname = trignames[ptcutbin];
    // eta bin
    TString etaslices[10] = {"-6","-4","-3","-2","-1","1","2","3","4","6"};
    TString absetaslices[5] = {"1","2","3","4","6"};
    TString etaslice = "";
    if(eta_or_abseta==0){
        etaslice = etaslices[etaslicebin]+" < #eta < "+etaslices[etaslicebin+1];
    }
    if(eta_or_abseta==1){
        etaslice = absetaslices[etaslicebin]+" < |#eta| < "+absetaslices[etaslicebin];
        if(etaslicebin==0){etaslice = "|#eta| < "+absetaslices[0];}
    }

    // clone of hist
    TGraphAsymmErrors *h_c = (TGraphAsymmErrors*)h->Clone(hname);

    // markers
    h_c->SetMarkerStyle(1);
    h_c->SetMarkerColor(kBlack);
    h_c->SetLineColor(kBlack);
    h_c->Draw("ap");

    // titles
    h_c->SetTitle("");
    h_c->SetName(hname);
    h_c->GetYaxis()->SetTitle(ytitle);
    h_c->GetYaxis()->CenterTitle(true);
    h_c->GetXaxis()->SetTitle(xtitle);
    h_c->GetXaxis()->CenterTitle(true);

    // // legend
    TLegend *l = new TLegend(0.1004,0.7,0.3,0.89);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry((TObject*)0, etaslice, "");
    l->AddEntry((TObject*)0, trigname, "");
    l->Draw("same");

    // y axis limits
    h_c->SetMinimum(0.0);
    h_c->SetMaximum(1.1);

    // saving
    c->Write();
    h_c->Write();
    c->SaveAs("../plots/"+hname+".png");

    // deleting
    delete c;
    delete h_c;
    delete l;
}

void save_tg_3(TGraphAsymmErrors *h[5][9], TGraphAsymmErrors *h1[5][5]){

    // trigger 
    TString trignames[5] = {"HLT_AK4PFJet40_v8","HLT_AK4PFJet60_v8","HLT_AK4PFJet80_v8","HLT_AK4PFJet100_v8","HLT_AK4PFJet120_v8"};
    // eta bin
    TString etaslices[10] = {"-6","-4","-3","-2","-1","1","2","3","4","6"};
    TString absetaslices[5] = {"1","2","3","4","6"};

    // making nice plots for jet turn on curves by abs eta
    for(unsigned int a=0; a<5; a++){
        
        // eta slice
        TString etaslice = "";
        if(a!=0){etaslice = absetaslices[a-1]+" < |#eta| < "+absetaslices[a];}
        if(a==0){etaslice = "|#eta| < "+absetaslices[a];}

        // canvas
        TCanvas *c = new TCanvas();
        c->SetTitle("");
        c->cd();

        // clone of hist
        TGraphAsymmErrors *h_c0 = (TGraphAsymmErrors*)h1[0][a]->Clone("h_c0");
        TGraphAsymmErrors *h_c1 = (TGraphAsymmErrors*)h1[1][a]->Clone("h_c1");
        TGraphAsymmErrors *h_c2 = (TGraphAsymmErrors*)h1[2][a]->Clone("h_c2");
        TGraphAsymmErrors *h_c3 = (TGraphAsymmErrors*)h1[3][a]->Clone("h_c3");
        TGraphAsymmErrors *h_c4 = (TGraphAsymmErrors*)h1[4][a]->Clone("h_c4");

        // setting marker type
        h_c0->SetMarkerStyle(1);
        h_c1->SetMarkerStyle(1);
        h_c2->SetMarkerStyle(1);
        h_c3->SetMarkerStyle(1);
        h_c4->SetMarkerStyle(1);

        // setting marker color
        h_c0->SetMarkerColor(kGray);
        h_c1->SetMarkerColor(kRed);
        h_c2->SetMarkerColor(kBlue);
        h_c3->SetMarkerColor(kGreen);
        h_c4->SetMarkerColor(kViolet);

        // setting line color
        h_c0->SetLineColor(kGray);
        h_c1->SetLineColor(kRed);
        h_c2->SetLineColor(kBlue);
        h_c3->SetLineColor(kGreen);
        h_c4->SetLineColor(kViolet);

        // legend
        // TLegend *leg = new TLegend(0.55,0.3,0.88,0.5);
        TLegend *leg = new TLegend(0.1004,0.6,0.3,0.89);
        leg->AddEntry(h_c0,trignames[0],"l");
        leg->AddEntry(h_c1,trignames[1],"l");
        leg->AddEntry(h_c2,trignames[2],"l");
        leg->AddEntry(h_c3,trignames[3],"l");
        leg->AddEntry(h_c4,trignames[4],"l");
        leg->AddEntry((TObject*)0, etaslice, "");
        leg->SetBorderSize(0);

        // y axis limits
        h_c0->SetMinimum(0.0);
        h_c0->SetMaximum(1.1);

        // setting up the titles
        h_c0->SetTitle("");
        h_c0->GetYaxis()->SetTitle("Trigger Efficiency");
        h_c0->GetXaxis()->SetTitle("p^{leading jet}_{T} [GeV]");
        h_c0->GetYaxis()->CenterTitle(true);
        h_c0->GetXaxis()->CenterTitle(true);

        // drawing the hists onto the canvas
        h_c0->Draw("ap");
        h_c1->Draw("p");
        h_c2->Draw("p");
        h_c3->Draw("p");
        h_c4->Draw("p");
        leg->Draw("same");

        // saving
        TString hname = Form("jet_turnoncurve_absEtabin%d",a);
        c->SetName(hname);
        c->Write();
        c->SaveAs("../plots/"+hname+".png");

        // deleting
        delete c;
        delete h_c0;
        delete h_c1;
        delete h_c2;
        delete h_c3;
        delete h_c4;
        delete leg;
    }

    // making nice plots for jet turn on curves by eta
    for(unsigned int a=0; a<9; a++){
        
        // eta slice
        TString etaslice = etaslices[a]+" < #eta < "+etaslices[a+1];

        // canvas
        TCanvas *c = new TCanvas();
        c->SetTitle("");
        c->cd();

        // clone of hist
        TGraphAsymmErrors *h_c0 = (TGraphAsymmErrors*)h[0][a]->Clone("h_c0");
        TGraphAsymmErrors *h_c1 = (TGraphAsymmErrors*)h[1][a]->Clone("h_c1");
        TGraphAsymmErrors *h_c2 = (TGraphAsymmErrors*)h[2][a]->Clone("h_c2");
        TGraphAsymmErrors *h_c3 = (TGraphAsymmErrors*)h[3][a]->Clone("h_c3");
        TGraphAsymmErrors *h_c4 = (TGraphAsymmErrors*)h[4][a]->Clone("h_c4");

        // setting marker type
        h_c0->SetMarkerStyle(1);
        h_c1->SetMarkerStyle(1);
        h_c2->SetMarkerStyle(1);
        h_c3->SetMarkerStyle(1);
        h_c4->SetMarkerStyle(1);

        // setting marker color
        h_c0->SetMarkerColor(kGray);
        h_c1->SetMarkerColor(kRed);
        h_c2->SetMarkerColor(kBlue);
        h_c3->SetMarkerColor(kGreen);
        h_c4->SetMarkerColor(kViolet);

        // setting line color
        h_c0->SetLineColor(kGray);
        h_c1->SetLineColor(kRed);
        h_c2->SetLineColor(kBlue);
        h_c3->SetLineColor(kGreen);
        h_c4->SetLineColor(kViolet);

        // legend
        // TLegend *leg = new TLegend(0.55,0.3,0.88,0.5);
        TLegend *leg = new TLegend(0.1004,0.6,0.3,0.89);
        leg->AddEntry(h_c0,trignames[0],"l");
        leg->AddEntry(h_c1,trignames[1],"l");
        leg->AddEntry(h_c2,trignames[2],"l");
        leg->AddEntry(h_c3,trignames[3],"l");
        leg->AddEntry(h_c4,trignames[4],"l");
        leg->AddEntry((TObject*)0, etaslice, "");
        leg->SetBorderSize(0);

        // y axis limits
        h_c0->SetMinimum(0.0);
        h_c0->SetMaximum(1.1);

        // setting up the titles
        h_c0->SetTitle("");
        h_c0->GetYaxis()->SetTitle("Trigger Efficiency");
        h_c0->GetXaxis()->SetTitle("p^{leading jet}_{T} [GeV]");
        h_c0->GetYaxis()->CenterTitle(true);
        h_c0->GetXaxis()->CenterTitle(true);

        // drawing the hists onto the canvas
        h_c0->Draw("ap");
        h_c1->Draw("p");
        h_c2->Draw("p");
        h_c3->Draw("p");
        h_c4->Draw("p");
        leg->Draw("same");

        // saving
        TString hname = Form("jet_turnoncurve_Etabin%d",a);
        c->SetName(hname);
        c->Write();
        c->SaveAs("../plots/"+hname+".png");

        // deleting
        delete c;
        delete h_c0;
        delete h_c1;
        delete h_c2;
        delete h_c3;
        delete h_c4;
        delete leg;
    }
}

void save_tg_4(TGraphAsymmErrors *h1[5][5]){

    // trigger 
    TString trignames[5] = {"HLT_AK4PFJet40_v8","HLT_AK4PFJet60_v8","HLT_AK4PFJet80_v8","HLT_AK4PFJet100_v8","HLT_AK4PFJet120_v8"};
    // eta bin
    TString etaslices[10] = {"-6","-4","-3","-2","-1","1","2","3","4","6"};
    TString absetaslices[5] = {"1","2","3","4","6"};

    // making nice plots for jet turn on curves by abs eta
    for(unsigned int a=0; a<5; a++){
        
        // eta slice
        TString absetaslice[5] = {"","","","",""};
        for(unsigned int b=0; b<5; b++){
            if(b!=0){absetaslice[b] = absetaslices[b-1]+" < |#eta| < "+absetaslices[b];}
            if(b==0){absetaslice[b] = "|#eta| < "+absetaslices[b];}
        }
        
        // canvas
        TCanvas *c = new TCanvas();
        c->SetTitle("");
        c->cd();

        // clone of hist
        TGraphAsymmErrors *h_c0 = (TGraphAsymmErrors*)h1[a][0]->Clone("h_c0");
        TGraphAsymmErrors *h_c1 = (TGraphAsymmErrors*)h1[a][1]->Clone("h_c1");
        TGraphAsymmErrors *h_c2 = (TGraphAsymmErrors*)h1[a][2]->Clone("h_c2");
        TGraphAsymmErrors *h_c3 = (TGraphAsymmErrors*)h1[a][3]->Clone("h_c3");
        TGraphAsymmErrors *h_c4 = (TGraphAsymmErrors*)h1[a][4]->Clone("h_c4");

        // setting marker type
        h_c0->SetMarkerStyle(1);
        h_c1->SetMarkerStyle(1);
        h_c2->SetMarkerStyle(1);
        h_c3->SetMarkerStyle(1);
        h_c4->SetMarkerStyle(1);

        // setting marker color
        h_c0->SetMarkerColor(kGray);
        h_c1->SetMarkerColor(kRed);
        h_c2->SetMarkerColor(kBlue);
        h_c3->SetMarkerColor(kGreen);
        h_c4->SetMarkerColor(kViolet);

        // setting line color
        h_c0->SetLineColor(kGray);
        h_c1->SetLineColor(kRed);
        h_c2->SetLineColor(kBlue);
        h_c3->SetLineColor(kGreen);
        h_c4->SetLineColor(kViolet);

        // legend
        // TLegend *leg = new TLegend(0.55,0.3,0.88,0.5);
        TLegend *leg = new TLegend(0.1004,0.6,0.3,0.89);
        leg->AddEntry(h_c0,absetaslice[0],"l");
        leg->AddEntry(h_c1,absetaslice[1],"l");
        leg->AddEntry(h_c2,absetaslice[2],"l");
        leg->AddEntry(h_c3,absetaslice[3],"l");
        // leg->AddEntry(h_c4,absetaslice[4],"l");
        leg->AddEntry((TObject*)0, trignames[a], "");
        leg->AddEntry((TObject*)0, "run 387607", "");
        leg->SetBorderSize(0);

        // y axis limits
        h_c0->SetMinimum(0.0);
        h_c0->SetMaximum(1.1);

        // setting up the titles
        h_c0->SetTitle("");
        h_c0->GetYaxis()->SetTitle("Trigger Efficiency");
        h_c0->GetXaxis()->SetTitle("p^{leading jet}_{T} [GeV]");
        h_c0->GetYaxis()->CenterTitle(true);
        h_c0->GetXaxis()->CenterTitle(true);

        // drawing the hists onto the canvas
        h_c0->Draw("ap");
        leg->Draw("same");
        h_c1->Draw("p");
        h_c2->Draw("p");
        h_c3->Draw("p");
        // h_c4->Draw("p");

        // saving
        TString hname = Form("jet_turnoncurve_absEtabin%d_bytrig",a);
        c->SetName(hname);
        c->Write();
        c->SaveAs("../plots/"+hname+".png");

        // deleting
        delete c;
        delete h_c0;
        delete h_c1;
        delete h_c2;
        delete h_c3;
        delete h_c4;
        delete leg;
    }
}

void save_h1d_2(TH1D *h, TString xtitle, TString ytitle, TString hname, TString trigname){

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
    // h_c->SetMarkerStyle(8);
    // h_c->SetMarkerColor(kBlack);
    // h_c->SetLineColor(kBlack);

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
    TLegend *l = new TLegend(0.7,0.7,0.9,0.9);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry(h_c,hname,"pl");
    l->AddEntry((TObject*)0, trigname, "");
    l->Draw("same");

    // // y axis limits
    // h_c->SetMinimum(0.8);
    // h_c->SetMaximum(1.1);

    // saving
    c->Write();
    h->Write();
    c->SaveAs("../plots/"+hname+".png");

    // deleting
    delete c;
    // delete line1;
    delete h_c;
    delete l;
}

void save_h1d_3(TH1D *h0, TH1D *h1, TH1D *h2, TH1D *h3, int bl){
    
    // canvas names
    TString cnames[4] = {"pT_atfulleff","eta_atfulleff","phi_atfulleff","asymms_atfulleff"};
    
    // trigger names
    TString trignames[4] = {"HLT_AK4PFJet40_v8","HLT_AK4PFJet80_v8","HLT_AK4PFJet100_v8","HLT_AK4PFJet120_v8"};
    
    // titles
    TString xtitles[4] = {"p_{T}^{leading jet}","#eta^{leading jet}","#phi^{leading jet}","(p_{T}^{leading jet} - p_{T}^{subleading jet})/(p_{T}^{leading jet} + p_{T}^{subleading jet})"};

    // canvas
    TCanvas *c = new TCanvas();
    c->SetTitle("");
    c->cd();
    c->SetName(cnames[bl]);

    // hists into array
    TH1D *h_c[4];
    h_c[0] = (TH1D*)h0->Clone();
    h_c[1] = (TH1D*)h1->Clone();
    h_c[2] = (TH1D*)h2->Clone();
    h_c[3] = (TH1D*)h3->Clone();

    // legend
    TLegend *l = new TLegend(0.1004,0.6,0.4,0.89);
    l->SetBorderSize(0);
    for(unsigned int a=0; a<4; a++){
        l->AddEntry(h_c[a],trignames[a],"l");
        normalizeh(h_c[a]);
        h_c[a]->SetMarkerStyle(1);
        if(a==0){
            h_c[a]->GetXaxis()->CenterTitle(true);
            h_c[a]->GetXaxis()->SetTitle(xtitles[bl]);
            if(bl==3){
                h_c[a]->SetAxisRange(0.0,1.0,"X");
                h_c[a]->SetAxisRange(0.0,0.1,"Y");
            }
            if(bl==1){h_c[a]->SetAxisRange(0.0,0.2,"Y");}
            if(bl==2){h_c[a]->SetAxisRange(0.0,0.05,"Y");}
            h_c[a]->SetLineColor(kBlue);
            h_c[a]->SetTitle("");
            h_c[a]->Draw("e1p");
        }
        if(a==1){h_c[a]->SetLineColor(kRed);}
        if(a==2){h_c[a]->SetLineColor(kGreen);}
        if(a==3){h_c[a]->SetLineColor(kViolet);}
        if(a!=0){h_c[a]->Draw("same");}
    }
    l->Draw("same");
    c->SaveAs("../plots/"+cnames[bl]+".png");

    // deleting
    delete c;
    delete l;
}

void save_h1d_4(TH1D *h0, TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4, int bl){
    
    // canvas names
    TString cnames[4] = {"pT_atfulleff","eta_atfulleff","phi_atfulleff","asymms_atfulleff"};
    
    // trigger names
    TString trignames[5] = {"HLT_HIPuAK4CaloJet40Eta5p1_MinBiasHF1AND_v6","HLT_HIPuAK4CaloJet60Eta5p1_MinBiasHF1AND_v6","HLT_HIPuAK4CaloJet80Eta5p1_v14","HLT_HIPuAK4CaloJet100Eta5p1_v14","HLT_HIPuAK4CaloJet120Eta5p1_v14"};
    
    // titles
    TString xtitles[4] = {"p_{T}^{leading jet}","#eta^{leading jet}","#phi^{leading jet}","(p_{T}^{leading jet} - p_{T}^{subleading jet})/(p_{T}^{leading jet} + p_{T}^{subleading jet})"};

    // canvas
    TCanvas *c = new TCanvas();
    c->SetTitle("");
    c->cd();
    c->SetName(cnames[bl]);

    // hists into array
    TH1D *h_c[5];
    h_c[0] = (TH1D*)h0->Clone();
    h_c[1] = (TH1D*)h1->Clone();
    h_c[2] = (TH1D*)h2->Clone();
    h_c[3] = (TH1D*)h3->Clone();
    h_c[4] = (TH1D*)h4->Clone();

    // legend
    TLegend *l = new TLegend(0.1004,0.6,0.5,0.89);
    l->SetBorderSize(0);
    for(unsigned int a=0; a<5; a++){
        l->AddEntry(h_c[a],trignames[a],"l");
        normalizeh(h_c[a]);
        h_c[a]->SetMarkerStyle(0);
        if(a==0){
            h_c[a]->GetXaxis()->CenterTitle(true);
            h_c[a]->GetXaxis()->SetTitle(xtitles[bl]);
            if(bl==3){
                h_c[a]->SetAxisRange(0.0,1.0,"X");
                h_c[a]->SetAxisRange(0.0,0.1,"Y");
            }
            if(bl==1){h_c[a]->SetAxisRange(0.0,0.2,"Y");}
            if(bl==2){h_c[a]->SetAxisRange(0.0,0.05,"Y");}
            h_c[a]->SetLineColor(kBlack);
            h_c[a]->SetTitle("");
            h_c[a]->Draw("e1p");
            l->Draw("same");
            h_c[a]->Draw("same");
        }
        if(a==1){h_c[a]->SetLineColor(kRed);}
        if(a==2){h_c[a]->SetLineColor(kBlue);}
        if(a==3){h_c[a]->SetLineColor(kGreen);}
        if(a==4){h_c[a]->SetLineColor(kViolet);}
        if(a!=0){h_c[a]->Draw("same");}
    }
    c->SaveAs("../plots/"+cnames[bl]+".png");

    // deleting
    delete c;
    delete l;
}

void maketurnons(TGraphAsymmErrors *gr_40, TGraphAsymmErrors *gr_60, TGraphAsymmErrors *gr_80, TGraphAsymmErrors *gr_100, TGraphAsymmErrors *gr_120, int a){
    // making the canvas
    TCanvas *c = new TCanvas();
    c->SetTitle("");
    c->SetName("JetTurnOnCurve_pT");
    c->cd();

    // setting marker stuff
    // marker color
    gr_40->SetMarkerColor(kBlack);
    gr_60->SetMarkerColor(kRed);
    gr_80->SetMarkerColor(kBlue);
    gr_100->SetMarkerColor(kGreen);
    gr_120->SetMarkerColor(kViolet);
    // line color
    gr_40->SetLineColor(kBlack);
    gr_60->SetLineColor(kRed);
    gr_80->SetLineColor(kBlue);
    gr_100->SetLineColor(kGreen);
    gr_120->SetLineColor(kViolet);
    // marker type
    gr_40->SetMarkerStyle(1);
    gr_60->SetMarkerStyle(1);
    gr_80->SetMarkerStyle(1);
    gr_100->SetMarkerStyle(1);
    gr_120->SetMarkerStyle(1);

    // legend
    TLegend *leg = new TLegend(0.4,0.15,0.88,0.5);
    // TLegend *leg = new TLegend(0.35,0.3,0.88,0.5);
    // TLegend *leg = new TLegend(0.1004,0.6,0.3,0.89);
    // PbPb
    leg->AddEntry(gr_40,"HLT_HIPuAK4CaloJet40Eta5p1_MinBiasHF1AND_v6","l");
    leg->AddEntry(gr_60,"HLT_HIPuAK4CaloJet60Eta5p1_MinBiasHF1AND_v6","l");
    leg->AddEntry(gr_80,"HLT_HIPuAK4CaloJet80Eta5p1_v14","l");
    leg->AddEntry(gr_100,"HLT_HIPuAK4CaloJet100Eta5p1_v14","l");
    leg->AddEntry(gr_120,"HLT_HIPuAK4CaloJet120Eta5p1_v14","l");
    leg->AddEntry((TObject*)0, "akCs4PF Jets", "");
    // leg->AddEntry((TObject*)0, "akFlowPuCs4PF Jets", "");
    // leg->AddEntry((TObject*)0, "1.31 #mub^{-1}", "");
    leg->AddEntry((TObject*)0, "run 387973", "");
    // leg->AddEntry((TObject*)0, "PromptRecoRawPrime0", "");
    // // ppRef
    // leg->AddEntry(gr_40,"HLT_AK4PFJet40_v8","l");
    // leg->AddEntry(gr_60,"HLT_AK4PFJet60_v8","l");
    // leg->AddEntry(gr_80,"HLT_AK4PFJet80_v8","l");
    // leg->AddEntry(gr_100,"HLT_AK4PFJet100_v8","l");
    // leg->AddEntry(gr_120,"HLT_AK4PFJet120_v8","l");
    leg->SetBorderSize(0);

    // setting up the titles
    gr_40->SetTitle("");
    gr_40->GetYaxis()->SetTitle("Trigger Efficiency");
    gr_40->GetXaxis()->SetTitle("p^{leading jet}_{T} [GeV]");
    if(a==0){gr_40->GetXaxis()->SetTitle("p^{avg}_{T} [GeV]");}
    gr_40->GetYaxis()->CenterTitle(true);
    gr_40->GetXaxis()->CenterTitle(true);
    gr_40->GetYaxis()->SetRangeUser(0.0,1.1);

    // drawing everything onto the canvas
    gr_40->Draw("ap");
    gr_60->Draw("p");
    gr_80->Draw("p");
    gr_100->Draw("p");
    gr_120->Draw("p");
    // leg->AddEntry((TObject*)0, "run 387607", "");
    leg->Draw("same");
    // for(unsigned int i=0; i<5; i++){
    //     fitty[i]->Draw("same");
    // } 

    // saving canvas
    c->Draw();
    // c->SaveAs("turnonCurves_PromptRecoRawPrime0_2024PbPb_11_16_2024_run387973_0.png");
    c->SaveAs("turnonCurves_StreamerRawPrime0_2024PbPb_11_16_2024_run387973_1.png");
    c->Write();
}

void TurnOnCurves_Assemble_v4(){
    
    // getting rid of legends in hists in root file
    gStyle->SetOptStat(0);

    // file of interest
    // TFile *fi = TFile::Open("turnon_PromptRecoRawPrime0_2024PbPb_11_16_2024_run387973.root","read");
    TFile *fi = TFile::Open("turnon_StreamerRawPrime0_2024PbPb_11_16_2024_run387973.root","read");
    
    // making new root file to store plots
    // TFile *fi1 = new TFile("TurnOnCurves_PromptRecoRawPrime0_2024PbPb_11_16_2024_run387973_0.root","recreate");
    TFile *fi1 = new TFile("TurnOnCurves_StreamerRawPrime0_2024PbPb_11_16_2024_run387973_1.root","recreate");
    fi1->cd();

// GETTING HISTS //

    // pt jet turn on hists
    // by pt avg
    TH1D *denom_a = (TH1D*) fi->Get("denom_a");
    TH1D *num_40_a = (TH1D*) fi->Get("num_40_pta");
    TH1D *num_60_a = (TH1D*) fi->Get("num_60_pta");
    TH1D *num_80_a = (TH1D*) fi->Get("num_80_pta");
    TH1D *num_100_a = (TH1D*) fi->Get("num_100_pta");
    TH1D *num_120_a = (TH1D*) fi->Get("num_120_pta");
    // by leading jet
    TH1D *denom = (TH1D*) fi->Get("denom");
    TH1D *num_40 = (TH1D*) fi->Get("num_40");
    TH1D *num_60 = (TH1D*) fi->Get("num_60");
    TH1D *num_80 = (TH1D*) fi->Get("num_80");
    TH1D *num_100 = (TH1D*) fi->Get("num_100");
    TH1D *num_120 = (TH1D*) fi->Get("num_120");

// // ETA SLICED JET TURN ONS //

//     // pt jet turn on hists by eta slices
//     TH1D *denom_byAbsEta[5];
//     TH1D *denom_byEta[9];
//     TH1D *num_byAbsEta[5][5];
//     TH1D *num_byEta[5][9];
//     TGraphAsymmErrors *r_byAbsEta[5][5];
//     TGraphAsymmErrors *r_byEta[5][9];
//     for(unsigned int p=0; p<5; p++){
//         for(unsigned int e=0; e<9; e++){
//             if(e<5){
//                 TString mtitle = Form("num_byAbsEta_ptbin%d_etabin%d",p,e);
//                 num_byAbsEta[p][e] = (TH1D*) fi->Get(mtitle);
//                 num_byAbsEta[p][e]->Write();
//                 if(p==0){
//                     TString ntitle = Form("denom_byAbsEta_etabin%d",e);
//                     denom_byAbsEta[e] = (TH1D*) fi->Get(ntitle);
//                     denom_byAbsEta[e]->Write();
//                 }
//                 TString mtitle1 = Form("r_byAbsEta_ptbin%d_etabin%d",p,e);
//                 r_byAbsEta[p][e] = new TGraphAsymmErrors(num_byAbsEta[p][e],denom_byAbsEta[e],"cl=0.683 b(1,1) mode");
//                 save_tg_2(r_byAbsEta[p][e], "p_{T}^{leading jet}", "Efficiency", mtitle1, p, e, 0);
//             }
//             TString stitle = Form("num_byEta_ptbin%d_etabin%d",p,e);
//             num_byEta[p][e] = (TH1D*) fi->Get(stitle);
//             num_byEta[p][e]->Write();
//             if(p==0){
//                 TString ttitle = Form("denom_byEta_etabin%d",e);
//                 denom_byEta[e] = (TH1D*) fi->Get(ttitle);
//                 denom_byEta[e]->Write();
//             }
//             TString stitle1 = Form("r_byEta_ptbin%d_etabin%d",p,e);
//             r_byEta[p][e] = new TGraphAsymmErrors(num_byEta[p][e],denom_byEta[e],"cl=0.683 b(1,1) mode");
//             save_tg_2(r_byEta[p][e], "p_{T}^{leading jet}", "Efficiency", stitle1, p, e, 0);
//         }
//     }

//     save_tg_3(r_byEta, r_byAbsEta);
//     save_tg_4(r_byAbsEta);

// // JET KINEMATICS AT 99% EFFICIENCIES //

//     // kinematics at 99% efficiency
//     // 0 to 4 is jet40 trig to jet120 trig
//     // 0, 1, 2 is pt, eta, phi respectively
//     TH1D *kin_99[5][4];
//     // trigger names
//     TString trignames[5] = {"HLT_AK4PFJet40_v8","HLT_AK4PFJet60_v8","HLT_AK4PFJet80_v8","HLT_AK4PFJet100_v8","HLT_AK4PFJet120_v8"};
//     Int_t kin[5] = {40, 60, 80, 100, 120};
//     for(unsigned int p=0; p<5; p++){
//         for(unsigned int q=0; q<4; q++){
//             if(q==0){
//                 TString hname_0 = Form("kin99_pt_jet%d",kin[p]);
//                 kin_99[p][q] = (TH1D*) fi->Get(hname_0);
//                 save_h1d_2(kin_99[p][q], "p_{T}^{leading jet}", "", hname_0, trignames[p]);
//             }
//             if(q==1){
//                 TString hname_1 = Form("kin99_eta_jet%d",kin[p]);
//                 kin_99[p][q] = (TH1D*) fi->Get(hname_1);
//                 save_h1d_2(kin_99[p][q], "#eta^{leading jet}", "", hname_1, trignames[p]);
//             }
//             if(q==2){
//                 TString hname_2 = Form("kin99_phi_jet%d",kin[p]);
//                 kin_99[p][q] = (TH1D*) fi->Get(hname_2);
//                 save_h1d_2(kin_99[p][q], "#phi^{leading jet}", "", hname_2, trignames[p]);
//             }
//             if(q==3){
//                 TString hname_3 = Form("kin99_asymmetry_jet%d",kin[p]);
//                 kin_99[p][q] = (TH1D*) fi->Get(hname_3);
//                 save_h1d_2(kin_99[p][q], "(p_{T}^{leading jet} - p_{subleading jet})/(p_{T}^{leading jet} + p_{subleading jet})", "", hname_3, trignames[p]);
//             }
//         }
//     }
//     for(unsigned int q=0; q<4; q++){
//         // save_h1d_3(kin_99[0][q], kin_99[2][q], kin_99[3][q], kin_99[4][q], q);
//         save_h1d_4(kin_99[0][q], kin_99[1][q], kin_99[2][q], kin_99[3][q], kin_99[4][q], q);
//     }

// MAKING RATIOS //

    // th1d ratios by bayes divide
    // by leading jet
    TGraphAsymmErrors* gr_40 = new TGraphAsymmErrors(num_40,denom,"cl=0.683 b(1,1) mode");
    TGraphAsymmErrors* gr_60 = new TGraphAsymmErrors(num_60,denom,"cl=0.683 b(1,1) mode");
    TGraphAsymmErrors* gr_80 = new TGraphAsymmErrors(num_80,denom,"cl=0.683 b(1,1) mode");
    TGraphAsymmErrors* gr_100 = new TGraphAsymmErrors(num_100,denom,"cl=0.683 b(1,1) mode");
    TGraphAsymmErrors* gr_120 = new TGraphAsymmErrors(num_120,denom,"cl=0.683 b(1,1) mode");
    // by pt avg
    TGraphAsymmErrors* gr_40_a = new TGraphAsymmErrors(num_40_a,denom_a,"cl=0.683 b(1,1) mode");
    TGraphAsymmErrors* gr_60_a = new TGraphAsymmErrors(num_60_a,denom_a,"cl=0.683 b(1,1) mode");
    TGraphAsymmErrors* gr_80_a = new TGraphAsymmErrors(num_80_a,denom_a,"cl=0.683 b(1,1) mode");
    TGraphAsymmErrors* gr_100_a = new TGraphAsymmErrors(num_100_a,denom_a,"cl=0.683 b(1,1) mode");
    TGraphAsymmErrors* gr_120_a = new TGraphAsymmErrors(num_120_a,denom_a,"cl=0.683 b(1,1) mode");

// // FITTING TURN ON CURVES //
    
//     // making fitting function
//     // TF1* fitty = new TF1("fitty","1.0+[0]*TMath::Erf((x-[1])*[2])",0,300);
//     Double_t fitty_gev[5] = {40.0, 60.0, 80.0, 100.0, 120.0};
//     // Double_t fitty_p[5] = {1.0/4.0, 1.0/17.0, 1.0, 1.0, 1.0};
//     Double_t fitty_p[5] = {1.0/4.0, 1.0, 1.0, 1.0, 1.0};
//     TF1* fitty[5];
//     for(unsigned int i=0; i<5; i++){
//         // fitty[i] = new TF1("fitty","[3]+[0]*TMath::Erf((x-[1])*[2])",0,300);
//         fitty[i] = new TF1("fitty","0.5*([2]+TMath::TanH((x-[1])*[0]))",0,300);
//         fitty[i]->SetLineColor(kRed);
//         fitty[i]->SetLineStyle(0);
//         // fitty[i]->SetParLimits(0, 0.0, 0.5);
//         // fitty[i]->SetParLimits(1, 0.0, 300.0);
//         // fitty[i]->SetParLimits(0, 0.0, 1.0);
//         // fitty[i]->SetParameter(1, fitty_gev[i]);
//         // fitty[i]->SetParameter(2, fitty_p[i]);
//         // fitty[i]->FixParameter(3, 1.0-0.5*fitty_p[i]);
//         // fitty[i]->FixParameter(0,0.5*fitty_p[i]);
//         // fitty[i]->SetRange(0.0, 0.0, 300.0, 1.1);
//         fitty[i]->SetMinimum(0.0);
//         fitty[i]->SetMaximum(fitty_p[i]);
//     } 
//     // fitty[0]->SetParameter(2, 0.10);

//     // fitting function colors
//     fitty[0]->SetLineColor(kBlack);
//     fitty[1]->SetLineColor(kRed);
//     fitty[2]->SetLineColor(kBlue);
//     fitty[3]->SetLineColor(kGreen);
//     fitty[4]->SetLineColor(kViolet);

// // PRINTING FITTING INFO //

//     // fitting the hists
//     gr_40->Fit(fitty[0], "QR");
//     // gr_40_a->Fit(fitty[0], "QR");
//     cout<<fitty[0]->GetX(0.99*fitty[0]->GetMaximum(0.0,300.0),0.0,300.0)<<" [GeV] is the p_{T} that jet_40 reaches 99% effeciency"<<endl;
//     cout<<fitty[0]->GetMaximum()<<" is the jet_40 efficiency plateau"<<endl;
//     gr_60->Fit(fitty[1], "QR");
//     // gr_60_a->Fit(fitty[1], "QR");
//     cout<<fitty[1]->GetX(0.99*fitty[1]->GetMaximum(0.0,300.0),0.0,300.0)<<" [GeV] is the p_{T} that jet_60 reaches 99% effeciency"<<endl;
//     cout<<fitty[1]->GetMaximum()<<" is the jet_60 efficiency plateau"<<endl;
//     gr_80->Fit(fitty[2], "QR");
//     // gr_80_a->Fit(fitty[2], "QR");
//     cout<<fitty[2]->GetX(0.99*fitty[2]->GetMaximum(0.0,300.0),0.0,300.0)<<" [GeV] is the p_{T} that jet_80 reaches 99% effeciency"<<endl;
//     cout<<fitty[2]->GetMaximum()<<" is the jet_80 efficiency plateau"<<endl;
//     gr_100->Fit(fitty[3], "QR");
//     // gr_100_a->Fit(fitty[3], "QR");
//     cout<<fitty[3]->GetX(0.99*fitty[3]->GetMaximum(0.0,300.0),0.0,300.0)<<" [GeV] is the p_{T} that jet_100 reaches 99% effeciency"<<endl;
//     cout<<fitty[3]->GetMaximum()<<" is the jet_100 efficiency plateau"<<endl;
//     gr_120->Fit(fitty[4], "QR");
//     // gr_120_a->Fit(fitty[4], "QR");
//     cout<<fitty[4]->GetX(0.99*fitty[4]->GetMaximum(0.0,300.0),0.0,300.0)<<" [GeV] is the p_{T} that jet_120 reaches 99% effeciency"<<endl;
//     cout<<fitty[4]->GetMaximum()<<" is the jet_120 efficiency plateau"<<endl;

// MAKING PRETTY JET TURN ON PLOTS //

    maketurnons(gr_40, gr_60, gr_80, gr_100, gr_120, 1);
    // maketurnons(gr_40_a, gr_60_a, gr_80_a, gr_100_a, gr_120_a, 0);

// SAVING PLOTS //

    fi1->cd();
    // turn on curves by leading jet pt
    gr_40->Write();
    gr_60->Write();
    gr_80->Write();
    gr_100->Write();
    gr_120->Write();
    // turn on curves by pt avg
    gr_40_a->Write();
    gr_60_a->Write();
    gr_80_a->Write();
    gr_100_a->Write();
    gr_120_a->Write();
}
