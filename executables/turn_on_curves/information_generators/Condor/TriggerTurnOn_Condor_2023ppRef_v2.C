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

#include "/afs/cern.ch/user/n/nbarnett/public/header_files/JetUncertainty.h"
#include "/afs/cern.ch/user/n/nbarnett/public/header_files/JetCorrector.h"

// executed code 
void TriggerTurnOn_Condor_2023ppRef_v2(TString input, TString output){

    // making parameters for hists of interest
    const double pth1d0[3] = {30,0.0,300.0};
    const double etah1d0[3] = {25,-5.0,5.0};
    const double phih1d0[3] = {100,-1*TMath::Pi(),TMath::Pi()};
    const double ah1d0[3] = {100,-1,1};

    // making hists of interest 

    // 2d hists of pt vs eta
    TH2D *denom_pt_eta = new TH2D("denom_pt_eta","denom_pt_eta",pth1d0[0],pth1d0[1],pth1d0[2],etah1d0[0],etah1d0[1],etah1d0[2]);
    TH2D *num_40_pt_eta = new TH2D("num_40_pt_eta","num_40_pt_eta",pth1d0[0],pth1d0[1],pth1d0[2],etah1d0[0],etah1d0[1],etah1d0[2]);
    TH2D *num_60_pt_eta = new TH2D("num_60_pt_eta","num_60_pt_eta",pth1d0[0],pth1d0[1],pth1d0[2],etah1d0[0],etah1d0[1],etah1d0[2]);
    TH2D *num_80_pt_eta = new TH2D("num_80_pt_eta","num_80_pt_eta",pth1d0[0],pth1d0[1],pth1d0[2],etah1d0[0],etah1d0[1],etah1d0[2]);
    TH2D *num_100_pt_eta = new TH2D("num_100_pt_eta","num_100_pt_eta",pth1d0[0],pth1d0[1],pth1d0[2],etah1d0[0],etah1d0[1],etah1d0[2]);
    TH2D *num_120_pt_eta = new TH2D("num_120_pt_eta","num_120_pt_eta",pth1d0[0],pth1d0[1],pth1d0[2],etah1d0[0],etah1d0[1],etah1d0[2]);

    // denominator for ratios
    TH1D *denom = new TH1D("denom","denom",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *denom_a = new TH1D("denom_a","denom_a",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *denomEta = new TH1D("denomEta","denomEta",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *denomPhi = new TH1D("denomPhi","denomPhi",phih1d0[0],phih1d0[1],phih1d0[2]);
    
    // jet pt turn on curves for different eta slices
    Double_t pTtrig[5] = {40.0,60.0,80.0,100.0,120.0};
    Double_t byAbsEta[6] = {0.0,1.0,2.0,3.0,4.0,6.0};
    Double_t byEta[10] = {-6,-4,-3,-2,-1,1,2,3,4,6};
    TH1D *denom_byAbsEta[5];
    TH1D *denom_byEta[9];
    TH1D *num_byAbsEta[5][5];
    TH1D *num_byEta[5][9];
    for(unsigned int p=0; p<5; p++){
        for(unsigned int e=0; e<9; e++){
            if(e<5){
                TString mtitle = Form("num_byAbsEta_ptbin%d_etabin%d",p,e);
                num_byAbsEta[p][e] = new TH1D(mtitle,mtitle,pth1d0[0],pth1d0[1],pth1d0[2]);
                if(p==0){
                    TString ntitle = Form("denom_byAbsEta_etabin%d",e);
                    denom_byAbsEta[e] = new TH1D(ntitle,ntitle,pth1d0[0],pth1d0[1],pth1d0[2]);
                }
            }
            TString stitle = Form("num_byEta_ptbin%d_etabin%d",p,e);
            num_byEta[p][e] = new TH1D(stitle,stitle,pth1d0[0],pth1d0[1],pth1d0[2]);
            if(p==0){
                TString ttitle = Form("denom_byEta_etabin%d",e);
                denom_byEta[e] = new TH1D(ttitle,ttitle,pth1d0[0],pth1d0[1],pth1d0[2]);
            }
        }
    }

    // kinematics at 99% efficiency
    // 0 to 4 is jet40 trig to jet120 trig
    // 0, 1, 2 is pt, eta, phi respectively
    TH1D *kin_99[5][4];
    Int_t kin[5] = {40, 60, 80, 100, 120};
    for(unsigned int p=0; p<5; p++){
        for(unsigned int q=0; q<4; q++){
            if(q==0){
                TString hname_0 = Form("kin99_pt_jet%d",kin[p]);
                kin_99[p][q] = new TH1D(hname_0,hname_0,pth1d0[0],pth1d0[1],pth1d0[2]);
            }
            if(q==1){
                TString hname_1 = Form("kin99_eta_jet%d",kin[p]);
                kin_99[p][q] = new TH1D(hname_1,hname_1,etah1d0[0],etah1d0[1],etah1d0[2]);
            }
            if(q==2){
                TString hname_2 = Form("kin99_phi_jet%d",kin[p]);
                kin_99[p][q] = new TH1D(hname_2,hname_2,phih1d0[0],phih1d0[1],phih1d0[2]);
            }
            if(q==3){
                TString hname_3 = Form("kin99_asymmetry_jet%d",kin[p]);
                kin_99[p][q] = new TH1D(hname_3,hname_3,ah1d0[0],ah1d0[1],ah1d0[2]);
            }
        }
    }

    // eta analysis
    TH1D *denomEta_40 = new TH1D("denomEta_40","denomEta_40",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *denomEta_60 = new TH1D("denomEta_60","denomEta_60",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *denomEta_80 = new TH1D("denomEta_80","denomEta_80",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *denomEta_100 = new TH1D("denomEta_100","denomEta_100",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *denomEta_120 = new TH1D("denomEta_120","denomEta_120",etah1d0[0],etah1d0[1],etah1d0[2]);

    // jets > 40 
    // pT
    TH1D *num_40 = new TH1D("num_40","num_40",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *num_40_pta = new TH1D("num_40_pta","num_40_pta",pth1d0[0],pth1d0[1],pth1d0[2]);
    // eta
    TH1D *numEta_40 = new TH1D("numEta_40","numEta_40",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *feEta_40 = new TH1D("feEta_40","feEta_40",etah1d0[0],etah1d0[1],etah1d0[2]);
    // phi
    TH1D *numPhi_40 = new TH1D("numPhi_40","numPhi_40",phih1d0[0],phih1d0[1],phih1d0[2]);
    TH1D *fePhi_40 = new TH1D("fePhi_40","fePhi_40",phih1d0[0],phih1d0[1],phih1d0[2]);

    // jets > 60 
    // pT
    TH1D *num_60 = new TH1D("num_60","num_60",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *num_60_pta = new TH1D("num_60_pta","num_60_pta",pth1d0[0],pth1d0[1],pth1d0[2]);
    // eta
    TH1D *numEta_60 = new TH1D("numEta_60","numEta_60",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *feEta_60 = new TH1D("feEta_60","feEta_60",etah1d0[0],etah1d0[1],etah1d0[2]);
    // phi
    TH1D *numPhi_60 = new TH1D("numPhi_60","numPhi_60",phih1d0[0],phih1d0[1],phih1d0[2]);
    TH1D *fePhi_60 = new TH1D("fePhi_60","fePhi_60",phih1d0[0],phih1d0[1],phih1d0[2]);

    // jets > 80 
    // pT
    TH1D *num_80 = new TH1D("num_80","num_80",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *num_80_pta = new TH1D("num_80_pta","num_80_pta",pth1d0[0],pth1d0[1],pth1d0[2]);
    // eta
    TH1D *numEta_80 = new TH1D("numEta_80","numEta_80",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *feEta_80 = new TH1D("feEta_80","feEta_80",etah1d0[0],etah1d0[1],etah1d0[2]);
    // phi
    TH1D *numPhi_80 = new TH1D("numPhi_80","numPhi_80",phih1d0[0],phih1d0[1],phih1d0[2]);
    TH1D *fePhi_80 = new TH1D("fePhi_80","fePhi_80",phih1d0[0],phih1d0[1],phih1d0[2]);

    // jets > 100 
    // pt jet turn on
    TH1D *num_100 = new TH1D("num_100","num_100",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *num_100_pta = new TH1D("num_100_pta","num_100_pta",pth1d0[0],pth1d0[1],pth1d0[2]);
    // eta jet turn on info
    TH1D *numEta_100 = new TH1D("numEta_100","numEta_100",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *feEta_100 = new TH1D("feEta_100","feEta_100",etah1d0[0],etah1d0[1],etah1d0[2]);
    // phi jet turn on info
    TH1D *numPhi_100 = new TH1D("numPhi_100","numPhi_100",phih1d0[0],phih1d0[1],phih1d0[2]);
    TH1D *fePhi_100 = new TH1D("fePhi_100","fePhi_100",phih1d0[0],phih1d0[1],phih1d0[2]);

    // jets > 120 
    // pT
    TH1D *num_120 = new TH1D("num_120","num_120",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *num_120_pta = new TH1D("num_120_pta","num_120_pta",pth1d0[0],pth1d0[1],pth1d0[2]);
    // eta
    TH1D *numEta_120 = new TH1D("numEta_120","numEta_120",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *feEta_120 = new TH1D("feEta_120","feEta_120",etah1d0[0],etah1d0[1],etah1d0[2]);
    // phi
    TH1D *numPhi_120 = new TH1D("numPhi_120","numPhi_120",phih1d0[0],phih1d0[1],phih1d0[2]);
    TH1D *fePhi_120 = new TH1D("fePhi_120","fePhi_120",phih1d0[0],phih1d0[1],phih1d0[2]);
        
    // jet numbers
    Int_t njets=0;
    Int_t njets_=0;
    Int_t njets_0=0;
    Int_t njets_20=0;
    Int_t njets_40=0;
    Int_t njets_60=0;
    Int_t njets_80=0;
    Int_t njets_100=0;
    Int_t njets_120=0;

    // reading input file
    TFile *fi_input = TFile::Open(input,"read");
    fi_input->cd();

    // making the variables to get from the ttrees
    // event variables
    Int_t nref;
    Float_t vz;
    Int_t ppvf;
    // jet hlt hits
    Int_t hlt_40;
    Int_t hlt_60;
    Int_t hlt_80;
    Int_t hlt_100;
    Int_t hlt_120;
    // zero bias hlt hits
    Int_t hlt_zb;
    // jet momenta
    const unsigned int mj = 10000;
    Float_t jtpt[mj];
    Float_t jteta[mj];
    Float_t jtphi[mj];

    // getting the ttrees of interest
    // triggers
    TTree *hlt_tree = (TTree*)fi_input->Get("hltanalysis/HltTree");
    // jet momentum
    TTree *jet_tree = (TTree*)fi_input->Get("ak4PFJetAnalyzer/t");
    // vertex position
    TTree *hievt_tree = (TTree*)fi_input->Get("hiEvtAnalyzer/HiTree");
    // vertex filter
    TTree *pvf_tree = (TTree*)fi_input->Get("skimanalysis/HltTree");

    // turning all tree branches off
    hlt_tree->SetBranchStatus("*", 0);
    jet_tree->SetBranchStatus("*", 0);
    hievt_tree->SetBranchStatus("*", 0);
    pvf_tree->SetBranchStatus("*", 0);

    // turning tree branches of interest on

    // event trees
    jet_tree->SetBranchStatus("nref", 1);
    hievt_tree->SetBranchStatus("vz", 1);
    pvf_tree->SetBranchStatus("pprimaryVertexFilter", 1);
    // jet hlt hits
    hlt_tree->SetBranchStatus("HLT_AK4PFJet40_v1", 1);
    hlt_tree->SetBranchStatus("HLT_AK4PFJet60_v1", 1);
    hlt_tree->SetBranchStatus("HLT_AK4PFJet80_v1", 1);
    hlt_tree->SetBranchStatus("HLT_AK4PFJet100_v1", 1);
    hlt_tree->SetBranchStatus("HLT_AK4PFJet120_v1", 1);
    // zb hlt hits
    hlt_tree->SetBranchStatus("HLT_PPRefZeroBias_v1", 1);
    // jet momenta
    jet_tree->SetBranchStatus("jtpt", 1);
    jet_tree->SetBranchStatus("jteta", 1);
    jet_tree->SetBranchStatus("jtphi", 1);

    // turning tree branches of interest on
    // jet hlts
    hlt_tree->SetBranchAddress("HLT_AK4PFJet40_v1", &hlt_40);
    hlt_tree->SetBranchAddress("HLT_AK4PFJet60_v1", &hlt_60);
    hlt_tree->SetBranchAddress("HLT_AK4PFJet80_v1", &hlt_80);
    hlt_tree->SetBranchAddress("HLT_AK4PFJet100_v1", &hlt_100);
    hlt_tree->SetBranchAddress("HLT_AK4PFJet120_v1", &hlt_120);
    // zb hlt hit
    hlt_tree->SetBranchAddress("HLT_PPRefZeroBias_v1", &hlt_zb);
    // jet momenta
    jet_tree->SetBranchAddress("jtpt", jtpt);
    jet_tree->SetBranchAddress("jteta", jteta);
    jet_tree->SetBranchAddress("jtphi", jtphi);
    // event variables
    jet_tree->SetBranchAddress("nref", &nref);
    hievt_tree->SetBranchAddress("vz", &vz);
    pvf_tree->SetBranchAddress("pprimaryVertexFilter", &ppvf);

    // event loop
    for(unsigned int i=0; i<jet_tree->GetEntries(); i++){
    // for(unsigned int i=0; i<1000; i++){

        njets+=1;
        
        // event filtering
        pvf_tree->GetEntry(i);
        if(ppvf!=1){continue;}
        hievt_tree->GetEntry(i);
        if(TMath::Abs(vz)>15){continue;}

        njets_+=1;

        jet_tree->GetEntry(i);
        hlt_tree->GetEntry(i);

        if(nref==0){continue;}

        // getting the corrected jet pt
        vector<string> Files;
        Files.push_back("/afs/cern.ch/user/n/nbarnett/public/txt_files/L2L3_ppReco_2023ppRef/ParallelMC_L2Relative_AK4PF_pp_Reco_v0_12-21-2023.txt");
        JetCorrector JEC(Files);

        Float_t jtcorrpt[nref];
        
        // getting leading jet mom in event
        for(Int_t j = 0; j<nref; j++){
            JEC.SetJetPT(jtpt[j]);
            JEC.SetJetEta(jteta[j]);
            JEC.SetJetPhi(jtphi[j]);  
            Float_t jet_pt_corr = JEC.GetCorrectedPT();

            // saving the corrected jet pt
            jtcorrpt[j] = jet_pt_corr;
        }

        Int_t lj = 0;
        Int_t slj = 1;

        // getting leading jet 
        for(Int_t j = 0; j<nref; j++){
            if(jtcorrpt[j]>jtcorrpt[lj]){
                lj=j;
                if(lj==1){slj=0;}
            }
        }

        // getting subleading jet index
        for(Int_t j = 0; j<nref; j++){
            if((jtcorrpt[j]>jtcorrpt[slj])&&(j!=lj)){
                slj=j;
            }
        }

        // threshholds for the jet hlts when they reach 99% efficiency
        // // 2024 ppRef
        // Double_t feff_jethlts_thresh[5] = {58.5197,9999,99.2784,115.987,141.701};
        // 2024 PbPb
        Double_t feff_jethlts_thresh[5] = {9999,9999,9999,9999,9999};

        // filling hists with leading jet info 
        if((jtcorrpt[lj]>0)&&(hlt_zb==1)){
            Double_t jet_asymm = (jtcorrpt[lj] - jtcorrpt[slj])/(jtcorrpt[lj] + jtcorrpt[slj]);
            Double_t jtpt_a = (jtcorrpt[lj] + jtcorrpt[slj])/2;
            denom->Fill(jtcorrpt[lj]);
            denom_a->Fill(jtpt_a);
            denomEta->Fill(jteta[lj]);
            denomPhi->Fill(jtphi[lj]);
            denom_pt_eta->Fill(jtcorrpt[lj],jteta[lj]);

            for(unsigned int q=0; q<9; q++){
                if(q<5){
                    if(((TMath::Abs(jteta[lj])==byAbsEta[q])||(TMath::Abs(jteta[lj])>byAbsEta[q]))&&(TMath::Abs(jteta[lj])<byAbsEta[q+1])){
                        denom_byAbsEta[q]->Fill(jtcorrpt[lj]);
                    }
                }
                if(((jteta[lj]==byEta[q])||(jteta[lj]>byEta[q]))&&(jteta[lj]<byEta[q+1])){
                    denom_byEta[q]->Fill(jtcorrpt[lj]);
                }
            }

            njets_0+=1;
            if(jtcorrpt[lj]>20){njets_20+=1;}
            if(jtcorrpt[lj]>40){
                njets_40+=1;
                denomEta_40->Fill(jteta[lj]);
            }
            if(jtcorrpt[lj]>60){
                njets_60+=1;
                denomEta_60->Fill(jteta[lj]);
            }
            if(jtcorrpt[lj]>80){
                njets_80+=1;
                denomEta_80->Fill(jteta[lj]);
            }
            if(jtcorrpt[lj]>100){
                njets_100+=1;
                denomEta_100->Fill(jteta[lj]);
            }
            if(jtcorrpt[lj]>120){
                njets_120+=1;
                denomEta_120->Fill(jteta[lj]);
            }

            // if the triggers are on
            if(hlt_40==1){
                // filling hists
                num_40->Fill(jtcorrpt[lj]);
                num_40_pta->Fill(jtpt_a);
                numEta_40->Fill(jteta[lj]);
                numPhi_40->Fill(jtphi[lj]);
                if(jtcorrpt[lj]>40){
                    feEta_40->Fill(jteta[lj]);
                    fePhi_40->Fill(jtphi[lj]);
                    if(jtcorrpt[lj]>feff_jethlts_thresh[0]){
                        kin_99[0][0]->Fill(jtcorrpt[lj]);
                        kin_99[0][1]->Fill(jteta[lj]);
                        kin_99[0][2]->Fill(jtphi[lj]);
                        kin_99[0][3]->Fill(jet_asymm);
                    }
                }
                num_40_pt_eta->Fill(jtcorrpt[lj],jteta[lj]);
                // pt jet turn on curve by eta bins
                for(unsigned int q=0; q<9; q++){
                    if(q<5){
                        if(((TMath::Abs(jteta[lj])==byAbsEta[q])||(TMath::Abs(jteta[lj])>byAbsEta[q]))&&(TMath::Abs(jteta[lj])<byAbsEta[q+1])){
                            num_byAbsEta[0][q]->Fill(jtcorrpt[lj]);
                        }
                    }
                    if(((jteta[lj]==byEta[q])||(jteta[lj]>byEta[q]))&&(jteta[lj]<byEta[q+1])){
                        num_byEta[0][q]->Fill(jtcorrpt[lj]);
                    }
                }
            }

            if(hlt_60==1){
                // filling hists
                num_60->Fill(jtcorrpt[lj]);
                num_60_pta->Fill(jtpt_a);
                numEta_60->Fill(jteta[lj]);
                numPhi_60->Fill(jtphi[lj]);
                if(jtcorrpt[lj]>60){
                    feEta_60->Fill(jteta[lj]);
                    fePhi_60->Fill(jtphi[lj]);
                    if(jtcorrpt[lj]>feff_jethlts_thresh[1]){
                        kin_99[1][0]->Fill(jtcorrpt[lj]);
                        kin_99[1][1]->Fill(jteta[lj]);
                        kin_99[1][2]->Fill(jtphi[lj]);
                        kin_99[1][3]->Fill(jet_asymm);
                    }
                }
                num_60_pt_eta->Fill(jtcorrpt[lj],jteta[lj]);
                // pt jet turn on curve by eta bins
                for(unsigned int q=0; q<9; q++){
                    if(q<5){
                        if(((TMath::Abs(jteta[lj])==byAbsEta[q])||(TMath::Abs(jteta[lj])>byAbsEta[q]))&&(TMath::Abs(jteta[lj])<byAbsEta[q+1])){
                            num_byAbsEta[1][q]->Fill(jtcorrpt[lj]);
                        }
                    }
                    if(((jteta[lj]==byEta[q])||(jteta[lj]>byEta[q]))&&(jteta[lj]<byEta[q+1])){
                        num_byEta[1][q]->Fill(jtcorrpt[lj]);
                    }
                }
            }
            if(hlt_80==1){
                // filling hists
                num_80->Fill(jtcorrpt[lj]);
                num_80_pta->Fill(jtpt_a);
                numEta_80->Fill(jteta[lj]);
                numPhi_80->Fill(jtphi[lj]);
                if(jtcorrpt[lj]>80){
                    feEta_80->Fill(jteta[lj]);
                    fePhi_80->Fill(jtphi[lj]);
                    if(jtcorrpt[lj]>feff_jethlts_thresh[2]){
                        kin_99[2][0]->Fill(jtcorrpt[lj]);
                        kin_99[2][1]->Fill(jteta[lj]);
                        kin_99[2][2]->Fill(jtphi[lj]);
                        kin_99[2][3]->Fill(jet_asymm);
                    }
                }
                num_80_pt_eta->Fill(jtcorrpt[lj],jteta[lj]);
                // pt jet turn on curve by eta bins
                for(unsigned int q=0; q<9; q++){
                    if(q<5){
                        if(((TMath::Abs(jteta[lj])==byAbsEta[q])||(TMath::Abs(jteta[lj])>byAbsEta[q]))&&(TMath::Abs(jteta[lj])<byAbsEta[q+1])){
                            num_byAbsEta[2][q]->Fill(jtcorrpt[lj]);
                        }
                    }
                    if(((jteta[lj]==byEta[q])||(jteta[lj]>byEta[q]))&&(jteta[lj]<byEta[q+1])){
                        num_byEta[2][q]->Fill(jtcorrpt[lj]);
                    }
                }
            }
            if(hlt_100==1){
                // filling hists
                num_100->Fill(jtcorrpt[lj]);
                num_100_pta->Fill(jtpt_a);
                numEta_100->Fill(jteta[lj]);
                numPhi_100->Fill(jtphi[lj]);
                if(jtcorrpt[lj]>100){
                    feEta_100->Fill(jteta[lj]);
                    fePhi_100->Fill(jtphi[lj]);
                    if(jtcorrpt[lj]>feff_jethlts_thresh[3]){
                        kin_99[3][0]->Fill(jtcorrpt[lj]);
                        kin_99[3][1]->Fill(jteta[lj]);
                        kin_99[3][2]->Fill(jtphi[lj]);
                        kin_99[3][3]->Fill(jet_asymm);
                    }
                }
                num_100_pt_eta->Fill(jtcorrpt[lj],jteta[lj]);
                // pt jet turn on curve by eta bins
                for(unsigned int q=0; q<9; q++){
                    if(q<5){
                        if(((TMath::Abs(jteta[lj])==byAbsEta[q])||(TMath::Abs(jteta[lj])>byAbsEta[q]))&&(TMath::Abs(jteta[lj])<byAbsEta[q+1])){
                            num_byAbsEta[3][q]->Fill(jtcorrpt[lj]);
                        }
                    }
                    if(((jteta[lj]==byEta[q])||(jteta[lj]>byEta[q]))&&(jteta[lj]<byEta[q+1])){
                        num_byEta[3][q]->Fill(jtcorrpt[lj]);
                    }
                }
            }
            if(hlt_120==1){
                // filling hists
                num_120->Fill(jtcorrpt[lj]);
                num_120_pta->Fill(jtpt_a);
                numEta_120->Fill(jteta[lj]);
                numPhi_120->Fill(jtphi[lj]);
                if(jtcorrpt[lj]>120){
                    feEta_120->Fill(jteta[lj]);
                    fePhi_120->Fill(jtphi[lj]);
                    if(jtcorrpt[lj]>feff_jethlts_thresh[4]){
                        kin_99[4][0]->Fill(jtcorrpt[lj]);
                        kin_99[4][1]->Fill(jteta[lj]);
                        kin_99[4][2]->Fill(jtphi[lj]);
                        kin_99[4][3]->Fill(jet_asymm);
                    }
                }
                num_120_pt_eta->Fill(jtcorrpt[lj],jteta[lj]);
                // pt jet turn on curve by eta bins
                for(unsigned int q=0; q<9; q++){
                    if(q<5){
                        if(((TMath::Abs(jteta[lj])==byAbsEta[q])||(TMath::Abs(jteta[lj])>byAbsEta[q]))&&(TMath::Abs(jteta[lj])<byAbsEta[q+1])){
                            num_byAbsEta[4][q]->Fill(jtcorrpt[lj]);
                        }
                    }
                    if(((jteta[lj]==byEta[q])||(jteta[lj]>byEta[q]))&&(jteta[lj]<byEta[q+1])){
                        num_byEta[4][q]->Fill(jtcorrpt[lj]);
                    }
                }
            }
        }
    }
    fi_input->Close();

    // making a new file to store all the histograms of interest in
    TFile *fi_output = new TFile(output,"recreate");
    fi_output->cd();

    // writing hists
    // denominator for all ratios
    denom->Write();
    denom_a->Write();
    denomEta->Write();
    denomPhi->Write();
    denom_pt_eta->Write();
    // 40
    num_40->Write();
    num_40_pta->Write();
    numEta_40->Write();
    feEta_40->Write();
    numPhi_40->Write();
    fePhi_40->Write();
    num_40_pt_eta->Write();
    denomEta_40->Write();
    // 60
    num_60->Write();
    num_60_pta->Write();
    numEta_60->Write();
    feEta_60->Write();
    numPhi_60->Write();
    fePhi_60->Write();
    num_60_pt_eta->Write();
    denomEta_60->Write();
    // 80
    num_80->Write();
    num_80_pta->Write();
    numEta_80->Write();
    feEta_80->Write();
    numPhi_80->Write();
    fePhi_80->Write();
    num_80_pt_eta->Write();
    denomEta_80->Write();
    // 100
    num_100->Write();
    num_100_pta->Write();
    numEta_100->Write();
    feEta_100->Write();
    numPhi_100->Write();
    fePhi_100->Write();
    num_100_pt_eta->Write();
    denomEta_100->Write();
    // 120
    num_120->Write();
    num_120_pta->Write();
    numEta_120->Write();
    feEta_120->Write();
    numPhi_120->Write();
    fePhi_120->Write();
    num_120_pt_eta->Write();
    denomEta_120->Write();

    // pt jet turn on curve by eta bins
    for(unsigned int p=0; p<5; p++){
        for(unsigned int e=0; e<9; e++){
            if(e<5){
                num_byAbsEta[p][e]->Write();
                if(p==0){
                    denom_byAbsEta[e]->Write();
                }
            }
            num_byEta[p][e]->Write();
            if(p==0){
                denom_byEta[e]->Write();
            }
        }
    }

    // saving kinematics for jet triggers above 99% efficiency
    for(unsigned int p=0; p<5; p++){
        for(unsigned int q=0; q<4; q++){
            kin_99[p][q]->Write();
        }
    }
}
