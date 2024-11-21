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
#include <iostream>
#include <iostream>

void TriggerTurnOn_lxplus_2024PbPb_v5(){

    // desired name of output file
    // TString output = "turnon_PromptRecoRawPrime0_2024PbPb_11_15_2024_run387973.root";
    // TString output = "turnon_PromptRecoRawPrime1_2024PbPb_11_15_2024_run387973.root";
    // TString output = "turnon_PromptRecoRawPrime2_2024PbPb_11_15_2024_run387973.root";
    TString output = "turnon_StreamerRawPrime0_2024PbPb_11_15_2024_run387973.root";
    // TString output = "turnon_StreamerRawPrime2_2024PbPb_11_15_2024_run387973.root";

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
    // phi
    TH1D *numPhi_40 = new TH1D("numPhi_40","numPhi_40",phih1d0[0],phih1d0[1],phih1d0[2]);

    // jets > 60 
    // pT
    TH1D *num_60 = new TH1D("num_60","num_60",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *num_60_pta = new TH1D("num_60_pta","num_60_pta",pth1d0[0],pth1d0[1],pth1d0[2]);
    // eta
    TH1D *numEta_60 = new TH1D("numEta_60","numEta_60",etah1d0[0],etah1d0[1],etah1d0[2]);
    // phi
    TH1D *numPhi_60 = new TH1D("numPhi_60","numPhi_60",phih1d0[0],phih1d0[1],phih1d0[2]);

    // jets > 80 
    // pT
    TH1D *num_80 = new TH1D("num_80","num_80",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *num_80_pta = new TH1D("num_80_pta","num_80_pta",pth1d0[0],pth1d0[1],pth1d0[2]);
    // eta
    TH1D *numEta_80 = new TH1D("numEta_80","numEta_80",etah1d0[0],etah1d0[1],etah1d0[2]);
    // phi
    TH1D *numPhi_80 = new TH1D("numPhi_80","numPhi_80",phih1d0[0],phih1d0[1],phih1d0[2]);

    // jets > 100 
    // pt jet turn on
    TH1D *num_100 = new TH1D("num_100","num_100",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *num_100_pta = new TH1D("num_100_pta","num_100_pta",pth1d0[0],pth1d0[1],pth1d0[2]);
    // eta jet turn on info
    TH1D *numEta_100 = new TH1D("numEta_100","numEta_100",etah1d0[0],etah1d0[1],etah1d0[2]);
    // phi jet turn on info
    TH1D *numPhi_100 = new TH1D("numPhi_100","numPhi_100",phih1d0[0],phih1d0[1],phih1d0[2]);

    // jets > 120 
    // pT
    TH1D *num_120 = new TH1D("num_120","num_120",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *num_120_pta = new TH1D("num_120_pta","num_120_pta",pth1d0[0],pth1d0[1],pth1d0[2]);
    // eta
    TH1D *numEta_120 = new TH1D("numEta_120","numEta_120",etah1d0[0],etah1d0[1],etah1d0[2]);
    // phi
    TH1D *numPhi_120 = new TH1D("numPhi_120","numPhi_120",phih1d0[0],phih1d0[1],phih1d0[2]);
        
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

    // ifstream myfile("/afs/cern.ch/user/n/nbarnett/public/txt_files/filename_txt_files/2024PbPb_filenames/PromptRecoRawPrime0_run878973.txt");
    // ifstream myfile("/afs/cern.ch/user/n/nbarnett/public/txt_files/filename_txt_files/2024PbPb_filenames/PromptRecoRawPrime1_run878973.txt");
    // ifstream myfile("/afs/cern.ch/user/n/nbarnett/public/txt_files/filename_txt_files/2024PbPb_filenames/PromptRecoRawPrime2_run878973.txt");
    ifstream myfile("/afs/cern.ch/user/n/nbarnett/public/txt_files/filename_txt_files/2024PbPb_filenames/StreamerRawPrime0_run878973.txt");
    // ifstream myfile("/afs/cern.ch/user/n/nbarnett/public/txt_files/filename_txt_files/2024PbPb_filenames/StreamerRawPrime2_run878973.txt");
    string filename;

    // std::ofstream outputFile("ProblemEvents_2024PbPb_PromptRecoRawPrime0_11_15_2024_run387973.txt");
    // std::ofstream outputFile("ProblemEvents_2024PbPb_PromptRecoRawPrime1_11_15_2024_run387973.txt");
    // std::ofstream outputFile("ProblemEvents_2024PbPb_PromptRecoRawPrime2_11_15_2024_run387973.txt");
    std::ofstream outputFile("ProblemEvents_2024PbPb_StreamerRawPrime0_11_15_2024_run387973.txt");
    // std::ofstream outputFile("ProblemEvents_2024PbPb_StreamerRawPrime2_11_15_2024_run387973.txt");

    // loop over the files by file names
    while(getline(myfile, filename)){

        // reading input file
        TString input = filename;
        TFile *fi_input = TFile::Open(input,"read");
        fi_input->cd();

        // cout<<"processing file "<<input<<endl;

        // making the variables to get from the ttrees
        // event variables
        Int_t nref;
        Int_t hiBin;
        Float_t vz;
        UInt_t lumi;
        ULong64_t evt;
        UInt_t run;
        Int_t ppvf;
        // jet hlt hits
        Int_t hlt_40;
        Int_t hlt_60;
        Int_t hlt_80;
        Int_t hlt_100;
        Int_t hlt_120;
        // zero bias hlt hits
        Int_t hlt_mb;
        // jet momenta
        const unsigned int mj = 10000;
        Float_t jtpt[mj];
        Float_t jteta[mj];
        Float_t jtphi[mj];
        Float_t jtchf[mj];
        Float_t jtmf[mj];

        // getting the ttrees of interest
        // triggers
        TTree *hlt_tree = (TTree*)fi_input->Get("hltanalysis/HltTree");
        // jet momentum
        TTree *jet_tree = (TTree*)fi_input->Get("akCs4PFJetAnalyzer/t");
        // TTree *jet_tree = (TTree*)fi_input->Get("akFlowPuCs4PFJetAnalyzer/t");
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
        hievt_tree->SetBranchStatus("lumi", 1);
        hievt_tree->SetBranchStatus("evt", 1);
        hievt_tree->SetBranchStatus("run", 1);
        hievt_tree->SetBranchStatus("hiBin", 1);
        pvf_tree->SetBranchStatus("pprimaryVertexFilter", 1);
        // jet hlt hits
        hlt_tree->SetBranchStatus("HLT_HIPuAK4CaloJet40Eta5p1_MinBiasHF1AND_v6", 1);
        hlt_tree->SetBranchStatus("HLT_HIPuAK4CaloJet60Eta5p1_MinBiasHF1AND_v6", 1);
        hlt_tree->SetBranchStatus("HLT_HIPuAK4CaloJet80Eta5p1_v14", 1);
        hlt_tree->SetBranchStatus("HLT_HIPuAK4CaloJet100Eta5p1_v14", 1);
        hlt_tree->SetBranchStatus("HLT_HIPuAK4CaloJet120Eta5p1_v14", 1);
        // zb hlt hits
        hlt_tree->SetBranchStatus("HLT_HIMinimumBiasHF1AND_v7", 1);
        // jet momenta
        jet_tree->SetBranchStatus("jtpt", 1);
        jet_tree->SetBranchStatus("jteta", 1);
        jet_tree->SetBranchStatus("jtphi", 1);
        jet_tree->SetBranchStatus("jtPfCHF", 1);
        jet_tree->SetBranchStatus("jtPfMUF", 1);

        // turning tree branches of interest on
        // jet hlts
        hlt_tree->SetBranchAddress("HLT_HIPuAK4CaloJet40Eta5p1_MinBiasHF1AND_v6", &hlt_40);
        hlt_tree->SetBranchAddress("HLT_HIPuAK4CaloJet60Eta5p1_MinBiasHF1AND_v6", &hlt_60);
        hlt_tree->SetBranchAddress("HLT_HIPuAK4CaloJet80Eta5p1_v14", &hlt_80);
        hlt_tree->SetBranchAddress("HLT_HIPuAK4CaloJet100Eta5p1_v14", &hlt_100);
        hlt_tree->SetBranchAddress("HLT_HIPuAK4CaloJet120Eta5p1_v14", &hlt_120);
        // zb hlt hit
        hlt_tree->SetBranchAddress("HLT_HIMinimumBiasHF1AND_v7", &hlt_mb);
        // jet momenta
        jet_tree->SetBranchAddress("jtpt", jtpt);
        jet_tree->SetBranchAddress("jteta", jteta);
        jet_tree->SetBranchAddress("jtphi", jtphi);
        jet_tree->SetBranchAddress("jtPfCHF", jtchf);
        jet_tree->SetBranchAddress("jtPfMUF", jtmf);
        // event variables
        jet_tree->SetBranchAddress("nref", &nref);
        hievt_tree->SetBranchAddress("vz", &vz);
        hievt_tree->SetBranchAddress("lumi", &lumi);
        hievt_tree->SetBranchAddress("hiBin", &hiBin);
        hievt_tree->SetBranchAddress("evt", &evt);
        hievt_tree->SetBranchAddress("run", &run);
        pvf_tree->SetBranchAddress("pprimaryVertexFilter", &ppvf);

        // event loop
        for(unsigned int i=0; i<jet_tree->GetEntries(); i++){
        // for(unsigned int i=0; i<1000; i++){

            njets+=1;
            
            // event filtering
            pvf_tree->GetEntry(i);
            // if(ppvf!=1){continue;}
            hievt_tree->GetEntry(i);
            // if(TMath::Abs(vz)>15){continue;}
            // if(lumi<35){continue;}
            // if(lumi>63){continue;}

            njets_+=1;

            jet_tree->GetEntry(i);
            hlt_tree->GetEntry(i);

            // if(hlt_mb==1){
            //     cout<<hlt_mb<<" is the value of the minbias hlt"<<endl;
            //     cout<< i <<" is the entry number"<<endl;
            //     cout<< nref <<" is nref"<<endl;
            //     cout<< jtpt[0] <<" is leading jet pt"<<endl;
            // }

            // excluding events without jets or outside of the hlt eta range
            // if((nref==0)||(TMath::Abs(jteta[0])>5.1)){continue;}

            // excluding events where the leading jet has a high muon or charged hadron fraction
            // if(jtchf[0]>0.99){continue;}
            // if(jtmf[0]>0.8){continue;}

            Int_t lj = 0;
            Int_t slj = 1;

            if((evt==73778532)&&(lumi==114)){
                cout<<filename<<" and entry is "<<i<<endl;
            }

            // getting leading jet 
            for(Int_t j = 0; j<nref; j++){
                if(jtpt[j]>jtpt[lj]){
                    lj=j;
                    if(lj==1){slj=0;}
                }
            }

            // getting subleading jet index
            for(Int_t j = 0; j<nref; j++){
                if((jtpt[j]>jtpt[slj])&&(j!=lj)){
                    slj=j;
                }
            }

            // threshholds for the jet hlts when they reach 99% efficiency
            // // 2024 ppRef
            // Double_t feff_jethlts_thresh[5] = {58.5197,9999,99.2784,115.987,141.701};
            // 2024 PbPb
            Double_t feff_jethlts_thresh[5] = {84.8729,129.738,138.933,148.207,161.07};

            // filling hists with leading jet info 
            if((jtpt[lj]>0)&&(hlt_mb==1)){
                Double_t jet_asymm = (jtpt[lj] - jtpt[slj])/(jtpt[lj] + jtpt[slj]);
                Double_t jtpt_a = (jtpt[lj] + jtpt[slj])/2;
                denom->Fill(jtpt[lj]);
                denom_a->Fill(jtpt_a);
                denomEta->Fill(jteta[lj]);
                denomPhi->Fill(jtphi[lj]);
                denom_pt_eta->Fill(jtpt[lj],jteta[lj]);

                for(unsigned int q=0; q<9; q++){
                    if(q<5){
                        if(((TMath::Abs(jteta[lj])==byAbsEta[q])||(TMath::Abs(jteta[lj])>byAbsEta[q]))&&(TMath::Abs(jteta[lj])<byAbsEta[q+1])){
                            denom_byAbsEta[q]->Fill(jtpt[lj]);
                        }
                    }
                    if(((jteta[lj]==byEta[q])||(jteta[lj]>byEta[q]))&&(jteta[lj]<byEta[q+1])){
                        denom_byEta[q]->Fill(jtpt[lj]);
                    }
                }

                njets_0+=1;
                if(jtpt[lj]>20){njets_20+=1;}
                if(jtpt[lj]>40){
                    njets_40+=1;
                    denomEta_40->Fill(jteta[lj]);
                }
                if(jtpt[lj]>60){
                    njets_60+=1;
                    denomEta_60->Fill(jteta[lj]);
                }
                if(jtpt[lj]>80){
                    njets_80+=1;
                    denomEta_80->Fill(jteta[lj]);
                }
                if(jtpt[lj]>100){
                    njets_100+=1;
                    denomEta_100->Fill(jteta[lj]);
                }
                if(jtpt[lj]>120){
                    njets_120+=1;
                    denomEta_120->Fill(jteta[lj]);
                }
                int flags[5] = {0,0,0,0,0};
                if((jtpt[lj]>feff_jethlts_thresh[0]+20)&&(hlt_40==0)){flags[0]=1;}
                if((jtpt[lj]>feff_jethlts_thresh[1]+20)&&(hlt_60==0)){flags[1]=1;}
                if((jtpt[lj]>feff_jethlts_thresh[2]+20)&&(hlt_80==0)){flags[2]=1;}
                if((jtpt[lj]>feff_jethlts_thresh[3]+20)&&(hlt_100==0)){flags[3]=1;}
                if((jtpt[lj]>feff_jethlts_thresh[4]+20)&&(hlt_120==0)){flags[4]=1;}
                if((flags[0]==1)||(flags[1]==1)||(flags[2]==1)||(flags[3]==1)||(flags[4]==1)){
                    // cout<<"entry "<<i<<endl;
                    // cout<<"leading jet pt "<<jtpt[lj]<<endl;
                    outputFile<<"\n" << filename<<"\n";
                    outputFile << "run is "<<run<<", lumi is "<<lumi<<", event is "<<evt<<", entry is "<<i<<"\n";
                    outputFile<<"hiBin is "<<hiBin << ", leading jet pt is "<<jtpt[lj]<<", charged hadron fraction for that jet is "<<jtchf[lj]<<"\n";
                    if(hlt_40==0){
                        // cout<<"HLT_HIPuAK4CaloJet40Eta5p1_MinBiasHF1AND_v6 didn't fire"<<endl;
                        outputFile<<"HLT_HIPuAK4CaloJet40Eta5p1_MinBiasHF1AND_v6 didn't fire\n";
                    }
                    if(hlt_60==0){
                        // cout<<"HLT_HIPuAK4CaloJet60Eta5p1_MinBiasHF1AND_v6 didn't fire"<<endl;
                        outputFile<<"HLT_HIPuAK4CaloJet60Eta5p1_MinBiasHF1AND_v6 didn't fire\n";
                    }
                    if(hlt_80==0){
                        // cout<<"HLT_HIPuAK4CaloJet80Eta5p1_v14 didn't fire"<<endl;
                        outputFile<<"HLT_HIPuAK4CaloJet80Eta5p1_v14 didn't fire\n";
                    }
                    if(hlt_100==0){
                        // cout<<"HLT_HIPuAK4CaloJet100Eta5p1_v14 didn't fire"<<endl;
                        outputFile<<"HLT_HIPuAK4CaloJet100Eta5p1_v14 didn't fire\n";
                    }
                    if(hlt_120==0){
                        // cout<<"HLT_HIPuAK4CaloJet120Eta5p1_v14 didn't fire"<<endl;
                        outputFile<<"HLT_HIPuAK4CaloJet120Eta5p1_v14 didn't fire\n";
                    }
                }
                // if the triggers are on
                if(hlt_40==1){
                    // filling hists
                    num_40->Fill(jtpt[lj]);
                    num_40_pta->Fill(jtpt_a);
                    numEta_40->Fill(jteta[lj]);
                    numPhi_40->Fill(jtphi[lj]);
                    if(jtpt[lj]>40){
                        if(jtpt[lj]>feff_jethlts_thresh[0]){
                            kin_99[0][0]->Fill(jtpt[lj]);
                            kin_99[0][1]->Fill(jteta[lj]);
                            kin_99[0][2]->Fill(jtphi[lj]);
                            kin_99[0][3]->Fill(jet_asymm);
                        }
                    }
                    num_40_pt_eta->Fill(jtpt[lj],jteta[lj]);
                    // pt jet turn on curve by eta bins
                    for(unsigned int q=0; q<9; q++){
                        if(q<5){
                            if(((TMath::Abs(jteta[lj])==byAbsEta[q])||(TMath::Abs(jteta[lj])>byAbsEta[q]))&&(TMath::Abs(jteta[lj])<byAbsEta[q+1])){
                                num_byAbsEta[0][q]->Fill(jtpt[lj]);
                            }
                        }
                        if(((jteta[lj]==byEta[q])||(jteta[lj]>byEta[q]))&&(jteta[lj]<byEta[q+1])){
                            num_byEta[0][q]->Fill(jtpt[lj]);
                        }
                    }
                }

                if(hlt_60==1){
                    // filling hists
                    num_60->Fill(jtpt[lj]);
                    num_60_pta->Fill(jtpt_a);
                    numEta_60->Fill(jteta[lj]);
                    numPhi_60->Fill(jtphi[lj]);
                    if(jtpt[lj]>60){
                        if(jtpt[lj]>feff_jethlts_thresh[1]){
                            kin_99[1][0]->Fill(jtpt[lj]);
                            kin_99[1][1]->Fill(jteta[lj]);
                            kin_99[1][2]->Fill(jtphi[lj]);
                            kin_99[1][3]->Fill(jet_asymm);
                        }
                    }
                    num_60_pt_eta->Fill(jtpt[lj],jteta[lj]);
                    // pt jet turn on curve by eta bins
                    for(unsigned int q=0; q<9; q++){
                        if(q<5){
                            if(((TMath::Abs(jteta[lj])==byAbsEta[q])||(TMath::Abs(jteta[lj])>byAbsEta[q]))&&(TMath::Abs(jteta[lj])<byAbsEta[q+1])){
                                num_byAbsEta[1][q]->Fill(jtpt[lj]);
                            }
                        }
                        if(((jteta[lj]==byEta[q])||(jteta[lj]>byEta[q]))&&(jteta[lj]<byEta[q+1])){
                            num_byEta[1][q]->Fill(jtpt[lj]);
                        }
                    }
                }
                if(hlt_80==1){
                    // filling hists
                    num_80->Fill(jtpt[lj]);
                    num_80_pta->Fill(jtpt_a);
                    numEta_80->Fill(jteta[lj]);
                    numPhi_80->Fill(jtphi[lj]);
                    if(jtpt[lj]>80){
                        if(jtpt[lj]>feff_jethlts_thresh[2]){
                            kin_99[2][0]->Fill(jtpt[lj]);
                            kin_99[2][1]->Fill(jteta[lj]);
                            kin_99[2][2]->Fill(jtphi[lj]);
                            kin_99[2][3]->Fill(jet_asymm);
                        }
                    }
                    num_80_pt_eta->Fill(jtpt[lj],jteta[lj]);
                    // pt jet turn on curve by eta bins
                    for(unsigned int q=0; q<9; q++){
                        if(q<5){
                            if(((TMath::Abs(jteta[lj])==byAbsEta[q])||(TMath::Abs(jteta[lj])>byAbsEta[q]))&&(TMath::Abs(jteta[lj])<byAbsEta[q+1])){
                                num_byAbsEta[2][q]->Fill(jtpt[lj]);
                            }
                        }
                        if(((jteta[lj]==byEta[q])||(jteta[lj]>byEta[q]))&&(jteta[lj]<byEta[q+1])){
                            num_byEta[2][q]->Fill(jtpt[lj]);
                        }
                    }
                }
                if(hlt_100==1){
                    // filling hists
                    num_100->Fill(jtpt[lj]);
                    num_100_pta->Fill(jtpt_a);
                    numEta_100->Fill(jteta[lj]);
                    numPhi_100->Fill(jtphi[lj]);
                    if(jtpt[lj]>100){
                        if(jtpt[lj]>feff_jethlts_thresh[3]){
                            kin_99[3][0]->Fill(jtpt[lj]);
                            kin_99[3][1]->Fill(jteta[lj]);
                            kin_99[3][2]->Fill(jtphi[lj]);
                            kin_99[3][3]->Fill(jet_asymm);
                        }
                    }
                    num_100_pt_eta->Fill(jtpt[lj],jteta[lj]);
                    // pt jet turn on curve by eta bins
                    for(unsigned int q=0; q<9; q++){
                        if(q<5){
                            if(((TMath::Abs(jteta[lj])==byAbsEta[q])||(TMath::Abs(jteta[lj])>byAbsEta[q]))&&(TMath::Abs(jteta[lj])<byAbsEta[q+1])){
                                num_byAbsEta[3][q]->Fill(jtpt[lj]);
                            }
                        }
                        if(((jteta[lj]==byEta[q])||(jteta[lj]>byEta[q]))&&(jteta[lj]<byEta[q+1])){
                            num_byEta[3][q]->Fill(jtpt[lj]);
                        }
                    }
                }
                if(hlt_120==1){
                    // filling hists
                    num_120->Fill(jtpt[lj]);
                    num_120_pta->Fill(jtpt_a);
                    numEta_120->Fill(jteta[lj]);
                    numPhi_120->Fill(jtphi[lj]);
                    if(jtpt[lj]>120){
                        if(jtpt[lj]>feff_jethlts_thresh[4]){
                            kin_99[4][0]->Fill(jtpt[lj]);
                            kin_99[4][1]->Fill(jteta[lj]);
                            kin_99[4][2]->Fill(jtphi[lj]);
                            kin_99[4][3]->Fill(jet_asymm);
                        }
                    }
                    num_120_pt_eta->Fill(jtpt[lj],jteta[lj]);
                    // pt jet turn on curve by eta bins
                    for(unsigned int q=0; q<9; q++){
                        if(q<5){
                            if(((TMath::Abs(jteta[lj])==byAbsEta[q])||(TMath::Abs(jteta[lj])>byAbsEta[q]))&&(TMath::Abs(jteta[lj])<byAbsEta[q+1])){
                                num_byAbsEta[4][q]->Fill(jtpt[lj]);
                            }
                        }
                        if(((jteta[lj]==byEta[q])||(jteta[lj]>byEta[q]))&&(jteta[lj]<byEta[q+1])){
                            num_byEta[4][q]->Fill(jtpt[lj]);
                        }
                    }
                }
            }
        }
        fi_input->Close();
    }

    // closing the txt file made to store problem events
    outputFile.close();

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
    numPhi_40->Write();
    num_40_pt_eta->Write();
    denomEta_40->Write();
    // 60
    num_60->Write();
    num_60_pta->Write();
    numEta_60->Write();
    numPhi_60->Write();
    num_60_pt_eta->Write();
    denomEta_60->Write();
    // 80
    num_80->Write();
    num_80_pta->Write();
    numEta_80->Write();
    numPhi_80->Write();
    num_80_pt_eta->Write();
    denomEta_80->Write();
    // 100
    num_100->Write();
    num_100_pta->Write();
    numEta_100->Write();
    numPhi_100->Write();
    num_100_pt_eta->Write();
    denomEta_100->Write();
    // 120
    num_120->Write();
    num_120_pta->Write();
    numEta_120->Write();
    numPhi_120->Write();
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

    // printing jet numbers
    // cout<<njets<<" is the number of leading jets processed"<<endl;
    // cout<<njets_0<<" is the number of leading jets passed the zb"<<endl;
    // cout<<njets_20<<" is the number of leading jets with pt above 20 GeV"<<endl;
    // cout<<njets_40<<" is the number of leading jets with pt above 40 GeV"<<endl;
    // cout<<njets_60<<" is the number of leading jets with pt above 60 GeV"<<endl;
    // cout<<njets_80<<" is the number of leading jets with pt above 80 GeV"<<endl;
    // cout<<njets_100<<" is the number of leading jets with pt above 100 GeV"<<endl;
    // cout<<njets_120<<" is the number of leading jets with pt above 120 GeV"<<endl;
}
