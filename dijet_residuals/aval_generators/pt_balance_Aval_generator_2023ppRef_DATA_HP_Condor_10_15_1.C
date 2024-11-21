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

// PLOTTING FUNCTIONS

void save_h1d(TH1D *h, TString hname){
    h->SetName(hname);
    h->Write();
}

void save_h2d(TH2D *h, TString hname){
    h->SetName(hname);
    h->Write();
}

// the script all runs in this function
void pt_balance_Aval_generator_2023ppRef_DATA_HP_Condor_10_15_1(TString input, TString output)
{   
    // taking into account the appropriate errors
    TH1::SetDefaultSumw2(); 
    
    // INITIALIZING HISTOGRAMS

    // creating some binning parameters
    // _h1dN is the Nth set of binning parameters for 1d hists of _
    const double vzh1d0[3] = {40,-20,20};
    const double pth1d0[3] = {100,15,500};
    const double phih1d0[3] = {100,-1*TMath::Pi(),TMath::Pi()};
    const double etah1d0[3] = {50,-5.2,5.2};
    const double etah1d1[3] = {25,-1.7,1.7};
    const double ah1d0[3] = {100,-1,1};

    // a values
    // number of values a can assume
    const Int_t anum = 100;
    const Int_t anum1 = 101;
    // the a values 
    double as[anum1] = {0};
    for(unsigned int a=0; a<anum1; a++){
        double aa = a;
        double asval = -1.0 + aa*0.02;
        as[a] = asval;
    }

    // ptslices
    // number of pt slices
    const Int_t ptslicenum = 4;
    const Int_t ptslicenum1 = 5;
    // the low and high pt values for each pt slice
    double pts[ptslicenum1] = {55,80,120,170,1000};
    double ptlow[ptslicenum] = {55,80,120,170};
    double pthigh[ptslicenum] = {80,120,170,1000};
    
    // eta slices
    // number of eta slices
    const Int_t etaslicenum = 36;
    const Int_t etaslicenum1 = 37;
    // the low and high eta values for each eta slice
    // double oldetas0[?] = {-5.2,-3.9,-2.6,-1.3,0,1.3,2.6,3.9,5.2};
    // double oldetas1[?] = {-5.2,-4.5,-3.9,-3.5,-3.2,-2.9,-2.6,-2.2,-1.9,-1.6,-1.3,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.3,1.6,1.9,2.2,2.6,2.9,3.2,3.5,3.9,4.5,5.2};
    // double oldetas2[20 +1] = {-5.2, -4.2, -3.6, -3.0, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 3.6, 4.2, 5.2};
    // double oldetas3[14 +1] = {-5.2, -4.2, -3.3, -2.4, -1.8, -1.2, -0.6, 0.0, 0.6, 1.2, 1.8, 2.4, 3.3, 4.2, 5.2};
    // double etas[82 +1] = {-5.191,-4.889,-4.716,-4.538,-4.363,-4.191,-4.013,-3.839,-3.664,-3.489,-3.314,-3.139,-2.964,-2.853,-2.65,-2.5,-2.322,-2.172,-2.043,-1.93,-1.83,-1.74,-1.653,-1.566,-1.479,-1.392,-1.305,-1.218,-1.131,-1.044,-0.957,-0.879,-0.783,-0.696,-0.609,-0.522,-0.435,-0.348,-0.261,-0.174,-0.087,0,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696,0.783,0.879,0.957,1.044,1.131,1.218,1.305,1.392,1.479,1.566,1.653,1.74,1.83,1.93,2.043,2.172,2.322,2.5,2.65,2.853,2.964,3.139,3.314,3.489,3.664,3.839,4.013,4.191,4.363,4.538,4.716,4.889,5.191};
    // double etas[36 +1] = {-5.191,-3.839,-3.489,-3.139,-2.964,-2.853,-2.65,-2.5,-2.322,-2.172,-1.93,-1.653,-1.479,-1.305,-1.044,-0.783,-0.522,-0.261,0,0.261,0.522,0.783,1.044,1.305,1.479,1.653,1.93,2.172,2.322,2.5,2.65,2.853,2.964,3.139,3.489,3.839,5.191};
    // new binnning
    double etas[etaslicenum1] = {-5.191,-3.839,-3.489,-3.139,-2.964,-2.853,-2.65,-2.5,-2.322,-2.172,-1.93,-1.653,-1.479,-1.305,-1.044,-0.783,-0.522,-0.261,0,0.261,0.522,0.783,1.044,1.305,1.479,1.653,1.93,2.172,2.322,2.5,2.65,2.853,2.964,3.139,3.489,3.839,5.191};
    double etahigh[etaslicenum] = {-3.839,-3.489,-3.139,-2.964,-2.853,-2.65,-2.5,-2.322,-2.172,-1.93,-1.653,-1.479,-1.305,-1.044,-0.783,-0.522,-0.261,0,0.261,0.522,0.783,1.044,1.305,1.479,1.653,1.93,2.172,2.322,2.5,2.65,2.853,2.964,3.139,3.489,3.839,5.191};
    double etalow[etaslicenum] = {-5.191,-3.839,-3.489,-3.139,-2.964,-2.853,-2.65,-2.5,-2.322,-2.172,-1.93,-1.653,-1.479,-1.305,-1.044,-0.783,-0.522,-0.261,0,0.261,0.522,0.783,1.044,1.305,1.479,1.653,1.93,2.172,2.322,2.5,2.65,2.853,2.964,3.139,3.489,3.839};

    // alpha bins
    // number of alpha bins
    const Int_t alphabinnum = 9;
    double alphas[alphabinnum] = {0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45};

    // for generating random numbers
    TRandom2 *rand = new TRandom2(1);

    // Making some hists for pt balance studies
    // the actual A and R values
    // <A> as well as R vs pt and eta for eta and pt slices respectively
    // eta and pt bin slices for the pt slice hists of A vs eta_probe and pt_avg
    // pt slices
    TH2D *ptslicesA[alphabinnum][ptslicenum];
    TH1D *etaslices_of_ptslicesA[alphabinnum][ptslicenum][etaslicenum];
    // eta slices
    TH2D *etaslicesA[alphabinnum][etaslicenum];
    TH1D *ptslices_of_etaslicesA[alphabinnum][etaslicenum][ptslicenum];
    
    // further initializing the histograms
    // looping over alpha bins
    for(unsigned int a=0; a<alphabinnum; a++){
        // looping over pt slices
        for(unsigned int p=0; p<ptslicenum; p++){
            // A vs eta for each pt slice with title having the high and low pt for the slice
            TString chtitle0 = Form("A_ptslice_%.0f_%.0f__alpha_%.0f",ptlow[p],pthigh[p],alphas[a]*100);
            // ptslicesA[p] = new TH2D(chtitle0,chtitle0,etah1d0[0],etah1d0[1],etah1d0[2],anum,as);
            ptslicesA[a][p] = new TH2D(chtitle0,chtitle0,etaslicenum,etas,anum,as);
            for(unsigned int q=0; q<etaslicenum; q++){
                // avoiding looping through the etaslicenum separately by only doing it on the first p value
                if(p==0){
                    // A vs pt for each eta slice with title having the high and low eta for the slice
                    if(etalow[q]<0){
                        TString ahtitle0 = Form("A_etaslice_%.0f_%.0f_n__alpha_%.0f",TMath::Abs(etalow[q]*10),TMath::Abs(etahigh[q]*10),alphas[a]*100);
                        etaslicesA[a][q] = new TH2D(ahtitle0,ahtitle0,ptslicenum,pts,anum,as);
                    }
                    if(etalow[q]>0||etalow[q]==0){
                        TString ahtitle0 = Form("A_etaslice_%.0f_%.0f__alpha_%.0f",etalow[q]*10,etahigh[q]*10,alphas[a]*100);
                        etaslicesA[a][q] = new TH2D(ahtitle0,ahtitle0,ptslicenum,pts,anum,as);
                    }
                }
                // intializing hists that are the bins of the slices
                if(etalow[q]<0){
                    TString dhtitle0 = Form("A_ptslice_%.0f_%.0f__etabin_%.0f_%.0f_n__alpha_%.0f",ptlow[p],pthigh[p],TMath::Abs(etalow[q]*10),TMath::Abs(etahigh[q]*10),alphas[a]*100);
                    TString dhtitle1 = Form("A_etaslice_%.0f_%.0f_n__ptbin_%.0f_%.0f__alpha_%.0f",TMath::Abs(etalow[q]*10),TMath::Abs(etahigh[q]*10),ptlow[p],pthigh[p],alphas[a]*100);
                    etaslices_of_ptslicesA[a][p][q] = new TH1D(dhtitle0,dhtitle0,anum,as);
                    ptslices_of_etaslicesA[a][q][p] = new TH1D(dhtitle1,dhtitle1,anum,as);
                }
                if(etalow[q]>0||etalow[q]==0){
                    TString dhtitle0 = Form("A_ptslice_%.0f_%.0f__etabin_%.0f_%.0f__alpha_%.0f",ptlow[p],pthigh[p],etalow[q]*10,etahigh[q]*10,alphas[a]*100);
                    TString dhtitle1 = Form("A_etaslice_%.0f_%.0f__ptbin_%.0f_%.0f__alpha_%.0f",etalow[q]*10,etahigh[q]*10,ptlow[p],pthigh[p],alphas[a]*100);
                    etaslices_of_ptslicesA[a][p][q] = new TH1D(dhtitle0,dhtitle0,anum,as);
                    ptslices_of_etaslicesA[a][q][p] = new TH1D(dhtitle1,dhtitle1,anum,as);
                }
            }
        }
    }

    // initializing histograms for general parameters
    TH1D *htrig = new TH1D("htrig","htrig",2,0,2);
    TH1D *hvz = new TH1D("vz","vz",vzh1d0[0],vzh1d0[1],vzh1d0[2]);
    TH1D *hvz_ = new TH1D("hvz_","hvz_",vzh1d0[0],vzh1d0[1],vzh1d0[2]);
    TH1D *hjtcorrpt = new TH1D("hjtcorrpt","hjtcorrpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hrawpt = new TH1D("hrawpt","hrawpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hjteta = new TH1D("hjteta","hjteta",etah1d1[0],etah1d1[1],etah1d1[2]);
    TH1D *hjtphi = new TH1D("hjtphi","hjtphi",phih1d0[0],phih1d0[1],phih1d0[2]);

    // pt balance probe and tag momenta hists
    TH1D *htagjtpt = new TH1D("htagjtpt","htagjtpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *hprobejtpt = new TH1D("hprobejtpt","hprobejtpt",pth1d0[0],pth1d0[1],pth1d0[2]);
    TH1D *htagjteta = new TH1D("htagjteta","htagjteta",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *hprobejteta = new TH1D("hprobejteta","hprobejteta",etah1d0[0],etah1d0[1],etah1d0[2]);
    TH1D *htagjtphi = new TH1D("htagjtphi","htagjtphi",phih1d0[0],phih1d0[1],phih1d0[2]);
    TH1D *hprobejtphi = new TH1D("hprobejtphi","hprobejtphi",phih1d0[0],phih1d0[1],phih1d0[2]);

    // INITIALIZING PARAMETERS

    // declaring variables
    // vertex position
    Float_t vz;
    // whatever ptcut is decided
    Float_t ptcut = 30;
    // a big number to make my arrays such that they aren't too small
    // need to have more entries in the arrays than number of jets in the event with most jets
    const Int_t nm = 200000;
    // uncorrected jet pt
    Float_t rawpt[nm];
    // corrected jet pt
    Float_t jtcorrpt[nm];
    // jet phi and pseudorapadity
    Float_t jtphi[nm];
    Float_t jteta[nm];
    // number of jets in event
    Int_t nref;
    // primary vertex filter declaration
    int pPVF;
    // jet trigger of interest declaration
    int jet_trigger;

    // reading and iterating through a list of files instead of a single file
    TFile *fi = TFile::Open(input,"read");

    // getting the TTrees from the input file
    // event info of interest
    TTree *t0 = (TTree*)fi->Get("ak4PFJetAnalyzer/t");

    // vertex position
    TTree *t1 = (TTree*)fi->Get("hiEvtAnalyzer/HiTree");

    // primary vertex filter
    TTree *t2 = (TTree*)fi->Get("skimanalysis/HltTree");

    // jet trigger
    TTree *t3 = (TTree*)fi->Get("hltanalysis/HltTree");
    
    // turning off all branches and then turning on only the ones I want
    t0->SetBranchStatus("*",0);
    t0->SetBranchStatus("jteta",1);
    t0->SetBranchStatus("jtphi",1);
    t0->SetBranchStatus("rawpt",1);
    t0->SetBranchStatus("nref",1);
    
    // doing the same for the other trees I want
    // t1
    t1->SetBranchStatus("*",0);
    t1->SetBranchStatus("vz",1);
    // t2
    t2->SetBranchStatus("*",0);
    t2->SetBranchStatus("pprimaryVertexFilter",1);
    // t3
    t3->SetBranchStatus("*",0);
    t3->SetBranchStatus("HLT_AK4PFJet60_v1",1);
    // t3->SetBranchStatus("L1_SingleJet24",1);   

    // general parameters
    t0->SetBranchAddress("jteta",jteta);
    t0->SetBranchAddress("jtphi",jtphi);
    t0->SetBranchAddress("rawpt",rawpt);
    t0->SetBranchAddress("nref",&nref);
    t1->SetBranchAddress("vz",&vz);
    t2->SetBranchAddress("pprimaryVertexFilter",&pPVF);
    t3->SetBranchAddress("HLT_AK4PFJet60_v1",&jet_trigger);

    // for loop going over events in the trees
    for(unsigned int i=0; i<t0->GetEntries(); i++){
    // for(unsigned int i=0; i<10000; i++){

        // EVENT CUT
        // event needs to have |vz|<15
        t1->GetEntry(i);
        hvz_->Fill(vz);
        if(TMath::Abs(vz)<15){
            // event needs to pass primary vertex filter
            t2->GetEntry(i);
            // and the trigger
            t3->GetEntry(i);
	        htrig->Fill(jet_trigger);
            if((pPVF==1)&&(jet_trigger==1)){
                t0->GetEntry(i);
                
                // not looking at events without a dijet
                if(nref<2){continue;}

                // filling the vertex position hist
                hvz->Fill(vz);

                // looping through all jets in each event
                for(unsigned int j=0; j<nref; j++){
                
                    hrawpt->Fill(rawpt[j]);

                    // getting the corrected jet pt
                    vector<string> Files;
                    Files.push_back("/afs/cern.ch/user/n/nbarnett/public/txt_files/L2L3_ppReco_2023ppRef/ParallelMC_L2Relative_AK4PF_pp_Reco_v0_12-21-2023.txt");
                    JetCorrector JEC(Files);
                    JEC.SetJetPT(rawpt[j]);
                    JEC.SetJetEta(jteta[j]);
                    JEC.SetJetPhi(jtphi[j]);  
                    Float_t jet_pt_corr = JEC.GetCorrectedPT();

                    // saving the corrected jet pt
                    jtcorrpt[j] = jet_pt_corr;

                    // only look at pt balance studies if jtcorrpt > ptcut
                    if(jtcorrpt[j]>ptcut){
                        // filling histograms that have variables with more than one value per event
                        hjteta->Fill(jteta[j]);
                        hjtphi->Fill(jtphi[j]);
                        hjtcorrpt->Fill(jtcorrpt[j]);
                    }
                }
                
                // PT BALANCE
                // making the iterator values for tag, probe, and third leading jet (if there is one)
                // tag jet must be leading, and probe jet must be subleading
                // then adjusting the iterator values depending on pt order of jets in the event
                int leaditer = 0;
                int subleaditer = 1;
                int tagiter = 0;
                int probeiter = 1;
                // finding the leading jet, which is the possible tag jet
                // looping through the jets in the event
                for(unsigned int j=0; j<nref; j++){
                    // only when the jet pt is highest will it replace the current iterator
                    // in the case this is never true the original leading jet would still be the leading jet and iter would be 0
                    if(jtcorrpt[j]>jtcorrpt[leaditer]){
                        leaditer=j;
                        // if now the leading is the original subleading jet
                        // then the subleading jet is changed to another iter before being found
                        // this will be true iff tagiter = 1 or the leading jet is the original subleading jet index
                        if(subleaditer==leaditer){
                            subleaditer=0;
                        }
                    }
                }
                // finding the subleading jet, which is the possible probe jet
                // looping through the jets in the event
                for(unsigned int j=0; j<nref; j++){
                    // only when the jet pt is larger than current subleading jet and smaller than leading jet
                    // in the case this is never true the original subleading, or possibly leading, jet would be the subleading jet and iter would be 1, or possibly 0
                    if((jtcorrpt[j]<jtcorrpt[leaditer])&&(jtcorrpt[j]>jtcorrpt[subleaditer])){
                        subleaditer=j;
                    }
                }
                // pt balance study for the case nref < 3
                if(nref==2){
                    // defining tag and probe iters based on leading and subleading jet iters
                    // tag and probe must be either leading or subleading jet
                    // if the leading jet is in the barrel, tag it
                    if(TMath::Abs(jteta[leaditer])<1.3){
                        tagiter = leaditer;
                        probeiter = subleaditer;
                    }
                    // if the subleading jet is in the barrel, tag it
                    if(TMath::Abs(jteta[subleaditer])<1.3){
                        tagiter = subleaditer;
                        probeiter = leaditer;
                    }
                    // in the case both jets are in the barrel we make a random number 
                    // if the random number is even or odd the tag jet is the leading or subleading jet respectively
                    if((TMath::Abs(jteta[leaditer])<1.3)&&(TMath::Abs(jteta[subleaditer])<1.3)){
                        int checkval1 = rand->Integer(100);
                        if((checkval1%2==0)&&(nref<3)){
                            tagiter = leaditer;
                            probeiter = subleaditer;
                        }
                        if((checkval1%2!=0)&&(nref<3)){
                            tagiter = subleaditer;
                            probeiter = leaditer;
                        }
                    }
                }
                // finding the third leading jet, iff there are at least three jets
                if(nref>2){
                    int thirditer = 2;
                    // start assuming the third leading jet is the third leading jet still, unless either the leading jet or subleading jet is already
                    // only looking at the first three jets
                    for(unsigned int q=0; q<3; q++){
                        // one of the first three jets isn't the leading or subleading jet and we initialize the third jet to be that one
                        if((q!=subleaditer)&&(q!=leaditer)){
                            thirditer = q;
                        }
                    }
                    // looping through the jets in the event
                    for(unsigned int j=0; j<nref; j++){
                        // determining if another jet between index 3 and nref-1 exists with higher pt than the current third jet, but only if it is less pt than and isn't the probe or tag jet
                        if((jtcorrpt[j]>jtcorrpt[thirditer])&&(jtcorrpt[j]<jtcorrpt[subleaditer])&&(jtcorrpt[j]<jtcorrpt[leaditer])&&(j!=subleaditer)&&(j!=leaditer)){
                            thirditer = j;
                        }
                    }
                    // still working if there is a third jet
                    // doing the whole pt balance study in the case there is a third jet
                    if((TMath::ACos(TMath::Cos(jtphi[leaditer]-jtphi[subleaditer]))>2.7)&&(jtcorrpt[subleaditer]>ptcut)&&((TMath::Abs(jteta[subleaditer])<1.3)||(TMath::Abs(jteta[leaditer])<1.3))){
                        if(TMath::Abs(jteta[leaditer])<1.3){
                            tagiter = leaditer;
                            probeiter = subleaditer;
                        }
                        if(TMath::Abs(jteta[subleaditer])<1.3){
                            tagiter = subleaditer;
                            probeiter = leaditer;
                        }
                        if((TMath::Abs(jteta[leaditer])<1.3)&&(TMath::Abs(jteta[subleaditer])<1.3)){
                            int checkval = rand->Integer(100);
                            if(checkval%2==0){
                                tagiter = leaditer;
                                probeiter = subleaditer;
                            }
                            if(checkval%2!=0){
                                tagiter = subleaditer;
                                probeiter = leaditer;
                            }
                        }

                        // making and saving the A value
                        double Aval = (jtcorrpt[probeiter]-jtcorrpt[tagiter])/(jtcorrpt[probeiter]+jtcorrpt[tagiter]);
                        // removing outlier A values outside of -0.7 to 0.7
                        if(TMath::Abs(Aval)>0.7){continue;}

                        // saving probe and tag momenta
                        htagjtpt->Fill(jtcorrpt[tagiter]);
                        hprobejtpt->Fill(jtcorrpt[probeiter]);
                        htagjteta->Fill(jteta[tagiter]);
                        hprobejteta->Fill(jteta[probeiter]);
                        htagjtphi->Fill(jtphi[tagiter]);
                        hprobejtphi->Fill(jtphi[probeiter]);
                        // alpha bins
                        for(unsigned int a=0; a<alphabinnum; a++){
                            if(((alphas[a]/2)*(jtcorrpt[probeiter]+jtcorrpt[tagiter]))>jtcorrpt[thirditer]){
                                // pt slices A value filling
                                // each p is a different slice of pt
                                for(unsigned int p=0; p<ptslicenum; p++){
                                    // average momentum between the probe and tag jet 
                                    // these are sliced originally to get A vs eta for different pt slices
                                    double ptavg = (jtcorrpt[probeiter]+jtcorrpt[tagiter])/2;
                                    // if the pt avg is within a certain slice range then add it to the pt slice hists
                                    if((ptavg>ptlow[p])&&(ptavg<pthigh[p])){
                                        // ptslicesA are A vs eta_probe hists for different pt ranges or slices
                                        ptslicesA[a][p]->Fill(jteta[probeiter],Aval);
                                    }
                                }
                                // eta slices A value filling
                                // each q is a different slice of eta
                                for(unsigned int q=0; q<etaslicenum; q++){
                                    // ptavg is the x axis in one desired type of plot
                                    double ptavg = (jtcorrpt[probeiter]+jtcorrpt[tagiter])/2;
                                    // if the eta is within a certain slice range then add it to the eta slice hists
                                    if((jteta[probeiter]>etalow[q])&&(jteta[probeiter]<etahigh[q])){
                                        // etaslicesA are A vs pt_avg hists for different eta ranges
                                        etaslicesA[a][q]->Fill(ptavg,Aval);
                                    }
                                }
                            }
                        }
                    }              
                }
                // finding the A values iff the leading or subleading jet has eta < 1.3 and subleading jet passes the pt cut
                if((TMath::Abs(jteta[tagiter])<1.3)&&(jtcorrpt[subleaditer]>ptcut)&&(nref==2)&&(TMath::ACos(TMath::Cos(jtphi[leaditer]-jtphi[subleaditer]))>2.7)){

                    // printing A value and saving it
                    double Aval = (jtcorrpt[probeiter]-jtcorrpt[tagiter])/(jtcorrpt[probeiter]+jtcorrpt[tagiter]);
                    // removing outlier A values outside of -0.7 to 0.7
                    if(TMath::Abs(Aval)>0.7){continue;}

                    // saving probe and tag momenta
                    htagjtpt->Fill(jtcorrpt[tagiter]);
                    hprobejtpt->Fill(jtcorrpt[probeiter]);
                    htagjteta->Fill(jteta[tagiter]);
                    hprobejteta->Fill(jteta[probeiter]);
                    htagjtphi->Fill(jtphi[tagiter]);
                    hprobejtphi->Fill(jtphi[probeiter]);

                    // put Aval in every alpha bin because no third jet pt is alpha of zero
                    for(unsigned int a=0; a<alphabinnum; a++){
                        // pt slices A value filling
                        // each p is a different slice of pt
                        for(unsigned int p=0; p<ptslicenum; p++){
                            // average momentum between the probe and tag jet 
                            // these are sliced originally to get A vs eta for different pt slices
                            double ptavg = (jtcorrpt[probeiter]+jtcorrpt[tagiter])/2;
                            // if the pt avg is within a certain slice range then add it to the pt slice hists
                            if((ptavg>ptlow[p])&&(ptavg<pthigh[p])){
                                // ptslicesA are A vs eta_probe hists for different pt ranges or slices
                                ptslicesA[a][p]->Fill(jteta[probeiter],Aval);
                            }
                        }
                        // eta slices A value filling
                        // each q is a different slice of eta
                        for(unsigned int q=0; q<etaslicenum; q++){
                            // ptavg is the x axis in one desired type of plot
                            double ptavg = (jtcorrpt[probeiter]+jtcorrpt[tagiter])/2;
                            // if the eta is within a certain slice range then add it to the eta slice hists
                            if((jteta[probeiter]>etalow[q])&&(jteta[probeiter]<etahigh[q])){
                                // etaslicesA are A vs pt_avg hists for different eta ranges
                                etaslicesA[a][q]->Fill(ptavg,Aval);
                            }
                        }
                    }
                }
            }  
        }
    }
    // closing the file I'm getting the information from
    fi->Close();

    // making a new file to store all the histograms of interest in
    TFile *f1 = new TFile(output,"recreate");
    f1->cd();

    // writing the base variable hists to this new file
    save_h1d(hvz, "hvz");
    save_h1d(hvz_, "hvz_");
    save_h1d(htrig, "htrig");
    save_h1d(hjtcorrpt, "hjtcorrpt");
    save_h1d(hjtphi, "hjtphi");
    save_h1d(hjteta, "hjteta");
    save_h1d(hrawpt, "hrawpt");

    // writing probe and tag momenta
    save_h1d(htagjtpt, "htagjtpt");
    save_h1d(hprobejtpt, "hprobejtpt");
    save_h1d(htagjteta, "htagjteta");
    save_h1d(hprobejteta, "hprobejteta");
    save_h1d(htagjtphi, "htagjtphi");
    save_h1d(hprobejtphi, "hprobejtphi");

    // looping through the alpha bins
    for(unsigned int a=0; a<alphabinnum; a++){

        // looping through the pt slices
        for(unsigned int p=0; p<ptslicenum; p++){

            // saving the A vs eta plots for each pt slice
            TString bhtitle0 = Form("A_ptslice_%.0f_%.0f__alpha_%.0f",ptlow[p],pthigh[p],alphas[a]*100);
            save_h2d(ptslicesA[a][p], bhtitle0);

            // looping through the eta slices
            for(unsigned int q=0; q<etaslicenum; q++){

                // conditional below acts like an separated etanum loop 
                if(p==0){
                    // saving eta slice hists of A vs pt
                    if(etalow[q]<0){
                        TString bhtitle1 = Form("A_etaslice_%.0f_%.0f_n__alpha_%.0f",TMath::Abs(etalow[q]*10),TMath::Abs(etahigh[q]*10),alphas[a]*100);
                        save_h2d(etaslicesA[a][q], bhtitle1);
                    }
                    if(etalow[q]>0||etalow[q]==0){
                        TString bhtitle1 = Form("A_etaslice_%.0f_%.0f__alpha_%.0f",etalow[q]*10,etahigh[q]*10,alphas[a]*100);
                        save_h2d(etaslicesA[a][q], bhtitle1);
                    }
                }

                // getting y projection or slice of each eta bin for each pt slice
                etaslices_of_ptslicesA[a][p][q] = ptslicesA[a][p]->ProjectionY("",q,q,"");

                // getting y projection or slice of each pt bin for each eta slice
                ptslices_of_etaslicesA[a][q][p] = etaslicesA[a][q]->ProjectionY("",p,p,"");
                
                // saving the projection hists
                // eta bins of pt slices
                if(etalow[q]<0){
                    TString htitle1 = Form("A_ptslice_%.0f_%.0f_etabin_%.0f_%.0f_n__alpha_%.0f",ptlow[p],pthigh[p],TMath::Abs(etalow[q]*10),TMath::Abs(etahigh[q]*10),alphas[a]*100);
                    TString htitle2 = Form("A_etaslice_%.0f_%.0f_n_ptbin_%.0f_%.0f__alpha_%.0f",TMath::Abs(etalow[q]*10),TMath::Abs(etahigh[q]*10),ptlow[p],pthigh[p],alphas[a]*100);
                    save_h1d(etaslices_of_ptslicesA[a][p][q], htitle1);
                    save_h1d(ptslices_of_etaslicesA[a][q][p], htitle2);
                }
                if(etalow[q]>0||etalow[q]==0){
                    TString htitle1 = Form("A_ptslice_%.0f_%.0f_etabin_%.0f_%.0f__alpha_%.0f",ptlow[p],pthigh[p],etalow[q]*10,etahigh[q]*10,alphas[a]*100);
                    TString htitle2 = Form("A_etaslice_%.0f_%.0f_ptbin_%.0f_%.0f__alpha_%.0f",etalow[q]*10,etahigh[q]*10,ptlow[p],pthigh[p],alphas[a]*100);
                    save_h1d(etaslices_of_ptslicesA[a][p][q], htitle1);
                    save_h1d(ptslices_of_etaslicesA[a][q][p], htitle2);
                }
            }
        }
    }
    f1->Close();
}
