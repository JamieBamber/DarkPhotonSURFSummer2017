#define BkgTree_cxx
#include "BkgTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "Riostream.h"

/*
##########
Script to produce a Histogram of the Gamma+Jet background weighted selected events, after applying the ggH cuts detailed in Gabrielli et al. (2016)
vs the lower bound of a Higgs pT cut also applied. (i.e. a cut Higgs pT > [x value] is applied)

used with Higgs_PT_maker.sh

© Jamie Bamber 2017
##########
*/

void BkgTree::Loop(Int_t index=0)
{
//   In a ROOT session, you can do:
//      root> .L BkgTree.C
//      root> BkgTree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   // Create Variables
   Float_t PT;				// photon transverse momentum 
   Float_t E;				// Missing transverse energy (MET)
   Float_t Lumin=100;			// luminosity in fb^-1 (standard is 1pb^-1)
   Float_t E_0=50;
   Float_t PT_0;
   Float_t Eta;				// Psuedo-rapidity of the photon
   Float_t Phi;				// difference in azimuthal angle from the photon to the MET
   Float_t MT;				// Transverse momentum variable
   Float_t DR;				// Delta R displacement variable 
   Float_t Pi = TMath::Pi();
   Int_t N = nentries; 		// number of BkgTrees
   Int_t nbins = 50; // number of bins
   Float_t bin_width = 10;
   //
   //Define Histogram 
   TH1F* h_bkg = new TH1F("h_bkg", "Higgs p_{T}, signal and #gamma+jet background", nbins, 0, nbins*bin_width);
   //
   PT_0 = E_0;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;
		
		// My calculations 
		E = MET;
		PT = PhotonPt; 
		Eta = TMath::Abs(PhotonEta);
		//
		TVector3 photon_v;
		photon_v.SetPtEtaPhi(PT,Eta,PhotonPhi);
		TVector3 MET_v;
		MET_v.SetPtEtaPhi(E,0,METPhi);
		//
		TVector3 Higgs_v = photon_v + MET_v;
		Float_t HPT = Higgs_v.Perp();
		
		if ( ((Eta > 1.4442) && (Eta < 1.566)) || (Eta > 2.5)) continue; // least stringent eta cut
		
		if (Eta > 1.44) continue; // more stringent eta cut
		
		if (PT < PT_0) continue;	// PT criteria 
	
		if (E < E_0) continue;	// MET criteria 	
		
		Phi = PhotonPhi - METPhi;
		 if (Phi > Pi) {
			 Phi = Phi - 2*Pi;
		 } else if (Phi < -Pi) {
			 Phi = Phi + 2*Pi;
		 }
		MT = sqrt(2*PT*E*(1-cos(Phi)));	
		
		if ((MT < 100) || (MT > 130)) continue;	// Transverse Mass cut 
		
		if (Iso_lepton2==1) continue; 				// Isolated lepton criterion
			 
		//h_bkg->Fill(HPT, weight*1000*Lumin);						// Fill Histogram
		Int_t min_bin;
		if ((HPT >= 0) && (HPT <= nbins*bin_width)) {
			min_bin = TMath::FloorNint(HPT/bin_width)+1;
		} else if (HPT > nbins*bin_width) {
			min_bin = nbins + 1;
		}
		for(Int_t i=1; i<=nbins+1; i++) {
			if (i <= min_bin) {	
				h_bkg->SetBinContent(i,h_bkg->GetBinContent(i)+1000*Lumin*weight);
			}
		}
		// progress indicator
		 Float_t percent;
		 Int_t number;
		 percent = static_cast<Float_t>((jentry+1)*100)/nentries;
		 number = int(10*percent);
		 if ((number % 10) == 0) {
			 printf("\r%.0f%% complete",percent);
			 cout.flush();
		}
   }
   cout << endl;
   //
   h_bkg->GetXaxis()->SetTitle("Higgs p_{t} / GeV");
   h_bkg->GetXaxis()->CenterTitle();
   h_bkg->GetYaxis()->SetTitle("frequency");
   h_bkg->GetYaxis()->CenterTitle();
   //
   /*cout << "Background:" << endl;
   for(Int_t i=1; i<=(nbins+1); i++) {
   	   cout << h_bkg->GetBinContent(i) << endl;
   }*/
   //
   TFile *hfile = new TFile("Higgs_PT.root", "RECREATE");
   h_bkg->Write();
   hfile->Close();
}
