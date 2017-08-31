#define Delphes_cxx
#include "Delphes.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "Riostream.h"

void Delphes::Loop(const Int_t index)
{
//   In a ROOT session, you can do:
//      root> .L event.C
//      root> event t
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
   //
   Float_t Lumin=100;
   Float_t sigma=0.8824*1000; // cross section in fb for the ZH channel
   Float_t BR=0.01;
   Float_t Event_weight = Lumin*sigma*BR/50000;
   //
   Float_t Eta;				// Psuedo-rapidity of the photon
   Float_t Phi;				// difference in azimuthal angle from the photon to the MET
   Float_t MT;				// Transverse momentum variable
   Float_t Phi_l;			// lepton phi
   Float_t Eta_l; 			// lepton eta
   Float_t DR;				// Delta R displacement variable 
   //
   Float_t Copy_Electron_PT[kMaxElectron]; // copy of Electron_PT array
   Float_t Copy_Muon_PT[kMaxMuon]; // copy of Muon_PT array
   Int_t e1_index;					// electron1 index
   Int_t e2_index;					// electron2 index
   Int_t m1_index;					// muon1 index
   Int_t m2_index;					// muon2 index
   Float_t Me=0.5109989461*0.001;	// electron invarient mass in GeV
   Float_t Mm=105.6583745*0.001;	// muon invarient mass in GeV
   Float_t Mz=91.1876;				// Z boson invarient mass in GeV
   //TLorentzVector l1;
   //TLorentzVector l2;
   TLorentzVector ZCandidate;
   TVector3 photon_vec;
   TVector3 MET_vec;
   Float_t dPhi_ll_EE;				// azimuthal angle between the sum of the lepton Pt and the sum of the MET and photon Pt
   Float_t dPhi_ll;					// azimuthal angle between leptons
   Float_t p_variable;
   Float_t DeltaM;
   //
   Int_t lep1Type = 0;
	Int_t lep1PassSelection = 0;
	Int_t lep1Pt = -999;
	Float_t lep1Eta = -999;
	Float_t lep1Phi = -999;  
	Int_t lep2Type = 0;
	Int_t lep2PassSelection = 0;
	Float_t lep2Pt = -999;
	Float_t lep2Eta = -999;
	Float_t lep2Phi = -999;
	Float_t dileptonMass = -999;
	Float_t lep1MT = -999;
	Float_t lep1GenMetMT = -999;
   //
   Float_t Pi = TMath::Pi();
   Int_t N = nentries; 		// number of events 
   Float_t N_s[12];			// number of selected events (array for different cuts)
   Float_t Eff;				// efficiency
   //
   for(int n=0; n<=11; n++) {
   	   N_s[n] = 0;
   }
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (ientry < 0) break;
		// if (Cut(ientry) < 0) continue;
		
		// My calculations 
		E = MissingET_MET[0];
		PT = TMath::MaxElement(Photon_,Photon_PT);
		Eta = TMath::Abs(Photon_Eta[TMath::LocMax(Photon_,Photon_PT)]);
		
		// #### Trigger cuts
		//if ( (((Eta > 1.4442) && (Eta < 1.566)) || (Eta > 2.5)) && (index==0) ) continue; // least stringent eta cut
		
		//if ( (Eta > 1.44) && (index>=1) ) continue; // barrel region eta cut
		
		if ( (PT < 20) && (index>=2) ) continue;	// Photon PT cut ### 20 value taken from CMS paper
	
		if ( (E < 60) && (index>=3) ) continue;	// MET cut	### 60 value taken from CMS paper
		
		Phi = Photon_Phi[TMath::LocMax(Photon_,Photon_PT)] - MissingET_Phi[0];
		if (Phi > Pi) {
			Phi = Phi - 2*Pi;
		} else if (Phi < -Pi) {
			 Phi = Phi + 2*Pi;
		}
		MT = sqrt(2*PT*E*(1-cos(Phi)));	
		
		if ( ((MT < 100) || (MT > 130)) && (index>=4) ) continue;	// Transverse momentum criterion ### kept over from Gabrielli et al. paper
		
		//
		
		// ## further cuts also from CMS paper
		
		if( index>=5 ) {
			// ## lepton selection cuts
			int razorbox = 0;
			//-------------------------------
			//1) Look for Zmm Candidate
			//-------------------------------
			double bestDiMuon_PT = -1;
			for( int i = 0; i < Muon_; i++ )	{
				//if(!isVetoMuon(i)) continue;  
				if(Muon_PT[i] < 20) continue;
				if(abs(Muon_Eta[i]) > 2.4) continue;
				for( int j = i+1; j < Muon_; j++ )	{
					//if(!isVetoMuon(j)) continue;  
					if(Muon_PT[j] < 20) continue;
					if(abs(Muon_Eta[j]) > 2.4) continue;
				
					TLorentzVector tmpMuon1;
					tmpMuon1.SetPtEtaPhiM(Muon_PT[i],Muon_Eta[i], Muon_Phi[i],0.1057);
					TLorentzVector tmpMuon2;
					tmpMuon2.SetPtEtaPhiM(Muon_PT[j],Muon_Eta[j], Muon_Phi[j],0.1057);
					double tmpMass = (tmpMuon1+tmpMuon2).M();	    
					double tmpDileptonPt = (tmpMuon1+tmpMuon2).Pt();
					
					//if ( _debug ) cout << "Zmm candidate: " << tmpMass << " " << tmpDileptonPt << "\n";
					
					if ( tmpMass > (Mz-15) && tmpMass < (Mz+15) && tmpDileptonPt > bestDiMuon_PT)  {
						bestDiMuon_PT = tmpDileptonPt;
						razorbox = 1; // 1 in place of Zmm
						lep1Type = 13 * -1 * Muon_Charge[i];
						lep1Pt = Muon_PT[i];
						lep1Eta = Muon_Eta[i];
						lep1Phi = Muon_Phi[i];
						//lep1PassSelection = 1 + 2 * isTightMuon(i);
						lep2Type = 13 * -1 * Muon_Charge[j];
						lep2Pt = Muon_PT[j];
						lep2Eta = Muon_Eta[j];
						lep2Phi = Muon_Phi[j];
						//lep2PassSelection = 1 + 2 * isTightMuon(j);
						dileptonMass = tmpMass;
						ZCandidate = tmpMuon1 + tmpMuon2;
					
						//for MC apply lepton eff scale factor
						/*if (!isData ) {
							if ( matchesGenMuon(lep1Eta,lep1Phi)) leptonEffSF *=  helper->getVetoMuonScaleFactor( lep1Pt, lep1Eta, true);		
							if ( matchesGenMuon(lep2Eta,lep2Phi)) leptonEffSF *=  helper->getVetoMuonScaleFactor( lep2Pt, lep2Eta, true);			
						}*/
					}
				}
			}
		
		
			//-------------------------------
			//2) Look for Zee Candidate
			//-------------------------------
			if (razorbox == 0) {
				double bestDielectronPt = -1;
				for( int i = 0; i < Electron_; i++ )	{
					//if(!isVetoElectron(i)) continue;  
					if(Electron_PT[i] < 20) continue;
					if(abs(Electron_Eta[i]) > 2.4) continue;
					for( int j = i+1; j < Electron_; j++ )	{
						//if(!isVetoElectron(j)) continue;  
						if(Electron_PT[j] < 20) continue;
						if(abs(Electron_Eta[j]) > 2.4) continue;
					
						TLorentzVector tmpElectron1;
						tmpElectron1.SetPtEtaPhiM(Electron_PT[i],Electron_Eta[i], Electron_Phi[i],0.000511);
						TLorentzVector tmpElectron2;
						tmpElectron2.SetPtEtaPhiM(Electron_PT[j],Electron_Eta[j], Electron_Phi[j],0.000511);
						double tmpMass = (tmpElectron1+tmpElectron2).M();	    
						double tmpDileptonPt = (tmpElectron1+tmpElectron2).Pt();
				
						//if ( _debug ) cout << "Zee candidate: " << tmpMass << " " << tmpDileptonPt << "\n";
				
						if ( tmpMass > (Mz-15) && tmpMass < (Mz+15) && tmpDileptonPt > bestDielectronPt)  {
							bestDielectronPt = tmpDileptonPt;
							razorbox = 2; // 2 in place of Zee
							lep1Type = 11 * -1 * Electron_Charge[i];
							lep1Pt = Electron_PT[i];
							lep1Eta = Electron_Eta[i];
							lep1Phi = Electron_Phi[i];
							//lep1PassSelection = 1 + 2 * isTightElectron(i);
							lep2Type = 11 * -1 * Electron_Charge[j];
							lep2Pt = Electron_PT[j];
							lep2Eta = Electron_Eta[j];
							lep2Phi = Electron_Phi[j];
							//lep2PassSelection = 1 + 2 * isTightElectron(j);
							dileptonMass = tmpMass;
							ZCandidate = tmpElectron1 + tmpElectron2;
				
							//for MC apply lepton eff scale factor
							/*if (!isData ) {
								if ( matchesGenElectron(lep1Eta,lep1Phi)) leptonEffSF *=  helper->getVetoElectronScaleFactor( lep1Pt, lep1Eta, true);		
								if ( matchesGenElectron(lep2Eta,lep2Phi)) leptonEffSF *=  helper->getVetoElectronScaleFactor( lep2Pt, lep2Eta, true);			
							}*/
						}
					}
				}
			}
			if(razorbox == 0) continue;
		}
		
		if ( (index>=6) && (TMath::MaxElement(Jet_,Jet_PT) > 30) ) continue; // no jets with PT greater than 30GeV cut
		
		if (index>=7) {														// first Delta Phi cut
			photon_vec.SetPtEtaPhi(PT,Photon_Eta[TMath::LocMax(Photon_,Photon_PT)],Photon_Phi[TMath::LocMax(Photon_,Photon_PT)]);
			MET_vec.SetPtEtaPhi(E,0,MissingET_Phi[0]);
			dPhi_ll_EE = TMath::Abs( (photon_vec + MET_vec).Phi() - ZCandidate.Phi() );
			if ( dPhi_ll_EE > Pi) dPhi_ll_EE = 2*Pi - dPhi_ll_EE;
			if ( dPhi_ll_EE < 2.7) continue;
		}	
		
		if (index>=8) {		// pt variable cut
			p_variable = TMath::Abs( (photon_vec+MET_vec).Pt() - ZCandidate.Pt() )/ZCandidate.Pt();
			if (p_variable > 0.5) continue;
		}
		
		if (index>=9) {		// azimuthal angle between leptons cut
			dPhi_ll = TMath::Abs(lep1Phi-lep2Phi);
			if ( dPhi_ll > Pi ) dPhi_ll = 2*Pi - dPhi_ll;
			if ( dPhi_ll > 2.25) continue;
		}
		
		if (index>=10) {		// lepton combined pT cut
			if ( ZCandidate.Pt() < 60) continue;
		}
		
		N_s[index] = N_s[index] + 1;								// increment element by weight
		fail:;													// goto this point if the event fails the test
		// progress indicator
		if (((jentry+1)*100) % nentries == 0) {
			Float_t percent;
			percent = static_cast<Float_t>((jentry+1)*100)/nentries;
			printf("\r%.0f%% complete",percent);
			cout.flush();
		}
   }
   cout << endl;
   // Print out data in csv format
   // using eff = no. selected / no. events
   for(int n=index; n<=index; n++) {
   	   cout << N_s[n] << ", " << N_s[n]*Event_weight << ", " << 100*N_s[n]/N << ", " << endl;
   }
}
