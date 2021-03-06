#include<iostream>
#include<TChain.h>
#include<TH1D.h>
#include<TPad.h>
#include<TCanvas.h>
#include<TString.h>
#include<TRegexp.h>
#include<TStyle.h>
#include<TFile.h>
#include"feio.h"

using namespace std;
using namespace RICHfrontend;

class goodRICHEvent:public RICHEvent{
  public:
	goodRICHEvent(int);
	~goodRICHEvent();
	void Fill(rawEvent&);
	void FindPedestals();
	void FillPedSubtracted(rawEvent&);
  private:
	TH1D* h1[NCHANNELS], *h1PedSub[NCHANNELS];
	uint8_t ipix[64];
	int nasic;
	int ped[NCHANNELS];
};

goodRICHEvent::goodRICHEvent(int _nasic):nasic(_nasic){
	std::copy(chan2pix, chan2pix+64, ipix);

     for(int ich=0;ich<NCHANNELS;ich++) {
		h1[ich]  = new TH1D(Form("h1_%03d",ich), Form("Channel %d, pixel %d; ADC",ich,ipix[ich%64]), 4100, -0.5, 4099.5 );
		h1PedSub[ich]  = new TH1D(Form("h1PedSub_%03d",ich), Form("Pedestal Subtracted Channel %d, pixel %d; ADC - pedestal",ich,ipix[ich%64]), 4100, -50.5, 4049.5 );
     }
}

void goodRICHEvent::Fill(rawEvent &rev)
{

	RICHEvent::Fill(rev);
	for(int ich=0;ich<NCHANNELS;ich++) {
		h1[ich]->Fill(fadc[ich]);
	}

}

void goodRICHEvent::FindPedestals()
{
	// find pedestal for subtraction
	for(int ich=0;ich<NCHANNELS;ich++) {
	        ped[ich] = h1[ich]->GetBinCenter(h1[ich]->GetMaximumBin());
	}
}

void goodRICHEvent::FillPedSubtracted(rawEvent &rev)
{

	RICHEvent::Fill(rev);	
        for(int ich=0;ich<NCHANNELS;ich++) {
                h1PedSub[ich]->Fill(fadc[ich]-ped[ich]);
        }

}


goodRICHEvent::~goodRICHEvent(){
	TCanvas* c1 = new TCanvas("c1","c1",800,800);
	c1->SetGrid();
	c1->Divide(1,2,.0001,.0001);
	c1->Print("adc_plots.pdf[");

	for(int ich=0;ich<NCHANNELS;ich++){
		c1->cd(1)->SetLogy();
		gPad->SetBottomMargin(0);
		gPad->SetGrid();
	     h1[ich]->GetXaxis()->SetRange(h1[ich]->FindFirstBinAbove(1), h1[ich]->FindLastBinAbove(1));
		h1[ich]->Draw();

		c1->cd(2)->SetTopMargin(0);
		gPad->SetGrid();
		TH1* hh = (TH1*) h1[ich]->Clone("hh");
	     hh->GetXaxis()->SetRange(h1[ich]->GetMaximumBin()+50, h1[ich]->FindLastBinAbove(1));
		hh->SetMaximum(hh->GetMaximum()*1.15);
	     hh->GetXaxis()->SetRange(h1[ich]->FindFirstBinAbove(1), h1[ich]->FindLastBinAbove(1));
		hh->Draw();

//		cout << NCHANNELS << endl;
		
//		exit(0);

		c1->Print("adc_plots.pdf");
		delete hh;

		if(ich==63 && nasic==2) ich+=64;
	}
	//gStyle->SetOptStat(0);
	c1->Print("adc_plots.pdf]");

	TFile *fout = new TFile("adc_plots.hist.root", "recreate");
	for(int ich=0;ich<NCHANNELS;ich++) {
                h1[ich]->Write();
		h1PedSub[ich]->Write();
	}
	fout->Close();

	for(int ich=0;ich<NCHANNELS;ich++) {
		delete h1[ich];
		delete h1PedSub[ich];
	}
}



//////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
 rawEvent rawEv;
 goodRICHEvent ev(TString(argv[1]).Contains("2ASIC") ? 2 : 3);




 TChain *tt = new TChain("h22");
 for(int iarg=1;iarg<argc;iarg++)
	if(TString(argv[iarg]).Contains(TRegexp(".root$"))) tt->AddFile(argv[iarg]);

 tt->SetBranchAddress("trigID", &rawEv.trigID);
 tt->SetBranchAddress("timeStamp", &rawEv.timeStamp);

 tt->SetBranchAddress("fadc", rawEv.fadc);

 tt->SetBranchAddress("fnedge", &rawEv.fnedge);
 tt->SetBranchAddress("ftdc", rawEv.ftdc);
 tt->SetBranchAddress("fchan", rawEv.fchan);
 tt->SetBranchAddress("fpolar", rawEv.fpolar);


 // fill initial histograms
 int nen = tt->GetEntries();
 for(int ien=0; ien<nen; ien++){
	tt->GetEntry(ien);
	ev.Fill(rawEv);
 }

 ev.FindPedestals();

 // fill pedestal subtracted histograms
 for(int ien=0; ien<nen; ien++){
        tt->GetEntry(ien);
        ev.FillPedSubtracted(rawEv);
 }

//	cout << "Here at all  ??? " << endl;



 return 0;
}

