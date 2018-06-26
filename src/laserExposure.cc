#include<iostream>
#include<TFile.h>
#include<TTree.h>
#include<TH2I.h>
#include<TGraph.h>
#include<TCanvas.h>
#include<TStyle.h>
#include"feio.h"

using namespace std;
using namespace RICHfrontend;

class goodRICHEvent:public RICHEvent{
  public:
	goodRICHEvent(int);
	~goodRICHEvent();
	void Fill(rawEvent&);
	void Print();

  private:
	TH1I* h1[NCHANNELS];
	TH2I* htot;
	TCanvas* c1;
	int nasic;
};

goodRICHEvent::goodRICHEvent(int _nasic):nasic(_nasic){
	c1 = new TCanvas("c1", "c1", nasic*400, 400);
	c1->Print("laser_exposure.pdf[");
     for(int ich=0;ich<NCHANNELS;ich++)
		h1[ich]  = new TH1I(Form("h1_%03d",ich), Form("Channel %d; ADC",ich), 4100, -0.5, 4099.5 );

	htot = new TH2I("htot", "laser exposure", nasic*8, 0.5,8*nasic+.5, 8, 0.5, 8.5);
}

void goodRICHEvent::Fill(rawEvent &rev)
{
	RICHEvent::Fill(rev);
	for(int ich=0;ich<NCHANNELS;ich++)
		h1[ich]->Fill(fadc[ich]);
}

void goodRICHEvent::Print()
{
	TH2I* hperfile = new TH2I("hperfile", "laser exposure", nasic*8, 0.5,8*nasic+.5, 8, 0.5, 8.5);

	Float_t value_min = 1000000000;

	for(int ichan=0;ichan<NCHANNELS;ichan++){
		int ipmt = ichan/64;

		if(nasic==2 && ipmt==1) continue;
		if(nasic==2 && ipmt==2) ipmt--;

		int ipx = chan2pix[ichan%64]-1;
		int icol = ipx%8+1 + (nasic-ipmt-1)*8;
		int irow = ipx/8+1;

	    double ww = h1[ichan]->Integral(h1[ichan]->GetMaximumBin()+50, h1[ichan]->GetNbinsX());

		// hperfile->Fill(icol, irow, ww);

		///*--------------------------------------------------*/
		/// Plot the inverted image of the laser map
		int irow_inverted = 9 - irow;


		if (ww > 0) {
			hperfile->Fill(icol, irow_inverted, ww);
//			cout << icol << "    "  <<  irow_inverted << "    " <<  ww << endl;

			if(ww<value_min) {
				value_min = ww;
			}

		}


		h1[ichan]->Reset();
	}
	htot->Add(hperfile);


	gStyle->SetOptStat(0);
//	gStyle->SetPalette(1);

//	hperfile->SetMinimum(1000);

	hperfile->GetZaxis()->SetRangeUser(value_min, hperfile->GetMaximum());
	htot->GetZaxis()->SetRangeUser(value_min, hperfile->GetMaximum());

	hperfile->Draw("colz");
	c1->Print("laser_exposure.pdf");

	delete hperfile;
}

goodRICHEvent::~goodRICHEvent(){
	htot->Draw("colz");
	c1->Print("laser_exposure.pdf");
	c1->Print("laser_exposure.pdf]");

	delete htot, c1;
}



//////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
 rawEvent rawEv;
 goodRICHEvent ev(TString(argv[1]).Contains("2ASIC") ? 2:3);

 for(int iarg=1;iarg<argc;iarg++){
	TFile* ff = new TFile(argv[iarg]);
	TTree* tt = (TTree*) ff->Get("h22");

	tt->SetBranchAddress("trigID", &rawEv.trigID);
	tt->SetBranchAddress("timeStamp", &rawEv.timeStamp);

	tt->SetBranchAddress("fadc", rawEv.fadc);

	tt->SetBranchAddress("fnedge", &rawEv.fnedge);
	tt->SetBranchAddress("ftdc", rawEv.ftdc);
	tt->SetBranchAddress("fchan", rawEv.fchan);
	tt->SetBranchAddress("fpolar", rawEv.fpolar);

	int nen = tt->GetEntries();
	for(int ien=0; ien<nen; ien++){
		tt->GetEntry(ien);
		ev.Fill(rawEv);
	}
	ev.Print();

	delete tt, ff;
 }

 return 0;
}

