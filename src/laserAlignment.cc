#include<iostream>
#include<TFile.h>
#include<TTree.h>
#include<TH1I.h>
#include<TGraphErrors.h>
#include<TCanvas.h>
#include<TStyle.h>
#include"feio.h"

using namespace RICHfrontend;

class goodRICHEvent:public RICHEvent{
  public:
	goodRICHEvent();
	~goodRICHEvent();
	void Fill(rawEvent&);
	double getCenterCounts();

  private:
	TH1I* h1[NCHANNELS];
	uint8_t ipix[64];
};

goodRICHEvent::goodRICHEvent(){
	std::copy(chan2pix, chan2pix+64, ipix);

     for(int ich=0;ich<NCHANNELS;ich++)
		h1[ich]  = new TH1I(Form("h1_%03d",ich), Form("Channel %d, pixel %d; ADC",ich,ipix[ich%64]), 4100, -0.5, 4099.5 );
}

void goodRICHEvent::Fill(rawEvent &rev)
{
	RICHEvent::Fill(rev);
	for(int ich=0;ich<NCHANNELS;ich++)
		h1[ich]->Fill(fadc[ich]);
}

double goodRICHEvent::getCenterCounts()
{
	double tot = 0;
	int centerchan[] = {16, 44, 12, 48};
	for(int ichan=0;ichan<4;ichan++){
		int ich = centerchan[ichan];
	     h1[ich]->GetXaxis()->SetRange(h1[ich]->GetMaximumBin()+20, h1[ich]->GetNbinsX());
		tot += h1[ich]->Integral();
	}
	return tot;
}

goodRICHEvent::~goodRICHEvent(){
	for(int ich=0;ich<NCHANNELS;ich++)
		delete h1[ich];
}



//////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
 rawEvent rawEv;
 TGraphErrors gr;

 for(int iarg=1;iarg<argc;iarg++){
	goodRICHEvent ev;

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
	gr.SetPoint(gr.GetN(), 50+(iarg-1)*10, ev.getCenterCounts());
	gr.SetPointError(gr.GetN()-1, 0, sqrt(ev.getCenterCounts()));

	delete tt, ff;
 }
 TCanvas* c1 = new TCanvas("c1","c1",800,400);
 c1->SetGrid();
 gr.SetTitle("Counts in center MAPMT vs laser Y position;laser Y position;center counts");
 gr.SetMarkerStyle(21);
 gr.SetMarkerSize(1);
 gr.SetMarkerColor(kBlue);
 gr.SetLineColor(kBlue);
 gr.Draw("AP");
 c1->Print("laser_alignment.pdf");

 return 0;
}

