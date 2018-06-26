#include<iostream>
#include<list>
#include<TFile.h>
#include<TTree.h>
#include<TF1.h>
#include<TH2F.h>
#include<TH3F.h>
#include<TPad.h>
#include<TCanvas.h>
#include<TString.h>
#include<TRegexp.h>
#include<TStyle.h>
#include<TProfile.h>
#include<TCutG.h>
#include"feio.h"

using namespace RICHfrontend;
using namespace std;

class goodRICHEvent:public RICHEvent{
  public:
	goodRICHEvent(int);
	~goodRICHEvent();
	void Fill(rawEvent&);

	void Timewalk_Calabration(TTree*);

  private:
	TH3F* h_dt_t0_adc[NCHANNELS];

	TFile* f1;

	TTree* hh_tree[NCHANNELS];
	TTree* hh_tree_all;

	int nasic;

	Int_t adc, delta_t, t0;



};



Double_t fitf(Double_t *x,Double_t *par);

goodRICHEvent::goodRICHEvent(int _nasic):nasic(_nasic){


	f1 = new TFile ("temp.root", "RECREATE");

	for(int ich=0;ich<NCHANNELS;ich++){

		h_dt_t0_adc[ich]  = new TH3F(Form("h_dt_t0_adc_%03d",ich), Form("Channel %d, pixel %d;TDC hit, duration;TDC hit, leading time;ADC",ich,chan2pix[ich%64]),
			100,0.5,100.5,		120, 90.5, 250.5,		380,250.5,2150.5);
//			100,0.5,100.5,		275, 100.5, 1200.5,		340,450.5,2150.5);


 		hh_tree[ich] = new TTree(Form("tree_%03d",ich), "Tree Tree");
// 
 		hh_tree[ich]->Branch("adc",     &adc,     "adc/I");
       	hh_tree[ich]->Branch("delta_t", &delta_t, "delta_t/I");
       	hh_tree[ich]->Branch("t0",      &t0,      "t0/I");	

	}	


// 
// //	    TFile *f = new TFile("temp.root");
// 
// 		hh_tree_all = new TTree("tree_all", "Tree Tree");
// 
//  		hh_tree_all->Branch("adc",     &adc,     "adc/I");
//        	hh_tree_all->Branch("delta_t", &delta_t, "delta_t/I");
//        	hh_tree_all->Branch("t0",      &t0,      "t0/I");	
//  
//  	    hh_tree_all->SetBranchStatus("*",1);
// 

}


void goodRICHEvent::Fill(rawEvent &rev)
{



	RICHEvent::Fill(rev);

//	for(int ichan=0; ichan<NCHANNELS; ichan++){

	for(int ichan=0; ichan<1; ichan++){

		for(int iedge=0; iedge<ftdc[ichan].size(); iedge++)

		if(fpolar[ichan][iedge]==fLeadingEdge
			&& iedge<ftdc[ichan].size()-1
			&& fpolar[ichan][iedge+1]==fTrailingEdge
			){

				ftime[ichan].push_back(ftdc[ichan][iedge]);
				fdur[ichan].push_back(ftdc[ichan][iedge+1] - ftdc[ichan][iedge]);
				iedge++;


//				cout << iedge << endl;
//				cout << ftdc[ichan][iedge+1] << endl;
//				exit(0);

		}


		if(ftime[ichan].size()>0) {
			h_dt_t0_adc[ichan]->Fill(fdur[ichan][0], ftime[ichan][0], fadc[ichan]);


			




//			cout << "asdadad " <<  ftime[ichan][0] << endl; 

//			exit(0);

 			t0  = ftime[ichan][0];
 			delta_t  = fdur[ichan][0];
 			adc = fadc[ichan];


		} else {
			h_dt_t0_adc[ichan]->Fill(-1, -1, fadc[ichan]);

  			adc = fadc[ichan];
 			t0  = -1;
 			delta_t  = -1;

		}

//   		cout << adc << "   " << delta_t << "   "<< t0 << endl; 


//		cout << adc << "   " << delta_t << "   "<< t0 << endl; 

//		cout << adc << "    " << t0 << "     " << delta_t << endl;

//		if (adc > 0 && t0 >0 && delta_t >0)

		hh_tree[ichan]->Fill();

//		hh_tree_all->Fill();


//		exit(0);

//		exit(0);

	}

}


goodRICHEvent::~goodRICHEvent(){
	TString pdfname("adc_vs_t0.pdf");
	TCanvas* c1 = new TCanvas("c1","c1",1200,1200);
//	c1->Divide(1,3,.0001,.0001);

	c1->Divide(2,2,.0001,.0001);

	c1->Print(pdfname+"[");


	TCanvas* c2 = new TCanvas();


	for(int ich=0; ich<1; ich++){


//	for(int ich=0; ich<NCHANNELS; ich++){

		TH3F* h3 = h_dt_t0_adc[ich];

		c1->cd(1);
		gPad->SetTopMargin(0);
		gPad->SetGrid();

		TH1* hyz = h3->Project3D("yz");


		hyz->Draw("colz");

//////////////
		c1->cd(2);
		gPad->SetTopMargin(0);
		gPad->SetGrid();

		TH1* hyx = h3->Project3D("yx");
		TF1* f1 = new TF1("f1", "[0]+expo(1)",0,100);
		f1->SetParameters(100,4,-0.04);
//		hyx->Fit(f1,"Q0");

		f1->SetParameter(0,f1->GetParameter(0)+7);
		hyx->Draw("colz");
//		f1->Draw("same");

//////////////
		c1->cd(3);
		gPad->SetTopMargin(0);
		gPad->SetGrid();

		TH1* hadc0 = (TH1*) h3->Project3D("hadc_z0");
		hadc0->SetLineStyle(2);
		hadc0->SetLineWidth(2);
		TH1* hadc1 = (TH1*) h3->Project3D("z1 NUF NOF");
		double hmax = hadc1->GetMaximum();

		hadc1->GetXaxis()->SetRange(hadc0->FindFirstBinAbove(hmax/50-50/hadc0->GetBinWidth(1)), hadc1->FindLastBinAbove(hmax/50));
		hyz->GetXaxis()->SetRange(hadc0->FindFirstBinAbove(hmax/50-50/hadc0->GetBinWidth(1)), hadc1->FindLastBinAbove(hmax/50));


		TH1* hadc2 = (TH1*) h3->Project3D("z2 NUF NOF");
		hadc2->Reset();
		hadc2->SetLineColor(kRed);
		for(int ibx=1;ibx<=h3->GetNbinsX();ibx++)
		for(int iby=1;iby<=h3->GetNbinsY();iby++)
			if(h3->GetYaxis()->GetBinCenter(iby) > f1->Eval(h3->GetXaxis()->GetBinCenter(ibx)))
			for(int ibz=1;ibz<=h3->GetNbinsZ();ibz++)
			for(int ien=0;ien<h3->GetBinContent(ibx,iby,ibz);ien++)
				hadc2->Fill(hadc2->GetBinCenter(ibz));
// 
// 		hadc1->Draw();
// //		for(int ibin=1;ibin<=hadc0->GetNbinsX();ibin++)
// //			std::cout<<"pavel:"<<ich<<":"<<ibin<<" "<<hadc0->GetBinContent(ibin)<<std::endl;
// 		hadc0->Draw("same");
// 		hadc1->Draw("same");
// 		hadc2->Scale(hmax/hadc2->GetMaximum());
// 		hadc2->Draw("same hist");



		c1->cd(4);
		gPad->SetTopMargin(0);
		gPad->SetGrid();

		TH1* hxz = h3->Project3D("xz");
		hxz->Draw("colz");

		hxz->GetXaxis()->SetRange(hadc0->FindFirstBinAbove(hmax/50-50/hadc0->GetBinWidth(1)), hadc1->FindLastBinAbove(hmax/50));


//		TProfile *prof = hxz->ProfileX();
	
//		TF1*f3 = new TF1("f3", fitf, 500, 1100, 3);

// 		f3->SetParameter(0, 500);
// 		f3->SetParameter(0, 500);
// 		f3->SetParameter(2, 60);
// 

// 	    hxz->Fit("f3", "R");

//   		prof->Fit("f3");

//////////////

		Timewalk_Calabration(hh_tree[ich]);



		c1->Print(pdfname);
		c1->Print("fit_test.root");

		if(ich==63 && nasic==2) ich+=64;

		delete hyz, hyx, hadc0, hadc1, hadc2, f1, h3;


//		break;


	}
	c1->Print(pdfname+"]");

//	c2->Print("fit_test.root");


	for(int ich=0; ich<NCHANNELS; ich++)
	if(h_dt_t0_adc[ich])
		delete h_dt_t0_adc[ich];
}


//////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
 rawEvent rawEv;
 goodRICHEvent ev(TString(argv[1]).Contains("2ASIC") ? 2 : 3);



 for(int iarg=1;iarg<argc;iarg++){
	if(!TString(argv[iarg]).Contains(TRegexp(".root$"))) continue;
	TString bname(TString(argv[iarg])(TRegexp("thr[0-9]*")));

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

	delete tt, ff;
 }

 return 0;
}




Double_t fitf(Double_t *x, Double_t *par) {
 
    Double_t xx = 0;

//    Double_t fitval = sqrt(xx-par[0])/par[1] + par[2]; 

    Double_t fitval =  par[0] * xx - par[1] + par[2]; 
 
    return fitval;

}








/*--------------------------------------------------*/
/*--------------------------------------------------*/


void goodRICHEvent::Timewalk_Calabration(TTree* t1) {

	TCanvas* c10 = new TCanvas("c10", "c10", 1200, 1200);
	c10-> Divide(3,2);

	

	TCutG *cutg = new TCutG("mycut", 10);

	cutg->SetVarX("adc");
	cutg->SetVarY("delta_t");

	cutg->SetPoint(0, 500,  15);
	cutg->SetPoint(1, 530,  32);
	cutg->SetPoint(2, 600,  44);
	cutg->SetPoint(3, 700,  52);
	cutg->SetPoint(4, 800,  55);
	cutg->SetPoint(5, 900,  57);
	cutg->SetPoint(6, 1150, 57);
	cutg->SetPoint(7, 1150, 50);
	cutg->SetPoint(8, 800,  40);
	cutg->SetPoint(9, 520,  15);
	cutg->SetPoint(10, 500, 15);



	c10->cd(1);

	/*--------------------------------------------------*/
	/// Distribution all

	t1->Draw("adc:delta_t", "adc < 1300 && adc > 400 && delta_t > 10 && delta_t < 70", "colz");
	TH1 *dt_adc = (TH1*)gPad->GetPrimitive("htemp"); 

	dt_adc->Draw("colz");
	

	/*--------------------------------------------------*/
	/// Distribution with narrow cut

	c10->cd(2);
	t1->Draw("adc:delta_t >> delta_adc(60, 10, 70, 900, 400, 1300)", "adc < 1300 && adc > 400 && delta_t > 10 && delta_t < 70 && mycut", "colz");
	TH1 *dt_adc_cut = (TH1*)gPad->GetPrimitive("delta_adc"); 

	TF1* ff = new TF1("ff", "[0]+expo(1)",15,55);
//	ff->SetParameter(0 , 60);

	TF1* f_linear = new TF1("f_linear", "[0]+[1]*x", 50, 60);

	dt_adc_cut->Fit("ff", "R");

	dt_adc_cut->Fit("f_linear", "R");

	ff->Draw("same");
	
	f_linear->SetLineColor(1);
	
	f_linear->Draw("same");


	
	/*--------------------------------------------------*/
	/// Distribution with narrow cut

	c10->cd(3);

	Double_t y_offset  = ff->GetParameter(0);
	Double_t expo_off  = ff->GetParameter(1);
	Double_t expo_muti = ff->GetParameter(2);

	Double_t corr_offset  = f_linear->GetParameter(0);
	Double_t corr_slope   = f_linear->GetParameter(1);


	TString cor_str;

	cor_str.Form("((adc - (%f + exp(%f + %f * delta_t))) + (delta_t*%f +%f))" , y_offset, expo_off, expo_muti, corr_slope, corr_offset);


//	cor_str.Form("((delta_t - (%f + exp(%f + %f * adc))) + (adc*%f +%f))" , y_offset, expo_off, expo_muti, corr_slope, corr_offset);

	t1->Draw( cor_str + ":delta_t >> low_corr(60, 10, 70, 1300, 0, 1300) ", "adc < 1300 && adc > 400 && delta_t > 10 && delta_t <= 50 && mycut", "colz");

	TH1 *dt_adc_low_cor = (TH1*)gPad->GetPrimitive("low_corr"); 


//	dt_adc_low_cor -> DrawCopy("colz");


	t1->Draw("adc:delta_t >> high_corr(60, 10, 70, 1300, 0, 1300)", "adc < 1300 && adc > 400 && delta_t > 50 && delta_t < 70 && mycut", "colz");

	TH1 *dt_adc_hi_cor = (TH1*)gPad->GetPrimitive("high_corr"); 


	dt_adc_low_cor->Draw("colz");

	dt_adc_hi_cor -> Draw("colz same");




	c10->cd(4);

//	t1->Draw("t0:adc >> time_corr(900, 400, 1300, 120, 120, 240)", "adc < 1300 && adc >= 500", "colz");
	t1->Draw("t0:adc >> time_corr(900, 400, 1300, 40, 130, 170)", "adc < 1300 && adc >= 500", "colz");
	


	c10->cd(5);



	TString cor_str_low;
	TString cor_str_hi;

	cor_str_low.Form( "(%f + exp(%f + %f * delta_t))", y_offset, expo_off, expo_muti);

	t1->Draw("t0:" + cor_str_low + " >> time_corr_low(900, 400, 1300, 40, 130, 170)", "adc < 1300 && adc >= 500", "colz");

//	t1->Draw("t0:adc >> time_corr_hi(900, 400, 1300, 40, 130, 170)", "adc < 1300 && adc >= 500", "colz");



	c10->Print("time_walk/test.png");



	delete c10;


}


