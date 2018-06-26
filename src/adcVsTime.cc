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
#include<TText.h>
#include<TProfile.h>
#include<TCutG.h>
#include<TMath.h>

#include"feio.h"

using namespace RICHfrontend;
using namespace std;

class goodRICHEvent:public RICHEvent{

  public:
	goodRICHEvent(int);
	~goodRICHEvent();
	void Fill(rawEvent&);

	void Timewalk_Calabration(TTree*);
	void Timewalk_Calabration_Raw(TTree*);

  private:
	TH3F* h_dt_t0_adc[NCHANNELS];

	TFile* f1;

	TTree* hh_tree[NCHANNELS];
	TTree* hh_tree_all;

	int nasic;

	Int_t adc, delta_t, t0;

};



///*--------------------------------------------------*/
/// Fitting functions

Double_t fitf(Double_t *x,Double_t *par);
Double_t poissonf(Double_t*x,Double_t*par); 

Double_t gauss_f(Double_t*x,Double_t*par);
Double_t gauss_n_expo(Double_t*x, Double_t*par);

Double_t fit_gaus_fun(Double_t*, Double_t*);
Double_t fit_landau(Double_t*, Double_t*);
Double_t fit_combo(Double_t*, Double_t*);


///*--------------------------------------------------*/
/// Slicing function

void Fit_t0_adc_slice(TTree*);

void Fit_dt_adc_slice(TTree*);

void Fit_adc_dt_slice(TTree* t1);

void Fit_t0_dt_slice(TTree* t1);




/*--------------------------------------------------*/

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

//for(int ichan=0; ichan<NCHANNELS; ichan++){

	for(int ichan=0; ichan<1; ichan++){

		for(int iedge=0; iedge<ftdc[ichan].size(); iedge++)

		if(fpolar[ichan][iedge]==fLeadingEdge
			&& iedge<ftdc[ichan].size()-1
			&& fpolar[ichan][iedge+1]==fTrailingEdge
			){

				ftime[ichan].push_back(ftdc[ichan][iedge]);
				fdur[ichan].push_back(ftdc[ichan][iedge+1] - ftdc[ichan][iedge]);
				iedge++;

		}


		if(ftime[ichan].size()>0) {
			h_dt_t0_adc[ichan]->Fill(fdur[ichan][0], ftime[ichan][0], fadc[ichan]);

 			t0  = ftime[ichan][0];
 			delta_t  = fdur[ichan][0];
 			adc = fadc[ichan];

		} else {
			h_dt_t0_adc[ichan]->Fill(-1, -1, fadc[ichan]);

  			adc = fadc[ichan];
 			t0  = -1;
 			delta_t  = -1;

		}

		hh_tree[ichan]->Fill();

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
		Timewalk_Calabration_Raw(hh_tree[ich]);



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





/*--------------------------------------------------*/
/*--------------------------------------------------*/


void goodRICHEvent::Timewalk_Calabration(TTree* t1) {

	TCanvas* c10 = new TCanvas("c10", "c10", 1200, 1200);
	c10-> Divide(3,3);

	

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
//	t1->Draw("t0:adc >> time_corr(900, 400, 1300, 40, 130, 170)", "adc < 1300 && adc >= 500", "colz");
	
//	t1->Draw("t0:adc >> time_corr(900, 400, 1300, 40, 130, 170)", "adc < 1300 && adc >= 500 && mycut", "colz");


	c10->cd(4);


	t1->Draw("t0:delta_t >> time_walk_org(70, 0, 70, 40, 130, 170)", "adc < 1300 && adc >= 500 && mycut", "colz");

//	t1->Draw("t0:adc >> time_corr_hi(900, 400, 1300, 40, 130, 170)", "adc < 1300 && adc >= 500", "colz");

	TH2D *time_walk_org = (TH2D*)gPad->GetPrimitive("time_walk_org");

	TF1* f_walk_org = new TF1("f_walk_org", "[0] + [1]/(x+[2])", 15, 60);
	time_walk_org->Fit("f_walk_org", "R");
	


	c10->cd(5);


	Float_t time_walk_offset_org = f_walk_org->GetParameter(0);
	Float_t time_top_org = f_walk_org->GetParameter(1);
	Float_t time_bottom_org = f_walk_org->GetParameter(2);

	TString timewalk_org_str;
	timewalk_org_str.Form( " (t0 - (%f + %f/(delta_t+%f)))", time_walk_offset_org, time_top_org, time_bottom_org);



	t1->Draw( timewalk_org_str + ": delta_t >> timewalk_before_org(70, 0, 70, 40, -20, 20)", "adc < 1300 && adc >= 500 && mycut", "colz");

	TH2D* timewalk_before_org = (TH2D*)gPad->GetPrimitive("timewalk_before_org");

	c10->cd(6);


	TH1D* y_proj_org = (TH1D*) timewalk_before_org->ProjectionY("before_py");
	y_proj_org->DrawCopy("hist");

	y_proj_org->GetXaxis()->SetRange(10, 30);

	y_proj_org->Draw("hist");


 	TF1* time_res_org = new TF1("time_res_org", "gaus", -5, 4);
 
 	y_proj_org->Fit("time_res_org", "R");
 
 	time_res_org->Draw("same");
 
 	TText *text1 = new TText();
	text1->DrawTextNDC( 0.65, 0.7, Form("Mean: %.3f", time_res_org->GetParameter(1)));
	text1->DrawTextNDC( 0.65, 0.6, Form("Sigma: %.3f", time_res_org->GetParameter(2)));






   /*--------------------------------------------------*/
   /*--------------------------------------------------*/
   /*--------------------------------------------------*/



	c10->cd(7);


	TString cor_str_low;
	TString cor_str_hi;

	cor_str_low.Form( "((log(adc - %f)- %f) / %f)", y_offset, expo_off, expo_muti);
//	cor_str_low.Form( "delta_t");

//	t1->Draw("t0:" + cor_str_low + " >> time_corr_low(70, 0, 70, 40, 130, 170)", "adc < 550 && mycut", "colz");

	t1->Draw("t0:" + cor_str_low + " >> time_corr_low(70, 0, 70, 40, 130, 170)", "adc < 1300 && adc >= 500 && mycut", "colz");



 	TH2D *time_walk_corr = (TH2D*)gPad->GetPrimitive("time_corr_low"); 
 	
 
 
 
 	TF1* f_walk_cor = new TF1("f_walk_cor", "[0] + [1]/(x+[2])", 10, 60);
 	time_walk_corr->Fit("f_walk_cor", "R");
 	
 
 
 	/*--------------------------------------------------*/
 	/*--------------------------------------------------*/
 
 	c10->cd(8);
 
 	Float_t time_walk_offset = f_walk_cor->GetParameter(0);
 	Float_t time_top = f_walk_cor->GetParameter(1);
 	Float_t time_bottom = f_walk_cor->GetParameter(2);
 
 
 	TString timewalk_str;
 	timewalk_str.Form( " (t0 - (%f + %f/(" + cor_str_low + " + %f)))", time_walk_offset, time_top, time_bottom);
 
 	t1->Draw( timewalk_str + ":" + cor_str_low + " >> timewalk_corr(70, 0, 70, 40, -20, 20)", "adc < 1300 && adc >= 500 && mycut", "colz");
 
 	TH2D* timewalk_corr = (TH2D*)gPad->GetPrimitive("timewalk_corr");
 
 
 	timewalk_corr->Draw("colz");
 
 	/*--------------------------------------------------*/
 	/*--------------------------------------------------*/
 
 //	c10->cd(8);
 
 	c10->cd(9);
 
 //	t1->Draw( timewalk_str + ":" + cor_str_low + " >> timewalk_corr(70, 0, 70, 40, -20, 20)", "adc < 1300 && adc >= 500 && mycut", "colz");
 
 	TH1D* y_proj = (TH1D*) timewalk_corr->ProjectionY("_py");
 	y_proj->DrawCopy("hist");
 
 
 
 	y_proj->GetXaxis()->SetRange(10, 30);
 //	y_proj->GetXaxis()->SetBit(TAxis::kAxisRange);
 
 	y_proj->Draw("hist");
 
 
 	TF1* time_res = new TF1("time_res", "gaus", -4, 4);
 
 	y_proj->Fit("time_res", "R");
 
 	time_res->Draw("same");
 
 	text1->DrawTextNDC( 0.65, 0.7, Form("Mean: %.3f", time_res->GetParameter(1)));
 	text1->DrawTextNDC( 0.65, 0.6, Form("Sigma: %.3f", time_res->GetParameter(2)));
 
 
 	c10->Update();

	c10->Print("time_walk/test.png");





	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/// ADC fitting for cross checking		

	TCanvas* c11 = new TCanvas("c11", "c11", 1200, 400);

	c11->Divide(3, 1);

	c11->cd(1);

//	t1->Draw("t0:adc >> timewalk_adc(900, 400, 1300, 40, 130, 170)", "adc < 1300 && adc >= 500 && mycut", "colz");
	t1->Draw("t0:adc >> timewalk_adc(900, 400, 1300, 70, 100, 170)", "adc < 1300 && adc >= 500", "colz");
	
	TH2D* timewalk_adc = (TH2D*)gPad->GetPrimitive("timewalk_adc");

	TF1* f_walk_adc = new TF1("f_walk_adc", "[0] + [1]/(x+[2])", 580, 900);

	f_walk_adc->SetParameter(0, 140);

	timewalk_adc->Fit("f_walk_adc", "R");


	/*--------------------------------------------------*/

	c11->cd(2);

	TString t0_adc_cor_str;	

	t0_adc_cor_str.Form("(t0 - (%f + %f/(adc + %f)))", f_walk_adc->GetParameter(0), f_walk_adc->GetParameter(1), f_walk_adc->GetParameter(2));

//	t1->Draw( t0_adc_cor_str + ":adc >> timewalk_adc_cor(900, 400, 1300, 40, -20, 20)", "adc < 1300 && adc >= 500 && mycut", "colz");

	t1->Draw( t0_adc_cor_str + ":adc >> timewalk_adc_cor(900, 400, 1300, 40, -20, 20)", "adc < 1300 && adc >= 500", "colz");


	TH2D* timewalk_adc_cor = (TH2D*)gPad->GetPrimitive("timewalk_adc_cor");

	TH1D* twalk_adc_y_proj = (TH1D*) timewalk_adc_cor->ProjectionY("cor_py");

	/*--------------------------------------------------*/

	c11->cd(3);

//	timewalk_adc_cor->Draw("colz");
	
	twalk_adc_y_proj->DrawCopy("hist");


// 
// 	twalk_adc_y_proj->GetXaxis()->SetRange(10, 30);
// //	y_proj->GetXaxis()->SetBit(TAxis::kAxisRange);
// 
// 	twalk_adc_y_proj->Draw("hist");
// 
// 
 	TF1* adc_time_res = new TF1("adc_time_res", "gaus", -4, 4);
// 
 	twalk_adc_y_proj->Fit("adc_time_res", "R");
// 
// 	adc_time_res->Draw("same");
// 
 	text1->DrawTextNDC( 0.65, 0.7, Form("Mean: %.3f", adc_time_res->GetParameter(1)));
 	text1->DrawTextNDC( 0.65, 0.6, Form("Sigma: %.3f", adc_time_res->GetParameter(2)));


	c11->Update();

	c11->Print("time_walk/test_1.png");

	delete c10;
	delete c11;


	TCanvas* c12 = new TCanvas("c12", "c12", 800, 400);

	c12->Divide(2, 1);

	c12->cd(1);

	t1->Draw("t0:adc >> timewalk_adc(900, 400, 1700, 70, 100, 170)", "", "colz");

	TH2D* timewalk_adc_raw = (TH2D*)gPad->GetPrimitive("timewalk_adc");

	TH1D* twalk_adc_raw_y_proj = (TH1D*) timewalk_adc_raw->ProjectionY("cor_py");

	c12->cd(2);

	twalk_adc_raw_y_proj->DrawCopy("hist");


 	TF1* adc_time_raw_f = new TF1("adc_time_raw", "gaus", twalk_adc_raw_y_proj->GetBinCenter(twalk_adc_raw_y_proj->GetMaximumBin())-2, twalk_adc_raw_y_proj->GetBinCenter(twalk_adc_raw_y_proj->GetMaximumBin())+2 );


 	twalk_adc_raw_y_proj->Fit("adc_time_raw", "R");


// 
 	adc_time_raw_f->Draw("same");
// 
 	text1->DrawTextNDC( 0.65, 0.7, Form("Mean: %.3f", adc_time_raw_f->GetParameter(1)));
 	text1->DrawTextNDC( 0.65, 0.6, Form("Sigma: %.3f", adc_time_raw_f->GetParameter(2)));

	c12->Print("time_walk/raw_t0_adc.png");


	/*--------------------------------------------------*/
	/// Slicing t0 in ADC
	
	// Fit_t0_adc_slice(t1);

	// Fit_dt_adc_slice(t1);

	// Fit_adc_dt_slice(t1);

    Fit_t0_dt_slice(t1);


}




/*--------------------------------------------------*/
/*--------------------------------------------------*/


void goodRICHEvent::Timewalk_Calabration_Raw(TTree* t1) {

	TCanvas* c10 = new TCanvas("c10", "c10", 1200, 1200);
	c10-> Divide(3,3);

	

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
//	t1->Draw("t0:adc >> time_corr(900, 400, 1300, 40, 130, 170)", "adc < 1300 && adc >= 500", "colz");
	
//	t1->Draw("t0:adc >> time_corr(900, 400, 1300, 40, 130, 170)", "adc < 1300 && adc >= 500 && mycut", "colz");


	c10->cd(4);


	t1->Draw("t0:delta_t >> time_walk_org(70, 0, 70, 40, 130, 170)", "adc < 1300 && adc >= 500", "colz");

//	t1->Draw("t0:adc >> time_corr_hi(900, 400, 1300, 40, 130, 170)", "adc < 1300 && adc >= 500", "colz");

	TH2D *time_walk_org = (TH2D*)gPad->GetPrimitive("time_walk_org");

	TF1* f_walk_org = new TF1("f_walk_org", "[0] + [1]/(x+[2])", 15, 60);
	time_walk_org->Fit("f_walk_org", "R");
	


	c10->cd(5);


	Float_t time_walk_offset_org = f_walk_org->GetParameter(0);
	Float_t time_top_org = f_walk_org->GetParameter(1);
	Float_t time_bottom_org = f_walk_org->GetParameter(2);

	TString timewalk_org_str;
	timewalk_org_str.Form( " (t0 - (%f + %f/(delta_t+%f)))", time_walk_offset_org, time_top_org, time_bottom_org);



	t1->Draw( timewalk_org_str + ": delta_t >> timewalk_before_org(70, 0, 70, 40, -20, 20)", "adc < 1300 && adc >= 500", "colz");

	TH2D* timewalk_before_org = (TH2D*)gPad->GetPrimitive("timewalk_before_org");

	c10->cd(6);


	TH1D* y_proj_org = (TH1D*) timewalk_before_org->ProjectionY("before_py");
	y_proj_org->DrawCopy("hist");

	y_proj_org->GetXaxis()->SetRange(10, 30);

	y_proj_org->Draw("hist");


 	TF1* time_res_org = new TF1("time_res_org", "gaus", -4, 5);
 
 	y_proj_org->Fit("time_res_org", "R");
 
 	time_res_org->Draw("same");
 
 	TText *text1 = new TText();
	text1->DrawTextNDC( 0.65, 0.7, Form("Mean: %.3f", time_res_org->GetParameter(1)));
	text1->DrawTextNDC( 0.65, 0.6, Form("Sigma: %.3f", time_res_org->GetParameter(2)));






   /*--------------------------------------------------*/
   /*--------------------------------------------------*/
   /*--------------------------------------------------*/



	c10->cd(7);


	TString cor_str_low;
	TString cor_str_hi;

	cor_str_low.Form( "((log(adc - %f)- %f) / %f)", y_offset, expo_off, expo_muti);
//	cor_str_low.Form( "delta_t");

//	t1->Draw("t0:" + cor_str_low + " >> time_corr_low(70, 0, 70, 40, 130, 170)", "adc < 550 && mycut", "colz");

	t1->Draw("t0:" + cor_str_low + " >> time_corr_low(70, 0, 70, 40, 130, 170)", "adc < 1300 && adc >= 500", "colz");



 	TH2D *time_walk_corr = (TH2D*)gPad->GetPrimitive("time_corr_low"); 
 	
 
 
 
 	TF1* f_walk_cor = new TF1("f_walk_cor", "[0] + [1]/(x+[2])", 10, 60);
 	time_walk_corr->Fit("f_walk_cor", "R");
 	
 
 
 	/*--------------------------------------------------*/
 	/*--------------------------------------------------*/
 
 	c10->cd(8);
 
 	Float_t time_walk_offset = f_walk_cor->GetParameter(0);
 	Float_t time_top = f_walk_cor->GetParameter(1);
 	Float_t time_bottom = f_walk_cor->GetParameter(2);
 
 
 	TString timewalk_str;
 	timewalk_str.Form( " (t0 - (%f + %f/(" + cor_str_low + " + %f)))", time_walk_offset, time_top, time_bottom);
 
 	t1->Draw( timewalk_str + ":" + cor_str_low + " >> timewalk_corr(70, 0, 70, 40, -20, 20)", "adc < 1300 && adc >= 500", "colz");
 
 	TH2D* timewalk_corr = (TH2D*)gPad->GetPrimitive("timewalk_corr");
 
 
 	timewalk_corr->Draw("colz");
 
 	/*--------------------------------------------------*/
 	/*--------------------------------------------------*/
 
 //	c10->cd(8);
 
 	c10->cd(9);
 
 //	t1->Draw( timewalk_str + ":" + cor_str_low + " >> timewalk_corr(70, 0, 70, 40, -20, 20)", "adc < 1300 && adc >= 500 && mycut", "colz");
 
 	TH1D* y_proj = (TH1D*) timewalk_corr->ProjectionY("_py");
 	y_proj->DrawCopy("hist");
 
 
 
 	y_proj->GetXaxis()->SetRange(10, 30);
 //	y_proj->GetXaxis()->SetBit(TAxis::kAxisRange);
 
 	y_proj->Draw("hist");
 
 
// 	TF1* time_res = new TF1("time_res", gauss_f, -4, 5);

// 	TF1* time_res = new TF1("time_res", "gaus(0)*expo(3)", -4, 5);

 	TF1* time_res = new TF1("time_res", "gaus", -4, 5);

 	y_proj->Fit("time_res", "R");
 
 	time_res->Draw("same");
 
 	text1->DrawTextNDC( 0.65, 0.7, Form("Mean: %.3f", time_res->GetParameter(1)));
 	text1->DrawTextNDC( 0.65, 0.6, Form("Sigma: %.3f", time_res->GetParameter(2)));
 
 
 	c10->Update();

	c10->Print("time_walk/test_raw.png");

	delete c10;

}



Double_t poissonf(Double_t*x,Double_t*par) {                                                                              
  return par[0]*TMath::Poisson(x[0] + par[2],par[1]);

}

/*--------------------------------------------------*/

Double_t gauss_n_expo(Double_t*x, Double_t*par) {

	Double_t xx = x[0];	

 	Double_t gaus_f = par[0]*exp(-0.5 * pow((xx-par[1])/par[2], 2));

 	Double_t expo_f = exp(par[3] - par[4]*(xx + par[5]));
	
  	return gaus_f + expo_f;
//  	return gaus_f ;

}
 

/*--------------------------------------------------*/
/// Slicing t0 in ADC

void Fit_t0_adc_slice(TTree* t1) {

	cout << "fuctuin" << endl;

	TGraph* g_mean = new TGraph();
	TGraph* g_RMS = new TGraph();

	TString cut_str;

	TCanvas *c1_temp = new TCanvas();
 
 	TText *text1 = new TText();
 

 	t1->Draw("t0:adc >> timewalk_slice(900, 400, 1700, 70, 100, 170)", cut_str, "colz");

 	c1_temp->Print("t0_adc_slicing/overall_t0_adc.png");
	

 	for (int i=500; i<1000; i=i+5) { 
 
 		cut_str.Form("adc > %i && adc <= %i", i, i+5 );
 
 		t1->Draw("t0:adc >> timewalk_slice(900, 400, 1700, 70, 100, 170)", cut_str, "colz");
 
  		TH2D* timewalk_adc_slice = (TH2D*)gPad->GetPrimitive("timewalk_slice");
  
 	 	TH1D* timewalk_adc_slice_y= (TH1D*) timewalk_adc_slice->ProjectionY("cor_py");
 
 // 		TF1* adc_time_slice_f = new TF1("adc_time_slice", "gaus", timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())-1, timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())+1);
 
 
 
  		TF1* adc_time_slice_f = new TF1("adc_time_slice", "gaus", timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())-1, timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())+1);
 
 
  		TF1* adc_time_slice_f_2 = new TF1("adc_time_slice_2", "gaus(0) + expo(3)", timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())-1, timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())+4);
 
 
  		TF1* adc_time_slice_f_3 = new TF1("adc_time_slice_3", "landau(0)", timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())-1, timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())+3);
 
 
  		TF1* adc_time_slice_f_4 = new TF1("adc_time_slice_4", gauss_n_expo, timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())-1, timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())+4, 6);
 
 
  		TF1* adc_time_slice_f_5 = new TF1("adc_time_slice_5", fit_combo, timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())-1, timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())+4, 6);
 
 
 
 
 
 
  // 		TF1* adc_time_slice_f_4 = new TF1("adc_time_slice_4", "[0]*exp(-((x-[1])*(x-[1])/[2]/[2])/2)", timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())-1, timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())+1);
  
  //		adc_time_slice_f_4->SetParLimits(0, 10, 100000);
  //		adc_time_slice_f_4->SetParameter(1, 113);
  //		adc_time_slice_f_4->SetParameter(2, 0.2);
  
  
  //		adc_time_slice_f_3->SetParLimits(2, 0.1, 1);
  		
  //	adc_time_slice_f_3->FixParameter(1, 0.3);
  
  
  		timewalk_adc_slice_y->DrawCopy("hist");
  
  //  		timewalk_adc_slice_y->Fit("adc_time_slice", "R");
  
  
    		timewalk_adc_slice_y->Fit("adc_time_slice", "R");
  
   		adc_time_slice_f_4->SetParameter(0, adc_time_slice_f->GetParameter(0));
   		adc_time_slice_f_4->SetParameter(1, adc_time_slice_f->GetParameter(1));
   		adc_time_slice_f_4->SetParameter(2, adc_time_slice_f->GetParameter(2));
  
  // 		adc_time_slice_f_4->SetParLimits(4, -1000000, 0);
  
   		adc_time_slice_f_4->SetParLimits(5, adc_time_slice_f->GetParameter(1), 20000);
  
  
    		timewalk_adc_slice_y->Fit("adc_time_slice", "R");
  
  
  // 		adc_time_slice_f_2->SetParameter(0, adc_time_slice_f->GetParameter(0));
  // 		adc_time_slice_f_2->SetParameter(1, adc_time_slice_f->GetParameter(1));
  // 		adc_time_slice_f_2->SetParameter(2, adc_time_slice_f->GetParameter(2));
  // 		adc_time_slice_f_2->SetParLimits(4, -10, 0);
  
  
  
  // 		timewalk_adc_slice_y->Fit("adc_time_slice_2", "R");
   
   		g_mean->SetPoint(g_mean->GetN(), i+2.5, adc_time_slice_f->GetParameter(1));
   		g_RMS->SetPoint(g_RMS->GetN(), i+2.5, adc_time_slice_f->GetParameter(2));
  
  //		delete adc_time_slice_f;
  
  
  
  
  // 
  //  		text1->DrawTextNDC( 0.65, 0.7, Form("Par 1: %.3f", adc_time_slice_f_3->GetParameter(0)));
  //  		text1->DrawTextNDC( 0.65, 0.6, Form("Par 2: %.3f", adc_time_slice_f_3->GetParameter(1)));
  //  		text1->DrawTextNDC( 0.65, 0.5, Form("Par 3: %.3f", adc_time_slice_f_3->GetParameter(2)));
  // 
  
  
   		text1->DrawTextNDC( 0.65, 0.7, Form("Par 1: %.3f", adc_time_slice_f->GetParameter(0)));
   		text1->DrawTextNDC( 0.65, 0.6, Form("Par 2: %.3f", adc_time_slice_f->GetParameter(1)));
   		text1->DrawTextNDC( 0.65, 0.5, Form("Par 3: %.3f", adc_time_slice_f->GetParameter(2)));
 
 
 
 
 // 		/*--------------------------------------------------*/
 // 		/*--------------------------------------------------*/
 // 		/// Fit: landau
 //
 
 
 //  		timewalk_adc_slice_y->DrawCopy("hist");
 // 
 // 
 // 		adc_time_slice_f_3->SetLineColor(4);
 // 
 // 
 //         adc_time_slice_f_3->SetParameter(0, timewalk_adc_slice_y->GetMaximum());
 //         adc_time_slice_f_3->SetParameter(1, timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin()) );
 //         adc_time_slice_f_3->SetParameter(2, 0.3);
 // 
 // 
 //    		timewalk_adc_slice_y->Fit("adc_time_slice_f_3", "R");
 // 		adc_time_slice_f_3->Draw("same");
 
 
 // 		/*--------------------------------------------------*/
 // 		/*--------------------------------------------------*/
 // 		/// Combined fit: gaussian + landau
 // 
 // 		
 // 		timewalk_adc_slice_y->DrawCopy("hist");
 // 
 //     	timewalk_adc_slice_y->Fit("adc_time_slice_5", "R");  
 //    
 //    
 //     	adc_time_slice_f_5->SetParameter(0,  1052  );
 //     	adc_time_slice_f_5->SetParameter(1,  110.8 );
 //     	adc_time_slice_f_5->SetParameter(2,  0.4   );
 //                                                    
 //     	adc_time_slice_f_5->SetParameter(3,  372   );
 //     	adc_time_slice_f_5->SetParameter(4,  112   );
 //     	adc_time_slice_f_5->SetParameter(5,  0.6   );
 
 
 
 
 //		if ( i == 800 )  {
 
 // 			TF1* adc_time_slice_f_5 = new TF1("adc_time_slice_5", fit_combo, timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())-1, timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())+4, 6);
 
 // 	    	adc_time_slice_f_5->SetParameter(0, 700  );
 // 	    	adc_time_slice_f_5->SetParameter(1, 110.5);
 // 	    	adc_time_slice_f_5->SetParameter(2, 0.4  );
 // 	    
 // 	    	adc_time_slice_f_5->SetParameter(3, 100);
 // 	    	adc_time_slice_f_5->SetParameter(4, 114 );
 // 	    	adc_time_slice_f_5->SetParameter(5, 1);
 
 //	    	timewalk_adc_slice_y->Fit("adc_time_slice_5", "R");  
 
 // 	        TF1 * f_combo = new TF1("f_combo", fit_combo, 108, 120, 6);
 
  	        TF1 * f_combo = new TF1("f_combo", fit_combo,  timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())-1,  timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())+6, 6);
             
  	        f_combo->SetParameter(0, timewalk_adc_slice_y->GetMaximum()  );
  	        f_combo->SetParameter(1, timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin()) );
  	        f_combo->SetParameter(2, 0.4   );
             
  	        f_combo->SetParameter(3, 372   );
  	        f_combo->SetParameter(4, timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin()) );
  	        f_combo->SetParameter(5, 0.6   );
 
 			timewalk_adc_slice_y->DrawCopy("hist");
 	
 	    	timewalk_adc_slice_y->Fit("f_combo", "R");  
 	   
 			f_combo->SetLineColor(4);
 
 			f_combo->Draw("same");
 
 
  
 	  		text1->DrawTextNDC( 0.65, 0.7, Form("Par 1: %.3f", f_combo->GetParameter(0)));
   			text1->DrawTextNDC( 0.65, 0.6, Form("Par 2: %.3f", f_combo->GetParameter(1)));
   			text1->DrawTextNDC( 0.65, 0.5, Form("Par 3: %.3f", f_combo->GetParameter(2)));
 
 
 
 


 // 			c1_temp->Print("fit_test_600.root");
 // 			c1_temp->Print("gauss_landau.png");
 
 //			exit(0);
 
 //		}
 
 
 
 		TString file_name;
 
 		file_name.Form("t0_adc_slicing/test_%i.png", i);
 
 		c1_temp->Print(file_name);
 		
 	}
 	
 	
 //		exit(0);
 
 	g_mean->Draw("A*");
 
 	TCanvas* c13 = new TCanvas("c13", "c13", 800, 400);
 
 	c13->Divide(2, 1);
 
 	c13->cd(1);
 
  	g_mean->Draw("A*");
 
 	c13->cd(2);
 
  	g_RMS->Draw("A*");
 
 	c13->Print("time_walk/sliced.png");


}




/*--------------------------------------------------*/
/// Slicing dt in ADC

void Fit_dt_adc_slice(TTree* t1) {

	cout << "fuctuin" << endl;

	TGraph* g_mean = new TGraph();
	TGraph* g_RMS = new TGraph();

	TString cut_str;

	TCanvas *c1_temp = new TCanvas();
 
 	TText *text1 = new TText();


 	t1->Draw("delta_t:adc", cut_str, "colz");

 	c1_temp->Print("dt_adc_slicing/overall_dt_adc.png");
	


 	for (int i=540; i<1000; i=i+5) { 
 
 		cut_str.Form("adc > %i && adc <= %i", i, i+5 );
 
 		t1->Draw("delta_t:adc >> timewalk_slice(900, 400, 1700, 100, 0, 100)", cut_str, "colz");
 
  		TH2D* timewalk_adc_slice = (TH2D*)gPad->GetPrimitive("timewalk_slice");

  	 	TH1D* timewalk_adc_slice_y= (TH1D*) timewalk_adc_slice->ProjectionY("cor_py");  

   		timewalk_adc_slice_y->DrawCopy("hist");
 
 		TString file_name;
 
 		file_name.Form("dt_adc_slicing/test_%i.png", i);
 
 		c1_temp->Print(file_name);
 		
 	}

 	g_mean->Draw("A*");
 
 	TCanvas* c13 = new TCanvas("c13", "c13", 800, 400);
 
 	c13->Divide(2, 1);
 
 	c13->cd(1);
 
  	g_mean->Draw("A*");
 
 	c13->cd(2);
 
  	g_RMS->Draw("A*");
 
 	c13->Print("time_walk/sliced.png");


}




/*--------------------------------------------------*/
/// Slicing ADC in dt

void Fit_adc_dt_slice(TTree* t1) {

	cout << "fuctuin" << endl;

	TGraph* g_mean = new TGraph();
	TGraph* g_RMS = new TGraph();

	TString cut_str;

	TCanvas *c1_temp = new TCanvas();
 
 	TText *text1 = new TText();

 	t1->Draw("delta_t:adc", cut_str, "colz");

 	c1_temp->Print("adc_dt_slicing/overall_dt_adc.png");



 	for (int i=0; i<90; i=i+1) { 
 
 		cut_str.Form("delta_t > %i && delta_t <= %i", i, i+1 );
 
 		t1->Draw("adc:delta_t >> timewalk_slice( 100, 0, 100, 900, 400, 1700)", cut_str, "colz");
 
  		TH2D* timewalk_adc_slice = (TH2D*)gPad->GetPrimitive("timewalk_slice");

  	 	TH1D* timewalk_adc_slice_y= (TH1D*) timewalk_adc_slice->ProjectionY("cor_py");  

   		timewalk_adc_slice_y->DrawCopy("hist");
 
 		TString file_name;
 
 		file_name.Form("adc_dt_slicing/test_%i.png", i);
 
 		c1_temp->Print(file_name);
 		
 	}

}


/*--------------------------------------------------*/
/// Slicing ADC in dt

void Fit_t0_dt_slice(TTree* t1) {

	cout << "fuction" << endl;

	TGraph* g_mean = new TGraph();
	TGraph* g_RMS = new TGraph();

	TString cut_str;

	TCanvas *c1_temp = new TCanvas();
 
 	TText *text1 = new TText();

 	t1->Draw("t0:delta_t >> timewalk_slice( 100, 0, 100, 70, 100, 170)", cut_str, "colz");

 	c1_temp->Print("t0_dt_slicing/overall_dt_adc.png");

 	for (int i=0; i<90; i=i+1) { 
 
 		cut_str.Form("delta_t > %i && delta_t <= %i", i, i+1 );
 
 		t1->Draw("t0:delta_t >> timewalk_slice( 100, 0, 100, 70, 100, 170)", cut_str, "colz");
 
  		TH2D* timewalk_adc_slice = (TH2D*)gPad->GetPrimitive("timewalk_slice");

  	 	TH1D* timewalk_adc_slice_y= (TH1D*) timewalk_adc_slice->ProjectionY("cor_py");  

   		timewalk_adc_slice_y->DrawCopy("hist");


  		TF1* adc_time_slice_f = new TF1("adc_time_slice", "gaus", timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())-1, timewalk_adc_slice_y->GetBinCenter(timewalk_adc_slice_y->GetMaximumBin())+1);



		timewalk_adc_slice_y->Fit("adc_time_slice");


   		text1->DrawTextNDC( 0.65, 0.7, Form("Par 1: %.3f", adc_time_slice_f->GetParameter(0)));
   		text1->DrawTextNDC( 0.65, 0.6, Form("Par 2: %.3f", adc_time_slice_f->GetParameter(1)));
   		text1->DrawTextNDC( 0.65, 0.5, Form("Par 3: %.3f", adc_time_slice_f->GetParameter(2)));


 
 		TString file_name;
 
 		file_name.Form("t0_dt_slicing/test_%i.png", i);
 
 		c1_temp->Print(file_name);
 		
 	}

}






















/*--------------------------------------------------*/

Double_t fit_gaus_fun(Double_t *x, Double_t*par) { 

	Double_t xx = x[0];

//	Double_t expo_res = exp(par[0] + par[1]*xx);  

	Double_t gaus_res = par[0] * exp( -0.5 * pow( (xx - par[1])/par[2], 2) );

	return gaus_res; 

}

/*--------------------------------------------------*/

Double_t fit_landau(Double_t *x, Double_t *par) {

	Double_t xx = x[0];

	Double_t res = par[0]*TMath::Landau(xx,par[1],par[2]);
	
	return res; 

}

/*--------------------------------------------------*/

Double_t fit_combo(Double_t *x, Double_t *par) {

      return fit_gaus_fun(x,par) + fit_landau(x,&par[3]);

}

/*--------------------------------------------------*/

Double_t fitf(Double_t *x, Double_t *par) {
 
    Double_t xx = 0;

    Double_t fitval =  par[0] * xx - par[1] + par[2]; 
 
    return fitval;

}


