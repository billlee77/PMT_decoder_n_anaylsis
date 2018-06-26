#include<iostream>
#include<TFile.h>
#include<TTree.h>
#include"feio.h"

using namespace RICHfrontend;

int main(int argc, char** argv)
{
 TString fname;
 TString dirname(".");
 for(int iarg=1;iarg<argc;iarg++){
	TString args(argv[iarg]);
	if(args=="-d"){
		iarg++;
		dirname = argv[iarg];
	}
	else fname = argv[iarg];
 }

 rawEvent rawEv;

 FILE *lbinaryFile = fopen(fname.Data(), "r");
 if(lbinaryFile==NULL){
	std::cerr<<fname.Data()<<" can not be opened ---------------------"<<std::endl;
	return 111;
 }

 TFile* ff = new TFile(Form("%s/%s.root", dirname.Data(), basename(fname.Data())), "recreate", "", 5);
 TTree* tt = new TTree("h22", "RICH FE tree");
 tt->Branch("trigID", &rawEv.trigID, "trigID/I");
 tt->Branch("timeStamp", &rawEv.timeStamp, "timeStamp/l");
 tt->Branch("hold", &rawEv.hold, "hold/b");

 tt->Branch("fadc", rawEv.fadc, Form("fadc[%d]/s", NCHANNELS));

 tt->Branch("fnedge", &rawEv.fnedge, "fnedge/s");
 tt->Branch("ftdc", rawEv.ftdc, "ftdc[fnedge]/s");
 tt->Branch("fchan", rawEv.fchan, "fchan[fnedge]/s");
 tt->Branch("fpolar", rawEv.fpolar, "fpolar[fnedge]/O");

 int tag_idx=0, tag=15;
 unsigned int val, asic;
 int nadc = 0;

 while(true){
	int nbytes = fread(&val, 1, sizeof(val), lbinaryFile);
	if(nbytes<=0) break;

	if (val & 0x80000000) {
		tag = (val >> 27) & 0xF;
		tag_idx = 0;
	} else {
		tag_idx++;
	}

	switch (tag) {
		case 0:
			break;

		case 1:
			break;

		case 2: // Event Header

			if(nadc==192)
				tt->Fill();
			nadc=0;

			rawEv.fnedge=0;
			std::fill(rawEv.fadc, rawEv.fadc+NCHANNELS, 0);

			rawEv.trigID = (val >> 0) & 0x3FFFFF;

			break;

		case 3: // Time stamp
			if (tag_idx == 0)
				rawEv.timeStamp = (val >> 0) & 0xFFFFFF;
			else if (tag_idx == 1) {
				unsigned long int timeStampH = (val >> 0) & 0xFFFFFF;
				rawEv.timeStamp += timeStampH<<24;
				rawEv.timeStamp <<= 3;
			}
			else{
				std::cerr<<"tag_idx>1 for trigger time"<<std::endl;
				return 666;
			}
			break;

		case 8: // TDC data
			rawEv.fpolar[rawEv.fnedge] = (val >> 26) & 0x1; // 0 = rising, 1 = falling seen by FPGA
			rawEv.fchan[rawEv.fnedge] = (val >> 16) & 0xFF;
			rawEv.ftdc[rawEv.fnedge] = (val >> 0) & 0xFFFF;
			rawEv.fnedge++;

			break;

		case 9: // ADC data
			if (tag_idx == 0) {
				asic = (val >> 0) & 0x03; // ASIC bit 1 and 0
				rawEv.hold = (val >> 8) & 0xFF; // HOLD1 bit[15..8]
				rawEv.hold = (val >> 16) & 0xFF; // HOLD2 bit[23..16]
			} else {
				int chLower = 2 * (tag_idx - 1) + 0;
				int chUpper = 2 * (tag_idx - 1) + 1;
				rawEv.fadc[asic*64 + chLower] = (val >> 0) & 0xFFF; // ADC [11.. 0]
				rawEv.fadc[asic*64 + chUpper] = (val >> 16) & 0xFFF; // ADC [27..16]

				nadc+=2;
			}
			break;

		default:
			printf("Junk data...Exit\n");
			fclose(lbinaryFile);
			tt->Write();
			ff->Close();
			return -1;
	}
 }
 fclose(lbinaryFile);

 tt->Write();
 ff->Close();

 return 0;
}

