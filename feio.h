#ifndef feio_hh
#define feio_hh

#include<stdint.h>
#include<vector>

namespace RICHfrontend{
 enum{	NCHANNELS = 192,
		NPIXELS = 64,
		MAXEDGES = 2000};

 uint8_t chan2pix[] = {60, 58, 59, 57, 52, 50, 51, 49, 44, 42, 43, 41, 36, 34, 35, 33, 28, 26, 27, 25, 20, 18, 19, 17, 12, 10, 11, 9, 4, 2, 3, 1, 5, 7, 6, 8, 13, 15, 14, 16, 21, 23, 22, 24, 29, 31, 30, 32, 37, 39, 38, 40, 45, 47, 46, 48, 53, 55, 54, 56, 61, 63, 62, 64};
 uint8_t pix2chan[] = {31, 29, 30, 28, 32, 34, 33, 35, 27, 25, 26, 24, 36, 38, 37, 39, 23, 21, 22, 20, 40, 42, 41, 43, 19, 17, 18, 16, 44, 46, 45, 47, 15, 13, 14, 12, 48, 50, 49, 51, 11, 9, 10, 8, 52, 54, 53, 55, 7, 5, 6, 4, 56, 58, 57, 59, 3, 1, 2, 0, 60, 62, 61, 63};


 struct rawEvent{
	Int_t trigID;
	ULong64_t timeStamp;
	UChar_t hold;

	UShort_t fadc[NCHANNELS];

	UShort_t fnedge;
	UShort_t ftdc[MAXEDGES];
	UShort_t fchan[MAXEDGES];
	Bool_t fpolar[MAXEDGES];
 };


 class RICHEvent{
  protected:
	bool fLeadingEdge;
	bool fTrailingEdge;

	ULong64_t timeStamp;
	Int_t trigID, hold;
	UShort_t fadc[NCHANNELS];

	std::vector<UShort_t> ftdc[NCHANNELS];
	std::vector<Bool_t> fpolar[NCHANNELS];

	std::vector<UShort_t> ftime[NCHANNELS];
	std::vector<UShort_t> fdur[NCHANNELS];

  public:
	RICHEvent():fLeadingEdge(true),fTrailingEdge(false){};

	void Clear(){
		trigID=0;
		timeStamp=0;

		for(int ich=0;ich<NCHANNELS;ich++){
			ftdc[ich].clear();
			fpolar[ich].clear();
			ftime[ich].clear();
			fdur[ich].clear();
		}
	};

	void Fill(rawEvent &rev){
		Clear();

		trigID = rev.trigID;
		timeStamp = rev.timeStamp;
		hold = rev.hold;

		std::copy(rev.fadc, rev.fadc+NCHANNELS, fadc);

		for(int iedge=0; iedge<rev.fnedge; iedge++){
			int ichan = rev.fchan[iedge];
			fpolar[ichan].push_back(rev.fpolar[iedge]);
			ftdc[ichan].push_back(rev.ftdc[iedge]);
		}
	};

	~RICHEvent(){};
 };
}

#endif
