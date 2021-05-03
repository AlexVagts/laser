#include <iostream>
#include <cmath> 
#include <string.h>
#include <stdio.h>  
#include <math.h> 
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>

#include "ReadRun.cc"
using namespace std;


void rtst1() // main
{
	int which = 18; //select meas

	// better create seperate rtst.c just for DC measurements to reduce if statements and improve readability and usability
	bool isDC = false;

	string path;

	switch (which) { //specify folder to run below, ALL bin files in this folder will be used
	case(0): {
		path = "C:/SHiP/data/76_calib_vb40_tune7840/"; // all SiPMs at once
		break;
	}
	case(1): {
		path = "C:/SHiP/data/70_dc_vb40/"; // all SiPMs at once
		isDC = true;
		break;
	}
	case(2): {
		path = "C:/SHiP/data/88_calib_vb40_tune7840/"; // one SiPM
		//isDC = true;
		break;
	}
	case(3): {
		path = "C:/SHiP/data/82_DC_vb40/"; // one SiPM
		isDC = true;
		break;
	}
	
	case(4): {
		path = "C:/SHiP/data/103_calib_vb40_tune8230_nopz/"; // one SiPM
		//isDC = true;
		break;
	}
	case(5): {
		path = "C:/SHiP/data/100_DC_vb40_nopz/"; // one SiPM
		isDC = true;
		break;
	}

	case(6): {
		path = "C:/SHiP/data/106_calib_vb40_tune8230_pz520/"; // one SiPM
		//isDC = true;
		break;
	}
	case(7): {
		path = "C:/SHiP/data/107_calib_vb40_tune8230_pz700/"; // one SiPM
		//isDC = true;
		break;
	}

	case(8): {
		path = "C:/SHiP/data/140_calib_vb44_tune8280_switch3/"; // one SiPM
		//isDC = true;
		break;
	}
	case(9): {
		path = "C:/SHiP/data/136_calib_vb40_tune8030_switch3/"; // one SiPM
		//isDC = true;
		break;
	}
	case(10): {
		path = "C:/SHiP/data/135_DC_vb44_switch3/"; // one SiPM
		isDC = true;
		break;
	}
	case(11): {
		path = "C:/SHiP/data/131_DC_vb40_switch3/"; // one SiPM 
		isDC = true;
		break;
	}
	case(12): {
		path = "C:/SHiP/data/161_calib_vb44_tune8280_switch5/"; // one SiPM 
		break; //
	}
	case(13): {
		path = "/home/alex/Dokumente/Studium/RootReader/data/151_calib_vb44_tune8280_switch4/"; // one SiPM 
		cout << "Der Pfad ist " << path << endl;
		break; //
	}
	case(14): {
		path = "C:/SHiP/data/149_calib_vb42_tune8220_switch4/"; // one SiPM 
		break; //
	}
	case(15): {
		path = "C:/SHiP/data/147_calib_vb40_tune8030_switch4/"; // one SiPM 
		break; //104_calib_vb40_tune8230_pz331
	}

	case(16): {
		path = "/home/alex/Dokumente/Studium/RootReader/data/104_calib_vb40_tune8230_pz331/"; // one SiPM 
		break; //
	}
	case(17): {
		path = "/home/alex/Dokumente/Studium/RootReader/data/103_calib_vb40_tune8230_nopz/"; // one SiPM 
		break; //
	}

	case(18): {
		path = "/home/alex/Dokumente/Studium/RootReader/data/196_calib_vb41_tune8180_switch45_nopz/"; // one SiPM 
		break; //
	}


	default: {
		path = "C:/SHiP/data/3_calib_vb56_tune8350/"; // ???
		break;
	}
	}

	bool saveplots = 0;
	
	ReadRun mymeas(path);

	//apply baseline correction to ALL waveforms
	//mymeas.CorrectBaseline(60.);	// <- slow when not compiled
	mymeas.CorrectBaselineMinSlopeRMS(80, true, 5, false, 300);	// <- VERY slow 

	////plotting

	//investigate individual waveforms
	//TCanvas* tstc = new TCanvas("tstc", "", 1600, 1000);
	//TH1F* histo = mymeas.Getwf(0, 0);			//
	//TH1F* histo2 = mymeas.Getwf(1, 186, 2);		//
	//TH1F* histo3 = mymeas.Getwf(3, 170, 3);	//
	//histo->Draw();
	//histo2->Draw("same");
	//histo3->Draw("same");
	//tstc->BuildLegend(0.85, 0.70, .99, .95);
	//if (saveplots) tstc->SaveAs("plots.pdf");
	
	// sums of all events per channel
	mymeas.PlotChannelSums(false);


	// investigate charge spectrum. should see photo electron peaks here
	float intwindowminus = 1.;	// assume avalanche signal starts 8 ns before maximum ...
	float intwindowplus = 1.;	// ... and end 30 ns afte maximum 
	float findmaxfrom = 110.;// 105.;	// assume signal from laser arrives between here ...
	float findmaxto = 127.;		// ... and here (150 - 220 ns, depends on trigger delay setting)	
	
	if (isDC) {
		findmaxfrom = 10 + intwindowminus;		// set = 5 + intwindowminus for dark count measurement
		findmaxto = 280. - intwindowplus;	// set = 300 - intwindowplus for dark count measurement
	}

	// plot all channels
	if (isDC) {
		mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -5, 100, 100);//-20, 180, 200);
	}
	else {
		mymeas.PrintChargeSpectrum(intwindowminus, intwindowplus, findmaxfrom, findmaxto, -10, 100, 200/*, 10, 100*/);
	}
	// plot waveforms for certain events with integration window
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, 170, -2., 10.);
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, 5100, -2, 10);
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, 15, -2, 10);
	mymeas.PrintChargeSpectrumWF(intwindowminus, intwindowplus, findmaxfrom, findmaxto, 4001, -2, 10);
}
