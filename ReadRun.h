
#ifndef _ReadRun 
#define _ReadRun

// many of these includes are probably not needed anymore

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TComplex.h"
//root
#include <TLine.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TF1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <THStack.h>
#include <THistPainter.h>
#include <TText.h>
#include <TSpectrum.h>   // peakfinder
#include <TPolyMarker.h> // peakfinder
#include <TError.h>      // root verbosity level
#include <TSystem.h>     // root verbosity level
#include <TLatex.h>      // root verbosity level

//#include <sys/resource.h>
//C, C++
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <numeric>
#include <tuple>
#include <map>


class ReadRun/* : public TObject*/ {

private:

	//int maxnchannels;				// number of channels (32)
	TClonesArray* rundata;		// data array
	//TClonesArray* blh;

	double** amplValuessum;		// collects sums of all waveforms for each channel

	#pragma pack(1) // padding suppression
	// struct copied from
	// WaveCatcher binary -> root converter
	// by manu chauveau@cenbg.in2p3.fr
	struct event_data
	{
		int EventNumber;
		double EpochTime;
		unsigned int Year;
		unsigned int Month;
		unsigned int Day;
		unsigned int Hour;
		unsigned int Minute;
		unsigned int Second;
		unsigned int Millisecond;
		unsigned long long int TDCsamIndex;
		int nchannelstored;
	};

	struct channel_data_without_measurement
	{
		int channel;
		int EventIDsamIndex;
		int FirstCellToPlotsamIndex;
		short waveform[1024];
	};

	struct channel_data_with_measurement
	{
		int channel;
		int EventIDsamIndex;
		int FirstCellToPlotsamIndex;
		float	MeasuredBaseline;
		float	AmplitudeValue;
		float	ComputedCharge;
		float	RiseTimeInstant;
		float	FallTimeInstant;
		float	RawTriggerRate;
		short waveform[1024];
	};
	#pragma pack() // padding suppression

public:
	
	// plots amplValuessum
	void PlotChannelSums(bool = true);
	
	// baseline correction (shifts all waveforms individually)
	void CorrectBaseline(float, float = -999);
	void CorrectBaselineMinSlopeRMS(int = 100, bool = false, double = 10, bool = false, int = 0);

	// average all waveforms to simplify peak ID
	void SmoothAll(double = 5);
	
	// functions for charge spectrum
	int* GetIntWindow(TH1F*, float, float, float, float, int);
	void PrintChargeSpectrumWF(float, float, float = 0, float = 300, int = 1, float = 0., float = 0.);
	TH1F* ChargeSpectrum(int, float, float, float = 0, float = 300, float = -50, float = 600, int = 750);
	void PrintChargeSpectrum(float, float, float = 0, float = 300, float = -50, float = 600, int = 750, float = 0., float = 0.);
	

	// helper functions
	stringstream list_files(const char*, const char*);	// find data files
	TH1F* Getwf(int, int, int = 1);						// channel, eventnr, color
	double* getx();										// x values
	double* gety(int, int);								// y values for waveform(ch, event)
	double* gety(TH1F*);								// y values for histogram

	void Convolute(double*&, double*, double*, int, int);	// convolution for filtering waveforms
	void SmoothArray(double*&, int, double = 1., bool = false);	// filtering


	ReadRun(string);

	virtual ~ReadRun();

	//int nbinsdata;
	int nevents;				// number of triggered events
	int nchannels;
	int nwf;					// number of waveforms (nchannels*nacquisitions)

	float SP;					// ns per bin
	float pe;					// mV*ns ????
	double coef;				// ?????
	int binNumber;				// 1024 samples per waveform

	int* maxSumBin;				// For fixed intergration window (triggered acqusition)

	vector<vector<float>> baseline_correction_result;

	ClassDef(ReadRun, 1)
};




#endif



class Fitf {
public:
	// use constructor to customize your function object

	double operator() (double* x, double* p) {
		//0 - N0: Normalization (~Number of events)
		//1 - mu: for generalized poisson distribution
		//2 - lambda: Borel-branching parameter for prompt crosstalk probability 1-exp(-lambda)

		//3,4 -sogma0, sigma1
		//5 - G: gain
		//6 - B: Pedestal
		double sum = 0;
		for (int kint = 0; kint <= 15; kint++) {
			double mu = p[1];
			double lambda = p[2];

			double sigma0 = p[3];
			double sigma1 = p[4];
			
			double G = p[5];
			double B = p[6];

			double k = static_cast<double>(kint);
			double sigmaK = sqrt(sigma0 * sigma0 + k * sigma1 * sigma1);
			//generalized poisson envelope
			double gp = mu * TMath::Power((mu + k * lambda), k - 1) * TMath::Exp(-(mu + k * lambda)) / TMath::Factorial(kint); 

			sum += p[0] * gp * (1. / sqrt(2. * TMath::Pi()) / sigmaK) * TMath::Exp(-TMath::Power(((x[0] - (k * G + B)) / sqrt(2) / sigmaK), 2));
		}
		return sum;
	};
};

class Fitf_biased {
public:
	// use constructor to customize your function object

	double operator() (double* x, double* p) {
		//0 - N0: Normalization (~Number of events)
		//1 - mu: for generalized poisson distribution
		//2 - lambda: Borel-branching parameter for prompt crosstalk probability 1-exp(-lambda)

		//3,4 -sogma0, sigma1
		//5 - G: gain
		//6 - B: Pedestal
		double sum = 0;
		for (int kint = 0; kint <= 15; kint++) {
			double mu = p[1];
			double lambda = p[2];

			double sigma0 = p[3];
			double sigma1 = p[4];

			double G = p[5];
			double B = p[6];

			double a_ped = p[7];
			double x_ped = p[8];

			double k = static_cast<double>(kint);
			double sigmaK = sqrt(sigma0 * sigma0 + k * sigma1 * sigma1);
			//generalized poisson envelope
			double gp = mu * TMath::Power((mu + k * lambda), k - 1) * TMath::Exp(-(mu + k * lambda)) / TMath::Factorial(kint);

			if (kint == 0) {
				sum += a_ped * p[0] * gp * (1. / sqrt(2. * TMath::Pi()) / sigmaK) * TMath::Exp(-TMath::Power(((x[0] - x_ped) / sqrt(2) / sigmaK), 2));
			}
			else {
				sum += p[0] * gp * (1. / sqrt(2. * TMath::Pi()) / sigmaK) * TMath::Exp(-TMath::Power(((x[0] - (k * G + B)) / sqrt(2) / sigmaK), 2));
			}
		}
		return sum;
	};
};