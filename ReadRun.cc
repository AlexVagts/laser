#include "ReadRun.h"

ClassImp(ReadRun)

ReadRun::ReadRun(string path) {
	// Wavecatcher hardware/software properties
	SP = 0.3125;					// ns per bin
	pe = 47.46;					//mV*ns ????
	coef = 2.5 / (4096 * 10);	//?????
	binNumber = 1024;				//default: 1024, hardcoded later on so it can be removed
	const int nChannelsWC = 64;			//max number of channels default: 32


	rundata = new TClonesArray("TH1F", 1e7); //raw data will be stored here as TH1F
	rundata->BypassStreamer(kFALSE);
	TClonesArray& testrundata = *rundata;

	//blh = new TClonesArray("TH1F", 1e7); 
	//blh->BypassStreamer(kFALSE);
	//TClonesArray& testblh = *blh;

	// verbosity
	bool debug_header = 0;
	bool debug_data = 0;


	unsigned short output_channel;
	unsigned int output_event;
	unsigned long long int output_tdc;
	unsigned short output_nbchannels;

	amplValuessum = new double* [nChannelsWC]; //sum of all wf for each channel
	for (int i = 0; i < nChannelsWC; i++) {//init
		amplValuessum[i] = new double[binNumber];
		for (int k = 0; k < binNumber; k++) amplValuessum[i][k] = 0.;
	}

	maxSumBin = new int[nChannelsWC];

	//Start reading the raw data from .bin files.
	stringstream inFileList = list_files(path.c_str(), ".bin"); //all *.bin* files in folder path
	int nitem = 1;
	string fileName;
	int fileCounter = 0;
	int currentPrint = -1;
	int wfcounter = 0;

	while (inFileList >> fileName) {
		// file loop

		fileName = path + fileName;
		ifstream input_file(fileName.c_str(), std::ios::binary | std::ios::in);

		bool has_measurement = false;

		if (!input_file.is_open()) {
			printf("*** failed to open '%s'\n", fileName.c_str());
			continue;
		}
		printf("+++ reading '%s' ...\n", fileName.c_str());

		// Header
		string header_line;
		// HEADER 1 //
		//
		// "=== DATA FILE SAVED WITH SOFTWARE VERSION: V?.??.? ==="
		//
		getline(input_file, header_line, '\n');
		if (debug_header) printf("%s\n", header_line.data());
		assert(header_line[0] == '=');

		size_t header_version_first = header_line.find_last_of('V');
		size_t header_version_last = header_line.find_first_of(' ', header_version_first);
		string software_version = header_line.substr(header_version_first, header_version_last - header_version_first);
		if (debug_header) printf("    |- data version = '%s'\n", software_version.data());

		//if (software_version == "V2.9.13")
		//	;
		//else if (software_version == "V2.9.15")
		//	;
		//else if (debug_header) printf("*** unsupported data version\n");


		// HEADER 2 // 
		// "=== WAVECATCHER SYSTEM OF TYPE ?? WITH ?? CHANNELS AND GAIN: ??? ==="
		getline(input_file, header_line, '\n');
		if (debug_header) printf("%s\n", header_line.data());
		assert(header_line[0] == '=');

		// HEADER 3 //
		// === Rate coincidence masks ... === Posttrig in ns for SamBlock ... ===
		getline(input_file, header_line, '\n');
		if (debug_header) printf("%s\n", header_line.data());
		assert(header_line[0] == '=');

		// HEADER 4 // 
		// V2.9.13: === DATA SAMPLES [1024] in Volts == NB OF CHANNELS ACQUIRED: 64 == Sampling Period: 312.5 ps == INL Correction: 1
		// V2.9.15: === DATA SAMPLES [1024] in Volts == NB OF CHANNELS ACQUIRED: 64 == Sampling Period: 312.5 ps == INL Correction: 1 == MEASUREMENTS: 0 ===
		getline(input_file, header_line, '\n');
		if (debug_header) printf("%s\n", header_line.data());
		assert(header_line[0] == '=');

		size_t nsamples_first = 1 + header_line.find_last_of('[');
		size_t nsamples_last = header_line.find_first_of(']', nsamples_first);
		string nsamples_str = header_line.substr(nsamples_first, nsamples_last - nsamples_first);

		int nsamples = atoi(nsamples_str.data());
		if (debug_header) printf("    |- data sample  = %d\n", nsamples);

		size_t nchannels_first = 10 + header_line.find("ACQUIRED: ", nsamples_first);
		size_t nchannels_last = header_line.find_first_of(' ', nchannels_first);
		string nchannels_str = header_line.substr(nchannels_first, nchannels_last - nchannels_first);

		nchannels = atoi(nchannels_str.data());
		if (debug_header) printf("    |- nchannels    = %d\n", nchannels);

		if (software_version == "V2.9.15" || software_version == "V2.9.16" || software_version == "V2.10.1") {
			size_t has_measurement_first = 14 + header_line.find("MEASUREMENTS: ", nsamples_first);
			size_t has_measurement_last = header_line.find_first_of(' ', has_measurement_first);
			string has_measurement_str = header_line.substr(has_measurement_first, has_measurement_last - has_measurement_first);
			has_measurement = atoi(has_measurement_str.data());
		} 
		else {
		//if (software_version == "V2.9.13") {
			// V2.9.13 has always measurement stored
			// (everything is set to 0 when disabled!)
			has_measurement = true;
		}

		if (debug_header) printf("    `- measurement  = %d\n", has_measurement);

		// end of header reader


		event_data an_event;

		while (input_file.read((char*)(&an_event), sizeof(an_event))) {
			//file loop

			if (debug_data) printf("%03d has %d channels\n", an_event.EventNumber, an_event.nchannelstored);

			output_event = an_event.EventNumber;
			output_tdc = an_event.TDCsamIndex;
			output_nbchannels = an_event.nchannelstored;

			if (debug_data && output_event % 200 == 0) printf("EventNr: %d, nCh: %d\n", output_event, output_nbchannels);
			
			for (int ch = 0; ch < nchannels ; ++ch) {
				// 
				channel_data_with_measurement a_channel_data;

				if (has_measurement) {
					// read with 'channel_data_with_measurement' struct
					input_file.read((char*)(&a_channel_data), sizeof(channel_data_with_measurement));
				}
				else {
					// read with 'channel_data_without_measurement' struct
					channel_data_without_measurement a_channel_data_without_measurement;
					input_file.read((char*)(&a_channel_data_without_measurement), sizeof(channel_data_without_measurement));

					// copy the content into 'channel_data_with_measurement' struct
					a_channel_data.channel = a_channel_data_without_measurement.channel;
					a_channel_data.EventIDsamIndex = a_channel_data_without_measurement.EventIDsamIndex;
					a_channel_data.FirstCellToPlotsamIndex = a_channel_data_without_measurement.FirstCellToPlotsamIndex;
					memcpy(a_channel_data.waveform, a_channel_data_without_measurement.waveform, 1024 * sizeof(short));
				}


				output_channel = a_channel_data.channel;
				if (debug_data) printf("- reading channel %d\n", output_channel);

				TString name(Form("channel_%02d, event %05d ", ch, an_event.EventNumber));
				TString title(Form("Channel %d, event %d data", ch, an_event.EventNumber));
				TH1F* hCh = (TH1F*)testrundata.ConstructedAt(wfcounter);
				hCh->SetName(name.Data());
				hCh->SetTitle(title.Data());
				hCh->SetBins(binNumber, -0.5 * SP, 1023.5 * SP);

				float val = 0.;
				for (int s = 0; s < binNumber; ++s) {
					val = a_channel_data.waveform[s] * coef * 1000.;
					hCh->SetBinContent(s + 1, val);
					//hCh->SetBinError(s, 0.5); //The error of each value in each bin is set to 0.5 mV -> Why??

					// channel sums
					amplValuessum[ch][s] += static_cast<double>(val);
				}

				//hCh->SetLineColor(ch + 1); // gets a bit too colorful
				//hCh->SetMarkerColor(ch + 1);

				wfcounter++;
			} // for ch

		} // while an_event

		input_file.close();
	} // for file_id


	// get bins where the sum spectrum has its maximum for runs with fixed trigger delay and fixed integration window relative to the max of the sum spectrum (not working for DC measurement)
	for (int ch = 0; ch < nchannels; ch++) {
		double max = 0.;
		for (int i = 0; i < binNumber; i++) {
			if (amplValuessum[ch][i] > max) {
				max = amplValuessum[ch][i];
				maxSumBin[ch] = i;
			}
		}
	}

	nevents = output_event;
	nwf = wfcounter;
}


ReadRun::~ReadRun() {
	// Destructor 
	rundata->Clear();
	//delete[] maxSumBin;
	//delete baseline_correction_result;
	cout << "deleting nothing currently..." << endl;
}




// plot sums of all waveforms for each channel



void ReadRun::PlotChannelSums(bool doaverage) {
	// doaverage: if true it will plot the running average +/- 4 bins

	double* xv = getx();
	TMultiGraph* mgsums = new TMultiGraph();
	mgsums->SetTitle("channel sums; t [ns]; amplitude [arb.]");

	for (int i = 0; i < nchannels; i++) {
		double* yv = amplValuessum[i];
		if (doaverage) SmoothArray(yv, binNumber, 4);
		TGraph* gr = new TGraph(binNumber, xv, yv);
		delete[] yv;

		TString name(Form("channel_%02d", i));
		TString title(Form("Channel %d", i));
		gr->SetName(name.Data());
		gr->SetTitle(title.Data());
		gr->SetLineColor(i + 1);
		gr->SetMarkerColor(i + 1);
		mgsums->Add(gr);

		cout << "\n at " << i;
	}
	delete[] xv;

	TCanvas* sumc = new TCanvas("Sums", "", 1600, 1000);
	mgsums->Draw("APL");
	mgsums->GetYaxis()->SetRangeUser(-1e4, 10e4);
	sumc->BuildLegend(0.85, 0.70, .99, .95);
	sumc->SaveAs("channelsums.pdf");
}



// averaging all waveforms (for testing)



void ReadRun::SmoothAll(double sigma) { //deprecated since it can be done with baseline correction??
	// just for testing, not very efficient
	cout << "\nsmoothing wfs";
	for (int j = 0; j < nwf; j++) {
		TH1F* his = ((TH1F*)rundata->At(j));
		double* yvals = gety(his);
		SmoothArray(yvals, binNumber, sigma);
		for (int i = 1; i < his->GetNbinsX(); i++) his->SetBinContent(i, yvals[i]);
		delete[] yvals;
		if (j % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j) / static_cast<float>(nwf) << "% -" << flush;
	}
}



// baseline correction



void ReadRun::CorrectBaseline(float tCut, float tCutEnd) {
	// corrects the baseline (DC offset) of all waveforms
	// tCut: time denoting the end of integration window (or the beginning if tCutEnd is set)

	int iCut, iCutEnd;
	float corr = 0;

	for (int j = 0; j < nwf; j++) {
		if (j == 0)  printf("Baseline correction (%d waveforms) :: ", nwf);
		TH1F* his = ((TH1F*)rundata->At(j));
		iCut = his->GetXaxis()->FindBin(tCut);

		if (tCutEnd <= 0) { //
			corr = his->Integral(1, iCut) / (iCut - 1);
		}
		else {
			iCutEnd = his->GetXaxis()->FindBin(tCutEnd);
			corr = his->Integral(iCut, iCutEnd) / (iCutEnd - iCut);
		}

		if (j % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j) / static_cast<float>(nwf) << "% -" << flush;
		
		// write corrected values to histograms
		for (int i = 1; i < his->GetNbinsX(); i++) his->SetBinContent(i, his->GetBinContent(i) - corr);
	}
}


void ReadRun::CorrectBaselineMinSlopeRMS(int nIntegrationWindow, bool doaverage, double sigma, bool convolution, int max_bin_for_baseline) {
	// corrects the baseline (DC offset) of all waveforms
	// determines the region of nIntegrationWindow bins where the root mean square of the slope of the smoothed waveform reaches its minimum
	// experimental, use with care and check output, 
	// it is very slow! only needed for DC measurement of for SiPMs with 
	// nIntegrationWindow: number of bins of integration window (should be ~10 % of spectrum) 

	float corr = 0;
	int binNumberSlope = binNumber - 1;
	double* slope = new double[binNumberSlope];

	int min_distance_from_max = 25 + nIntegrationWindow;
	

	for (int j = 0; j < nwf; j++) {
		if (j == 0)  printf("Baseline correction (%d waveforms) :: ", nwf);
		TH1F* his = ((TH1F*)rundata->At(j));
		double* yvals = gety(his); //find faster way
		SmoothArray(yvals, binNumber, sigma, convolution); // smoothing important to suppress variations in slope due to noise so the method is more sensitve to excluding peaks
		
		//calculate slope
		for (int i = 0; i < binNumberSlope; i++) slope[i] = yvals[i + 1] - yvals[i];
		
		float minchange = 1.e9;
		int iintwindowstart = 0;


		//find region where rms of the slope reaches minimum 
		int imax = his->GetMaximumBin();	// exclude 25 bins (~8 ns) before 
											// and NOT 100 bins (~32 ns) after maximum since it might catch ringing
		int search_before_max = imax - min_distance_from_max;

		float sum = 0.;
		float sqsum = 0.;
		for (int i = 0; i < binNumber - nIntegrationWindow - 1; i += 3) { // currently in steps of 3 bins (~1 ns) to make it faster
			if ((max_bin_for_baseline != 0 && i <= max_bin_for_baseline) || (max_bin_for_baseline == 0 && i <= search_before_max)) { // 
				sum = 0.;
				sqsum = 0.;
				for (int k = i; k < nIntegrationWindow + i; k++) {
					sum += slope[k];
					sqsum += (slope[k] * slope[k]);
				}
				if (sqsum + sum * sum < minchange) {
					minchange = sqsum + sum * sum;
					iintwindowstart = i;
				}
			}
		}

		corr = 0.;
		if (!doaverage) {
			corr = his->Integral(iintwindowstart, iintwindowstart + nIntegrationWindow) / static_cast<float>(nIntegrationWindow);
		}
		else {
			for (int i = iintwindowstart; i < iintwindowstart + nIntegrationWindow; i++) corr += yvals[i];
			corr /= static_cast<float>(nIntegrationWindow);
		}

		baseline_correction_result.push_back(vector<float>());
		baseline_correction_result[j].push_back(corr);
		baseline_correction_result[j].push_back(minchange);
		baseline_correction_result[j].push_back(static_cast<float>(iintwindowstart) * SP);
		baseline_correction_result[j].push_back(static_cast<float>(iintwindowstart + nIntegrationWindow) * SP);

		if (j % (nwf / 10) == 0) cout << " " << 100. * static_cast<float>(j) / static_cast<float>(nwf) << "% -" << flush;
		//if (j % 10000 == 0) printf("\nwf: %d, start: %d, min: %f, baseline correction: %f\n", j, iintwindowstart, minchange, corr);

		for (int i = 0; i < binNumber; i++) {
			if (!doaverage) his->SetBinContent(i, his->GetBinContent(i) - corr);
			else his->SetBinContent(i, yvals[i] - corr);
		}
		delete[] yvals; //delete slow
	}
	delete[] slope;
}



// functions for charge spectrum



int* ReadRun::GetIntWindow(TH1F* his, float windowlow, float windowhi, float start, float end, int channel) {
	// find maximum in range (start, end) and return bin numbers for [0] the max, [1] t_max - windowlow, and [2] t_max + windowhi
	// if (start < 0 || end < 0) doesn't return max and integration window is fixed t(max(sum_spectrum[channel])) +/- hi/lo
	// if (windowlow == start && windowwhi == end) doesn't return max and sets fixed integration window from start until end for all channels

	int istart, iend;
	int* foundindices = new int[3];// 
	foundindices[0] = 0;
	
	if (start < 0 || end < 0) {								// fixed integration window relative to maximum of sum spectrum for each channel
		foundindices[1] = his->GetXaxis()->FindBin(static_cast<float>(maxSumBin[channel]) * SP - windowlow);
		foundindices[2] = his->GetXaxis()->FindBin(static_cast<float>(maxSumBin[channel]) * SP + windowhi);
	}
	else if (windowlow == start && windowhi == end) {				// fixed integration window for all channels
		foundindices[1] = his->GetXaxis()->FindBin(windowlow);
		foundindices[2] = his->GetXaxis()->FindBin(windowhi);
	} 
	else {															// fixed integration window relative to maximum of each individual waveform
		istart = his->GetXaxis()->FindBin(start);
		iend = his->GetXaxis()->FindBin(end);

		if (istart<1 || iend>his->GetNbinsX()) {
			cout << "\nError: start or end out of range" << endl;
			return 0;
		}

		float max = -1e3;
		float val = 0;
		for (int i = istart; i < iend; i++) {
			val = his->GetBinContent(i);
			if (val > max) {
				max = val;
				foundindices[0] = i;
			}
		}

		foundindices[1] = his->GetXaxis()->FindBin(static_cast<float>(foundindices[0]) * SP - windowlow);
		foundindices[2] = his->GetXaxis()->FindBin(static_cast<float>(foundindices[0]) * SP + windowhi);
	}
	return foundindices;
}


void ReadRun::PrintChargeSpectrumWF(float windowlow, float windowhi, float start, float end, int eventnr, float ymin, float ymax) {
	// plot waveforms of all channels for a given event number eventnr and add the determined integration windwos to the plot

	TString name(Form("waveforms_event__%04d", eventnr));
	TCanvas* intwinc = new TCanvas(name.Data(), name.Data(), 1600, 1000);
	intwinc->Divide(ceil(nchannels / 3), 3, 0, 0);

	if (eventnr > 0) eventnr -= 1; //rundata->At() counter starts at 0 which contains event 1 (wavecatcher event numbering starts with 1)

	for (int i = 0; i < nchannels; i++) {
		TH1F* his;
		his = ((TH1F*)rundata->At(eventnr * nchannels + i));
		int* windowind = GetIntWindow(his, windowlow, windowhi, start, end, i);
		// create lines to indicate the integration window
		TLine* low = new TLine(his->GetXaxis()->GetBinCenter(windowind[1]), -5, his->GetXaxis()->GetBinCenter(windowind[1]), 10);
		low->SetLineColor(2);
		TLine* hi = new TLine(his->GetXaxis()->GetBinCenter(windowind[2]), -2, his->GetXaxis()->GetBinCenter(windowind[2]), 3);
		hi->SetLineColor(2);
		TLine* zero = new TLine(0, 0, 320, 0); // draw line at x=0 to check if baseline correction worked
		zero->SetLineColor(1);
		delete[] windowind;

		TLine* baselinel = new TLine(baseline_correction_result[eventnr * nchannels + i][2], -1, baseline_correction_result[eventnr * nchannels + i][2], 1);
		baselinel->SetLineColor(6);
		baselinel->SetLineWidth(2);
		TLine* baselineh = new TLine(baseline_correction_result[eventnr * nchannels + i][3], -1, baseline_correction_result[eventnr * nchannels + i][3], 1);
		baselineh->SetLineColor(6);
		baselineh->SetLineWidth(2);
		TLine* baseline = new TLine(baseline_correction_result[eventnr * nchannels + i][2], 0, baseline_correction_result[eventnr * nchannels + i][3], 0);
		baseline->SetLineColor(6);


		// draw to canvas
		intwinc->cd(i + 1);
		his->Draw();
		if (ymin != 0. && ymax != 0.) his->GetYaxis()->SetRangeUser(ymin, ymax); //for better comparison fix y range
		low->Draw("same");
		hi->Draw("same");
		zero->Draw("same");
		baselinel->Draw("same"); 
		baselineh->Draw("same");
		baseline->Draw("same");
	}
	intwinc->Update();

	stringstream namess;
	namess << name.Data() << ".pdf";
	intwinc->SaveAs(namess.str().c_str());
}


TH1F* ReadRun::ChargeSpectrum(int channel, float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins) {
	// integrate all pulses in range (start, end) from t_max - windowlow to t_max + windowhi for a given channel and return the charge histogram with x range (rangestart, rangeend) and the number of bins nbins 

	TString name(Form("channel__%02d", channel));
	TH1F* h1 = new TH1F(name.Data(), name.Data(), nbins, rangestart, rangeend);

	for (int j = 0; j < nevents; j++) {
		TH1F* his = ((TH1F*)rundata->At(j * nchannels + channel));
		int* windowind = GetIntWindow(his, windowlow, windowhi, start, end, channel);	// find integration window
		h1->Fill(his->Integral(windowind[1], windowind[2]));					// fill charge spectrum
		delete[] windowind;
	}

	//h1->Rebin();
	return h1;
}


void ReadRun::PrintChargeSpectrum(float windowlow, float windowhi, float start, float end, float rangestart, float rangeend, int nbins, float fitrangestart, float fitrangeend) {
	// print ReadRun::ChargeSpectrum for all channels 

	if (fitrangestart == 0.) fitrangestart = rangestart;
	if (fitrangeend == 0.) fitrangeend = rangeend;

	TCanvas* chargec = new TCanvas("charge spectra", "charge spectra", 1600, 1000);
	chargec->Divide(ceil(nchannels / 3), 3, 0, 0);
	for (int i = 0; i < nchannels; i++) {
		TH1F* his;
		his = ChargeSpectrum(i, windowlow, windowhi, start, end, rangestart, rangeend, nbins);
		chargec->cd(i + 1);

		//Fitf fitf;
		//TF1* f = new TF1("fitf", fitf, fitrangestart, fitrangeend, 7);
		//f->SetLineColor(3);
		//f->SetParName(0, "N0");					f->SetParameter(0, his->Integral()/3.);
		//f->SetParName(1, "#mu");				f->SetParameter(1, 2.);
		//f->SetParName(2, "#lambda");			f->SetParameter(2, .15); //0.2 or 3
		//f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 3.2);//3.6
		//f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, .12);		f->SetParLimits(4, 1.e-9, 1.e3);	//f->FixParameter(4, 0.1);
		//f->SetParName(5, "Gain");				f->SetParameter(5, 10.);	//f->FixParameter(5, 10.);
		//f->SetParName(6, "Pedestal");			f->SetParameter(6, 2.);//7.);									//f->FixParameter(6, 0.);
		

		Fitf_biased fitf_biased;
		TF1* f = new TF1("fitf", fitf_biased, fitrangestart, fitrangeend, 9);
		f->SetLineColor(3);

		//+-1ns 44V
		f->SetParName(0, "N0");					f->SetParameter(0, his->Integral());
		f->SetParName(1, "#mu");				f->SetParameter(1, 1.6);// 1.6);
		f->SetParName(2, "#lambda");			f->SetParameter(2, .1); //0.2 or 3
		f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 3.);//3.6
		f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, 1.);		f->SetParLimits(4, 1.e-9, 1.e3);	//f->FixParameter(4, 0.1);
		f->SetParName(5, "Gain");				f->SetParameter(5, 10);	//f->FixParameter(5, 10.);
		f->SetParName(6, "Pedestal");			f->SetParameter(6, 2.);
		f->SetParName(7, "norm_0");				f->SetParameter(7, 0.8); //f->FixParameter(7, 1.);
		f->SetParName(8, "x_0");				f->SetParameter(8, 5.);
		
		//+-2ns
		//f->SetParName(0, "N0");					f->SetParameter(0, his->Integral());
		//f->SetParName(1, "#mu");				f->SetParameter(1, 2.);
		//f->SetParName(2, "#lambda");			f->SetParameter(2, .1); //
		//f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 7.5);//7.5
		//f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, 2.);		f->SetParLimits(4, 1.e-9, 1.e3);	//f->FixParameter(4, 0.1);
		//f->SetParName(5, "Gain");				f->SetParameter(5, 18.);//25	//f->FixParameter(5, 10.);
		//f->SetParName(6, "Pedestal");			f->SetParameter(6, 4.);//2
		//f->SetParName(7, "norm_0");				f->SetParameter(7, 0.8); //f->FixParameter(7, 1.);
		//f->SetParName(8, "x_0");				f->SetParameter(8, 8.5);

		//+-3ns
		//f->SetParName(0, "N0");					f->SetParameter(0, his->Integral());
		//f->SetParName(1, "#mu");				f->SetParameter(1, 2.2);
		//f->SetParName(2, "#lambda");			f->SetParameter(2, .1); //
		//f->SetParName(3, "#sigma_{0}");			f->SetParameter(3, 9.5);//7.5
		//f->SetParName(4, "#sigma_{1}");			f->SetParameter(4, 1.);		f->SetParLimits(4, 1.e-9, 1.e3);	//f->FixParameter(4, 0.1);
		//f->SetParName(5, "Gain");				f->SetParameter(5, 26.);//25	//f->FixParameter(5, 10.);
		//f->SetParName(6, "Pedestal");			f->SetParameter(6, 6.);//2
		//f->SetParName(7, "norm_0");				f->SetParameter(7, 0.9); //f->FixParameter(7, 1.);
		//f->SetParName(8, "x_0");				f->SetParameter(8, 11.);

		if (i != 8) {
			cout << "\n\n---------------------- Fit for channel " << i << " ----------------------\n";
			his->Fit(f, "RM");
		}
		his->GetYaxis()->SetTitle("#Entries");
		his->GetXaxis()->SetTitle("integral in mV#timesns");
		

		his->Draw();
		
		//f->Draw("same");
	}
	chargec->Update();
	chargec->SaveAs("ChargeSpectra.pdf");
}






// helper functions



stringstream ReadRun::list_files(const char* dirname, const char* ext) {
	// helper creating list of all bin files in directory
	stringstream ss;
	TSystemDirectory dir(dirname, dirname);
	TList* files = dir.GetListOfFiles();
	if (files) {
		TSystemFile* file;
		TString fname;
		TIter next(files);
		while ((file = (TSystemFile*)next())) {
			fname = file->GetName();
			if (!file->IsDirectory() && fname.EndsWith(ext)) {
				ss << fname.Data() << "\n";
				cout << fname.Data() << "\n";
			}
		}
		TIter next2(files);
		while ((file = (TSystemFile*)next2())) {
			fname = file->GetName();
			if (!file->IsDirectory() && !fname.EndsWith(ext) && fname.Contains(ext)) {
				ss << fname.Data() << "\n";
				cout << fname.Data() << "\n";
			}
		}
	}
	return ss;
}



TH1F* ReadRun::Getwf(int channelnr, int eventnr, int color) {
	TH1F* his;
	if (eventnr > 0) eventnr -= 1; //rundata->At() counter starts at 0 which contains event 1 (wavecatcher event numbering starts with 1)
	his = (TH1F*)rundata->At(eventnr * nchannels + channelnr);
	his->SetLineColor(color);
	his->SetMarkerColor(color);
	return his;
}


double* ReadRun::getx() {
	double* xvals = new double[1024];
	for (int i = 0; i < 1024; i++) {
		xvals[i] = static_cast<double>(SP) * static_cast<double>(i);
	}
	return xvals;
}


double* ReadRun::gety(int channelnr, int eventnr) {
	TH1F* his = Getwf(channelnr, eventnr);
	double* yvals = new double[1024];
	for (int i = 0; i < his->GetNbinsX(); i++) {
		yvals[i] = his->GetBinContent(i);
	}
	return yvals;
}


double* ReadRun::gety(TH1F* his) {
	double* yvals = new double[1024];
	for (int i = 0; i < his->GetNbinsX(); i++) {
		yvals[i] = his->GetBinContent(i);
	}
	return yvals;
}


void ReadRun::Convolute(double* &result, double* first, double* second, int size1, int size2) {
	// faster if size1<size2
	//   double *result = new double[size2+1];
	for (int i = 0; i < size2; i++) {
		result[i] = 0.;
		for (int j = 0; j < TMath::Min(size1, i); j++) {
			result[i] += first[j] * second[i - j];
		}
	}
}



void ReadRun::SmoothArray(double* &ar, int nbins, double sigma, bool doconv) {
	//apply smoothing array of double with length nbins

	double* artmp = new double[nbins];
	for (int i = 0; i < nbins; i++) artmp[i] = ar[i];

	if (doconv) {
		// convolution with gauss with sigma (experimental, not yet tested)
		double* gauss = new double[nbins];
		
		double sum = 0.;

		for (int i = 0; i < nbins; i++) {
			gauss[i] = TMath::Exp(-1. * TMath::Power((static_cast<double>(i) * SP -5*sigma/*-> offset?*/), 2.) / (2. * sigma * sigma)) / (sigma * 2.506628);
			sum += gauss[i];
		}

		for (int i = 0; i < nbins; i++) {
			gauss[i] /= sum;
		}

		Convolute(ar, artmp, gauss, nbins, nbins);
		delete[] gauss;
	}
	else {
		// calculate running average from -sigma until +sigma
		for (int i = 0; i < nbins; i++) {
			double mean1 = 0.; 
			int nmn = 0;
			for (int k = -1 * static_cast<int>(floor(sigma)); k <= static_cast<int>(ceil(sigma)); k++) {
				if (i + k >= 0 && i + k < nbins) {
					mean1 += artmp[i + k];
					nmn++;
				}
			}
			if (nmn != 0.) {
				ar[i] = mean1 / static_cast<double>(nmn);
			}
		}
	}
	delete[] artmp;
}



//TODO: 
// 1: new function plotting all wfs with integration windows for event x				<- done
// 2: new function which takes th1f , finds peaks and returns integration windows		<- done
// 3: new function to analyze sum spectra and plot and use the results for integrals	<-
// 4: baseline fit with chi2 cut														<- sort of done min(RMS(slope{-range,range}))
// 5: Fit charge spectrum																<- 
// 6: DC probability																	<- 
// 7: implement method to discard individual waveforms (not sure where this is needed)	<-
// 8: compile as library																<-
// 9: implement method to determine number of "events" over threshold					<-

