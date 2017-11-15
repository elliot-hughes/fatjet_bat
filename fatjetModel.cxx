// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "fatjetModel.h"

#include <BAT/BCMath.h>
#include <BAT/BCMTFProcess.h>
#include <BAT/BCMTFChannel.h>
#include <BAT/BCMTFTemplate.h>

#include <string>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <THStack.h>
#include <TString.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <iostream>
#include <sstream>

// ---------------------------------------------------------
fatjetModel::fatjetModel(const std::string& name) :
	BCMTF(name)
{
    // Define parameters here in the constructor. For example:
    // AddParameter("mu",-2,1,"#mu");
    // And set priors, if using built-in priors. For example:
    // GetParamater("mu").SetPrior(new BCPriorGaus(-1, 0.25));
    this->fill_hists();
    this->fill_cdfs();
    
}

// ---------------------------------------------------------
fatjetModel::~fatjetModel()
{
    // destructor
    for (auto const& i : hists) delete i.second;
    for (auto const& i : cdfs) delete i.second;
}

// ---------------------------------------------------------
TH1D* fatjetModel::make_cdf(const TH1D* h, TString name) {
	if (!name[0]) name = h->GetName() + TString("_cdf");
	TH1D* cdf = (TH1D*) h->Clone(name);
	
	for (int ibin = 0; ibin <= h->GetNbinsX() + 1; ibin++) {
		double integral = 0;
		for (int jbin = 0; jbin < ibin; jbin ++) integral += h->GetBinContent(jbin);
		cdf->SetBinContent(ibin, integral);
		cdf->SetBinError(ibin, 0);
	}
//	std::cout << name << std::endl;
	return cdf;
}

// ---------------------------------------------------------
double fatjetModel::median_from_cdf(const TH1D* cdf) {
	double max = cdf->GetMaximum();
	double curval, preval = cdf->GetBinContent(0);
	int nbins = cdf->GetNbinsX() + 1;
	for(int i = 0; i <= nbins+1; ++i) {
		curval = cdf->GetBinContent(i);
		if(curval >= max/2.0 && preval < max/2.0) {
			return cdf->GetXaxis()->GetBinLowEdge(i);
		}
	}
//	std::cerr << "Could not find the median!!!" << std::endl;
	return cdf->GetXaxis()->GetBinUpEdge(nbins);
}

// ---------------------------------------------------------
double fatjetModel::evaluate_cdf(const TH1D* cdf, const double x) {
	int bin = cdf->GetXaxis()->FindBin(x);
	double edge_low = cdf->GetXaxis()->GetBinLowEdge(bin);
	double edge_high = cdf->GetXaxis()->GetBinUpEdge(bin);
	double bin_value = cdf->GetBinContent(bin);
	double nextbin_value = cdf->GetBinContent(bin + 1);
	
	double evaluation = (nextbin_value - bin_value)/(edge_high - edge_low)*(x - edge_low) + bin_value;
//	if(DEBUG) std::cout << "x=" << x << "; bin=" << bin << "; edge_low=" << edge_low
//		<< "; edge_high=" << edge_high << "; bin_value=" << bin_value
//		<< "; nextbin_value=" << nextbin_value << "; evaluation=" << evaluation << std::endl;
	
	return evaluation;
}

//string fatjetModel::get_process_name_from_template(BCMTFTemplate*) {
//	
//}

// ---------------------------------------------------------
void fatjetModel::fill_hists() {
	// Read file:
	TFile* tfile = TFile::Open(TString(input_file), "READ");
	if (!tfile || !tfile->IsOpen()) {
		BCLog::OutError(Form("Could not open file %s.", input_file.c_str()));
		return;
	}
	// Fetch histograms:
	for (unsigned i = 0; i < channels.size(); ++i) {
		std::string channel = channels[i];
		for (unsigned j = 0; j < processes.size(); ++j) {
			std::string process = processes[j];
			TString name = "mavg" + channel + "__" + process;
			std::string key = channel + "_" + process;
			TH1D* h = (TH1D*) tfile->Get(name);
			if (!h) {
				BCLog::OutError(Form("Could not find histogram named %s.", name.Data()));
			}
			else {
				hists[key] = h;
			}
		}
	}
//	return;
//	// Loop over channels: ("SR", "SB")
//	for (int ichannel = 0; ichannel < BCMTF::GetNChannels(); ++ichannel) {
//		BCMTFChannel* channel = BCMTF::GetChannel(ichannel);
//		std::string channel_name = channel->GetName();		// "SR" or "SB"
//		
//		// Assign DATA histogram:
//		BCMTFTemplate* temp_data = channel->GetData();
//		hists[channel_name + "_DATA"] = temp_data->GetHistogram();
//		
//		// Assign other process histograms:
//		for (int itemp = 0; itemp < BCMTF::GetNProcesses(); ++itemp) {
//			BCMTFTemplate* temp = channel->GetTemplate(itemp);
//			std::string process_name = "None";
//			// TTbar:
//			if (temp->GetProcessName() == "normTTbar") process_name = "TTbar";
//			// QCD:
//			else if (temp->GetProcessName() == "normQCDSR" || temp->GetProcessName() == "normQCDSB") process_name = "QCD";
//			// Signal:
//			else {
//				for (unsigned isig = 0; isig < processes_sig.size(); ++isig) {
//					if (temp->GetProcessName() == "norm" + processes_sig[isig]) {
//						process_name = processes_sig[isig];
//						break;
//					}
//				}
//			}
//			TH1D* h = temp->GetHistogram();
//			if (h) hists[channel_name + "_" + process_name] = h;
//			else std::cout << "This histogram is empty: " << channel_name << "  " << temp->GetProcessName() << std::endl;
//		}
//	}
	
	// Print stuff:
	unsigned n = hists.size();
	if (n > 0) {
		for (auto const& i : hists) {
//			BCLog::OutSummary(Form("  * %s", i.first.Data()));
			std::cout << i.first << std::endl;
		}
//		BCLog::OutSummary("");
		std::cout << std::endl;
	}
	else {
//		BCLog::OutError("hists is empty!");
		std::cout << "hists is empty!" << std::endl;
	}
	return;
}

// ---------------------------------------------------------
TH1D* fatjetModel::get_histogram(std::string channel_name, std::string process_name) {
	return hists[channel_name + "_" + process_name];
}

// ---------------------------------------------------------
void fatjetModel::fill_cdfs() {
	std::cout << "Filling cdfs" << std::endl;
	TFile* tfile_sb = TFile::Open("analysis_plots_sb_f1.root");
//	TFile* tfile_sb = TFile::Open("analysis_plots_sbb_f1.root");
	TFile* tfile_sr = TFile::Open("analysis_plots_sig_f1.root");
	cdfs["SR_QCD"] = (TH1D*) tfile_sr->Get("cdf_jetht");
	cdfs["SB_QCD"] = (TH1D*) tfile_sb->Get("cdf_jetht");
	cdfs["SR_TTbar"] = (TH1D*) tfile_sr->Get("cdf_ttbar");
	cdfs["SB_TTbar"] = (TH1D*) tfile_sb->Get("cdf_ttbar");
	for (unsigned i = 0; i < mss.size(); ++i) {
		cdfs["SB_Ms" + mss[i]] = (TH1D*) tfile_sb->Get(TString("cdf_sq") + mss[i] + "to4j");
		cdfs["SR_Ms" + mss[i]] = (TH1D*) tfile_sr->Get(TString("cdf_sq") + mss[i] + "to4j");
	}
	
//	std::cout << "Filling cdfs" << std::endl;
//	for (auto const& i : hists) {
//		std::string name = i.first;
//		if (name == "SR_DATA" || name == "SB_DATA") continue;
////		std::cout << name << std::endl;
//		cdfs[name] = make_cdf(i.second);
//	}
////	std::cout << "done" << std::endl;
	return;
}

// ---------------------------------------------------------
TH1D* fatjetModel::get_cdf(std::string channel_name, std::string process_name) {
	return cdfs[channel_name + "_" + process_name];
}

// ---------------------------------------------------------
void fatjetModel::differentiate_cdf(const TH1D* cdf, TH1D* out, double amplitude, double shift, double stretch) {
//	assert(stretch > 0.0);
	
	double peak = median_from_cdf(cdf);
//	std::cout << peak << "  " << amplitude << "  " << shift << "  " << stretch << std::endl;		// DEBUG
	
	for(int i = 0; i <= out->GetNbinsX() + 1; i++) {
		double lo = out->GetXaxis()->GetBinLowEdge(i);
		double hi = out->GetXaxis()->GetBinUpEdge(i);
		
		double loprime = 1./stretch*(lo - peak - shift) + peak;
		double hiprime = 1./stretch*(hi - peak - shift) + peak;
		
		double y1 = evaluate_cdf(cdf, loprime);
		double y2 = evaluate_cdf(cdf, hiprime);
		
		double val = amplitude*(y2 - y1);

//		if(DEBUG) cout << "lo=" << lo << "; hi=" << hi << "; loprime=" << loprime << "; hiprime=" << hiprime
//			<< "; y1=" << y1 << "; y2=" << y2 << "; val=" << val << endl;
		
		out->SetBinContent(i, val);
		out->SetBinError(i, 0.);
	}
	
	return;
}

// ---------------------------------------------------------
TH1D* fatjetModel::morph_template(std::string channel_name, std::string process_name, double amplitude, double shift, double stretch) {
	TH1D* cdf = get_cdf(channel_name, process_name);
	TH1D* h = get_histogram(channel_name, process_name);
	TH1D* h_morphed = (TH1D*) h->Clone(h->GetName() + TString("_morphed"));
	differentiate_cdf(cdf, h_morphed, amplitude, shift, stretch);
//	delete cdf;
	return h_morphed;
}

double fatjetModel::get_parameter_value(const std::vector<double>& parameters, std::string name) {
	for (unsigned i = 0; i < parameters.size(); ++i) {
//		double par_value = parameters[i];
//		std::string par_name = BCMTF::GetParameter(i).GetName();
//		if (par_name == name) return par_value;
		if (BCMTF::GetParameter(i).GetName() == name) return parameters[i];
	}
	std::cout << "ERROR getting parameter value for " << name << std::endl;
	return 0;
}

// ---------------------------------------------------------
double fatjetModel::LogLikelihood(const std::vector<double>& parameters)
{
	// This returns the log of the conditional probability p(data|pars)
	// This is where you define your model.
	// BCMath contains many functions you will find helpful
	
//	return BCMTF::LogLikelihood(parameters);
	
//	std::cout << "lik start" << std::endl;
	double result = 0;

	// Loop over channels:
	for (int ichannel = 0; ichannel < BCMTF::GetNChannels(); ++ichannel) {
		BCMTFChannel* channel = BCMTF::GetChannel(ichannel);
		std::string channel_name = channel->GetName();		// "SR" or "SB"
//		std::cout << channel_name << std::endl;
		
		// Get data template and histogram:
		BCMTFTemplate* temp_data = channel->GetData();
		TH1D* h_data = temp_data->GetHistogram();
//		
//		// Get ttbar and QCD templates:
//		BCMTFTemplate* temp_ttbar;
//		BCMTFTemplate* temp_qcd;
//		for (int itemp = 0; itemp < 3; ++itemp) {		// I don't know how to get ntemplates ...
//			BCMTFTemplate* temp = channel->GetTemplate(itemp);
////			std::cout << temp->GetProcessName() << std::endl;
//			if (temp->GetProcessName() == "normTTbar") temp_ttbar = temp;
//			else if (temp->GetProcessName() == "normQCD" + channel_name) temp_qcd = temp;
//		}
//		/// TODO: I really need to add a check here.
		
		// Get parameters:
		double normQCD = get_parameter_value(parameters, "normQCD" + channel_name);
		double normTTbar = get_parameter_value(parameters, "normTTbar");
		double shiftQCD = get_parameter_value(parameters, "shiftQCD");
		double shiftTTbar = get_parameter_value(parameters, "shiftTTbar");
		double stretchQCD = get_parameter_value(parameters, "stretchQCD");
		double stretchTTbar = get_parameter_value(parameters, "stretchTTbar");
		
//		std::cout << "here 268" << std::endl;
		
		// Shift and stretch the templates:
//		TH1D* h_ttbar = get_histogram(channel_name, "TTbar");
//		TH1D* h_qcd = get_histogram(channel_name, "QCD");
		TH1D* h_ttbar_morphed = morph_template(channel_name, "TTbar", normTTbar, shiftTTbar, stretchTTbar);
		TH1D* h_qcd_morphed = morph_template(channel_name, "QCD", normQCD, shiftQCD, stretchQCD);
		
		TH1D* h_sig_morphed;
		int nsig = 0;
		for (unsigned isig = 0; isig < processes_sig.size(); ++isig) {
			double normSig = get_parameter_value(parameters, "norm" + processes_sig[isig]);
			if (normSig > 0) {
				if (nsig == 0) h_sig_morphed = morph_template(channel_name, processes_sig[isig], normSig, shiftTTbar, stretchTTbar);
				else h_sig_morphed->Add(morph_template(channel_name, processes_sig[isig], normSig, shiftTTbar, stretchTTbar));
				nsig += 1;
			}
		}
		
//		std::cout << "here 277" << std::endl;
		
//		TH1D* cdf_ttbar = make_cdf(h_ttbar, "");
//		TH1D* cdf_qcd = make_cdf(h_qcd, "");
//		TH1D* h_ttbar_morphed = (TH1D*) h_ttbar->Clone(h_ttbar->GetName() + TString("_morphed"));
//		differentiate_cdf(cdf_ttbar, h_ttbar_morphed, normTTbar, shiftTTbar, stretchTTbar);
////		std::cout << h_ttbar->Integral() << "   " << h_ttbar_morphed->Integral() << std::endl;
//		TH1D* h_qcd_morphed = (TH1D*) h_qcd->Clone(h_qcd->GetName() + TString("_morphed"));
//		differentiate_cdf(cdf_qcd, h_qcd_morphed, normQCD, shiftQCD, stretchQCD);
		
		// Calculate likelihood component:
		int nbins = h_data->GetNbinsX();
		int nbins_flagged = 0;
		for (int ibin = 1; ibin <= nbins; ++ibin) {
			double m = h_data->GetXaxis()->GetBinUpEdge(ibin);
			if (m > 900) break;
			double obs = h_data->GetBinContent(ibin);
			double exp = h_qcd_morphed->GetBinContent(ibin) + h_ttbar_morphed->GetBinContent(ibin);
			if (nsig > 0) exp += h_sig_morphed->GetBinContent(ibin);
			if (exp <= 0) {
				nbins_flagged += 1;
				exp = 0.000001;
			}
			double result_component = obs*log(exp) - exp - TMath::LnGamma(obs + 1);
//			std::cout << "obs = " << obs << ", exp = " << exp << ", L = " << result_component << std::endl;		// DEBUG
			result += result_component;
//			std::cout << "obs = " << obs << ", exp = " << exp << ", L = " << BCMath::LogPoisson(obs, exp) << std::endl;		// THIS FUNCTION GIVES CRAZY RESULTS: compare BCMath::LogPoisson(228,899) = -361.832, BCMath::LogPoisson(228,899.001) = -254.732
//			result += BCMath::LogPoisson(obs, exp);
		}
		
//		if (bins_used > 0) result_component *= bin_frac;
//		else result_component = 0;
//		result += result_component;
		
//		std::cout << "here 311" << std::endl;
		
//		delete cdf_ttbar;
//		delete cdf_qcd;
		delete h_ttbar_morphed;
		delete h_sig_morphed;
		delete h_qcd_morphed;
		
////		std::cout << "normQCD = " << normQCD << std::endl;
////		std::cout << "normTTbar = " << normTTbar << std::endl;
////		std::cout << "shiftQCD = " << shiftQCD << std::endl;
////		std::cout << "shiftTTbar = " << shiftTTbar << std::endl;
////		std::cout << "stretchQCD = " << stretchQCD << std::endl;
////		std::cout << "stretchTTbar = " << stretchTTbar << std::endl;
////		std::cout << "nbins_flagged = " << nbins_flagged << std::endl;
//		std::cout << "L = " << result << std::endl;		// DEBUG
	}
//	std::cout << parameters.size() << std::endl;
//	std::cout << hists.size() << std::endl;
//	std::cout << cdfs.size() << std::endl;
	
//	std::cout << result << std::endl;
//	std::cout << "lik done" << std::endl;

	
	return result;
}

// ---------------------------------------------------------
TH1D* fatjetModel::make_toy(std::string channel) {
	TH1D* h_bkg = (TH1D*) get_histogram(channel, "QCD")->Clone(TString(channel) + "_bkg");
	h_bkg->Add(get_histogram(channel, "TTbar"));
	TH1D* h_toy = (TH1D*) h_bkg->Clone(TString(channel) + "_toy");
	h_toy->Reset();
	gRandom->SetSeed(0);
	h_toy->FillRandom(h_bkg, h_bkg->Integral());
	TCanvas* tc = new TCanvas();
	h_toy->Draw("e0");
	tc->SaveAs("toy.pdf");
	return h_toy;
}

// ---------------------------------------------------------
std::vector<double> fatjetModel::get_limit(bool observed, std::string signal_process, int n) {
	std::vector<double> results;
	std::vector<double> output;
	
	// Set signal strength of particular signal to 1:
	for (unsigned isig = 0; isig < processes_sig.size(); ++isig) {
		if (signal_process != processes_sig[isig]) this->GetParameter("norm" + processes_sig[isig]).Fix(0.0);
		else this->GetParameter("norm" + processes_sig[isig]).Unfix();
	}
	
	// Do n iterations:
	for (int i = 0; i < n; ++i) {
		// Set data or toy:
		if (observed) {
			for (unsigned i = 0; i < channels.size(); ++i) {
				this->SetData(channels[i], *(this->get_histogram(channels[i], "DATA")));
			}
		}
		else {
			for (unsigned i = 0; i < channels.size(); ++i) {
				this->SetData(channels[i], *(make_toy(channels[i])));
			}
		}
	
		// Get limit:
		this->SetRandomSeed(0);		// This changes the random seed, otherwise every marginalization would give the same result.
		this->MarginalizeAll();
		results.push_back(this->GetMarginalized("norm" + signal_process).GetQuantile(0.95));
	}
	
	// Calculate the output:
	/// Sort the results:
	std::sort(results.begin(), results.end());
	/// Calculate the central values and bands:
	//// This procedure is intended to replicate Theta (which replicates combined?).
	double median = results[n/2];
	if (observed) {
		double width = results[int(0.75*n)] - results[int(0.25*n)];
		std::vector<double> results_trimmed;
		for (unsigned i = 0; i < results.size(); ++i){
			if (results[i] >= median - width && results[i] <= median + width) results_trimmed.push_back(results[i]);
		}
		double mean = 0;
		for (unsigned i = 0; i < results_trimmed.size(); ++i) mean += results_trimmed[i];
		mean /= results_trimmed.size();
		double sigma = 0;
		if (results_trimmed.size() > 1) {
			for (unsigned i = 0; i < results_trimmed.size(); ++i) sigma += pow(results_trimmed[i] - mean, 2);
			sigma /= results_trimmed.size() - 1;
			sigma = sqrt(sigma);
		}
		output.push_back(mean);
		output.push_back(mean - sigma);
		output.push_back(mean + sigma);
		output.push_back(mean - 2*sigma);
		output.push_back(mean + 2*sigma);
	}
	else {
		output.push_back(median);
		output.push_back(results[int(0.16*n)]);		// -1 sigma
		output.push_back(results[int(0.84*n)]);		// +1 sigma
		output.push_back(results[int(0.025*n)]);		// -2 sigma
		output.push_back(results[int(0.975*n)]);		// +2 sigma
	}
	return output;
}

// ---------------------------------------------------------
std::vector<double> fatjetModel::get_observed_limit(std::string signal_process, int n) {
	return get_limit(true, signal_process, n);
}

// ---------------------------------------------------------
std::vector<double> fatjetModel::get_expected_limit(std::string signal_process, int n) {
	return get_limit(false, signal_process, n);
}



// ---------------------------------------------------------
void fatjetModel::PrintStack(int channelindex, const std::vector<double>& parameters, const std::string& filename, const std::string& options) {
//	return BCMTF::PrintStack(channelindex, parameters, filename, options);
	
	// This is a modified version of BCMTF::PrintStack.
	
    // check if parameters are filled
    if (parameters.empty())
        return;

    // check options
    bool flag_logx   = false; // plot x-axis in log-scale
    bool flag_logy   = false; // plot y-axis in log-scale
    bool flag_bw     = false; // plot in black and white

    bool flag_sum    = false; // plot sum of all templates
    bool flag_stack  = false; // plot stack of templates

    bool flag_e0     = false; // do not draw error bars on data
    bool flag_e1     = false; // draw sqrt(N) error bars on data

    bool flag_b0     = false; // draw an error band on the expectation
    bool flag_b1     = false; // draw an error band on the number of events

    if (options.find("logx") < options.size())
        flag_logx = true;

    if (options.find("logy") < options.size())
        flag_logy = true;

    if (options.find("bw") < options.size())
        flag_bw = true;

    if (options.find("sum") < options.size())
        flag_sum = true;

    if (options.find("stack") < options.size())
        flag_stack = true;

    if (options.find("e0") < options.size())
        flag_e0 = true;

    if (options.find("e1") < options.size())
        flag_e1 = true;

    if (options.find("b0") < options.size())
        flag_b0 = true;

    if (options.find("b1") < options.size())
        flag_b1 = true;

    if (!flag_e0)
        flag_e1 = true;

    // check if MCMC ran
    if (!(GetMarginalizationMethod() == BCIntegrate::kMargMetropolis)) {
        flag_b0 = false;
        flag_b1 = false;
        BCLog::OutWarning("fatjetModel::PrintStack : Did not run MCMC. Error bands are not available.");
    }

    // get channel
    BCMTFChannel* channel = GetChannel(channelindex);
    std::string channel_name = channel->GetName();

    // create canvas
    TCanvas* c1 = new TCanvas();
    c1->cd();

    // set log or linear scale
    if (flag_logx) c1->SetLogx();
    if (flag_logy) c1->SetLogy();

    // get data histogram
    TH1D* hist_data = channel->GetData()->GetHistogram();

    // get number of bins
    int nbins = hist_data->GetNbinsX();

    // define sum of templates
    TH1D* hist_sum = new TH1D(*hist_data);
    hist_sum->SetLineColor(kBlack);
    for (int i = 1; i <= nbins; ++i)
        hist_sum->SetBinContent(i, 0);

//    // define error band
//    TH1D* hist_error_band = new TH1D(*hist_data);
//    hist_error_band->SetFillColor(kBlack);
//    hist_error_band->SetFillStyle(3005);
//    hist_error_band->SetLineWidth(1);
//    hist_error_band->SetStats(kFALSE);
//    hist_error_band->SetMarkerSize(0);

//    // get histogram for uncertainty band
//    TH2D* hist_uncbandexp = channel->GetHistUncertaintyBandExpectation();

//    // fill error band
//    if (flag_b0) {
//        for (int i = 1; i <= nbins; ++i) {
//            TH1D* proj = hist_uncbandexp->ProjectionY("_py", i, i);
//            if (proj->Integral() > 0)
//                proj->Scale(1.0 / proj->Integral());
//            double quantiles[3];
//            double sums[3] = {0.16, 0.5, 0.84};
//            proj->GetQuantiles(3, quantiles, sums);
//            hist_error_band->SetBinContent(i, 0.5 * (quantiles[2] + quantiles[0]));
//            hist_error_band->SetBinError(i, 0, 0.5 * (quantiles[2] - quantiles[0]));
//            delete proj;
//        }
//    }

    // create stack
    THStack* stack = new THStack("", "");

    // create a container of temporary histograms
    std::vector<TH1D*> histcontainer;

    // get number of templates
    unsigned int ntemplates = GetNProcesses();

    // loop over all templates
    for (unsigned int i = 0; i < ntemplates; ++i) {
		BCMTFTemplate* temp = channel->GetTemplate(i);
		if (!temp->GetHistogram()) continue;
		std::string temp_name = temp->GetProcessName();
		std::string process_name = "none";
		if (temp_name == "normTTbar") process_name = "TTbar";
		else if (temp_name == "normQCDSR" || temp_name == "normQCDSB") process_name = "QCD";
		
        // get histogram
//        TH1D* temphist = get_histogram(channel_name, process_name);
//        if (!temphist) continue;
        
        // Morph histogram:
        double norm = get_parameter_value(parameters, temp_name);
        double shift = get_parameter_value(parameters, "shift" + process_name);
        double stretch = get_parameter_value(parameters, "stretch" + process_name);
        TH1D* hist = morph_template(channel_name, process_name, norm, shift, stretch);
        std::cout << hist->GetName() << "  " << hist->Integral() << " +/- " << sqrt(hist->Integral()) << " events" << std::endl;

//        // create new histogram
//        TH1D* hist(0);
//        hist = new TH1D( *(temphist) );

        // set histogram style
		if (process_name == "TTbar") {
			hist->SetFillColor(kRed - 4);
			hist->SetFillStyle(3003);
			hist->SetMarkerSize(0);
			
//			h_fit_ttbar->SetLineWidth(2);
//			h_fjp_jetht->Rebin(nrebin);
//			h_fjp_jetht->SetLineWidth(2);
//			h_fjp_jetht->SetMarkerSize(1.2);
//			h_fjp_jetht->GetXaxis()->SetRangeUser(0, xmax);
//			h_fjp_jetht->GetXaxis()->SetNdivisions(406);
//			h_fjp_jetht->GetXaxis()->SetTitle(get_xtitle("mavg"));
		}
		else if (process_name == "QCD") {
			hist->SetFillColorAlpha(kBlue-10, 0.5);
			hist->SetFillStyle(1001);
			hist->SetMarkerSize(0);
			hist->SetLineStyle(2);
		}

        if (flag_bw) {
            hist->SetFillColor(0);
        }
		
		
        // add histogram to container (for memory management)
        histcontainer.push_back(hist);

        // add histogram to stack
        stack->Add(hist);
    }

    //draw data
	hist_data->GetXaxis()->SetRangeUser(0, 900);
	hist_data->GetXaxis()->SetNdivisions(405);
	hist_data->GetXaxis()->SetTitle("#bar{#it{m}}  [GeV]");
	hist_data->Draw("P0");

    // define variable for maximum in y-direction
    double ymax = 0;;

    if (flag_e1)
        ymax = hist_data->GetMaximum() + sqrt(hist_data->GetMaximum());
    else
        ymax = hist_data->GetMaximum();

    // set range user
    hist_data->GetYaxis()->SetRangeUser(channel->GetRangeYMin(), channel->GetRangeYMax());

    // draw stack
    if (flag_stack) {
        stack->Draw("SAMEHIST");
        if (stack->GetMaximum() > ymax)
            ymax = stack->GetMaximum();
    }

//    // draw error band on number of observed events
//    if (flag_b1) {
//        channel->CalculateHistUncertaintyBandPoisson();
//        TH1D* hist_temp = channel->CalculateUncertaintyBandPoisson(0.001, 0.999, kRed);
//        hist_temp->Draw("SAMEE2");
//        channel->CalculateUncertaintyBandPoisson(0.023, 0.977, kOrange)->Draw("SAMEE2");
//        channel->CalculateUncertaintyBandPoisson(0.159, 0.841, kGreen)->Draw("SAMEE2");

//        // get bin with maximum
//        int ymaxbin = hist_temp->GetMaximumBin();

//        if (hist_temp->GetBinContent(ymaxbin) + hist_temp->GetBinError(ymaxbin) > ymax)
//            ymax = hist_temp->GetBinContent(ymaxbin) + hist_temp->GetBinError(ymaxbin);
//    }

//    // draw error band on expectation
//    if (flag_b0) {
//        hist_error_band->Draw("SAMEE2");
//    }

    if (flag_sum)
        hist_sum->Draw("SAME");

    //draw data again
    if (flag_e0)
        hist_data->Draw("SAMEP0");

    if (flag_e1)
        hist_data->Draw("SAMEP0E");

    hist_data->GetYaxis()->SetRangeUser(0., 1.1 * ymax);

    // redraw the axes
    gPad->RedrawAxis();
    
    // Write parameters:
	std::ostringstream oss1;
	oss1 << "t#bar{t} norm = " << std::fixed << std::setprecision(2) << parameters[1];
	std::ostringstream oss2;
	oss2 << "t#bar{t} shift = " << std::fixed << std::setprecision(1) << parameters[3] << " GeV";
	std::ostringstream oss3;
	oss3 << "t#bar{t} stretch = " << std::fixed << std::setprecision(1);
	if (parameters[5] > 1) oss3 << "+";
	oss3 << 100*(parameters[5] - 1) << " \%";
	std::ostringstream oss4;
	oss4 << "QCD shift = " << std::fixed << std::setprecision(1) << parameters[2] << " GeV                                                    #bar{t}";
	std::ostringstream oss5;
	oss5 << "QCD stretch = " << std::fixed << std::setprecision(1);
	if (parameters[4] > 1) oss5 << "+";
	oss5 << 100*(parameters[4] - 1) << " \%                                                  #bar{t}";
	
	std::vector<TString> texts_par;
	texts_par.push_back("#bf{Fit parameters:}");
	texts_par.push_back(oss1.str());
	texts_par.push_back(oss2.str());
	texts_par.push_back(oss3.str());
	texts_par.push_back(oss4.str());
	texts_par.push_back(oss5.str());
	
	TString text = texts_par[texts_par.size() - 1];
	for (unsigned i = texts_par.size() - 1; i-- > 0;) {
		text = "#splitline{" + texts_par[i] + "}{" + text + "}";
	}
	TLatex* ttext = new TLatex(0, 0, text);
	ttext->SetNDC();		// Set text position to NDC coordinates.
	ttext->SetX(0.57);
	ttext->SetY(0.55);
	ttext->SetTextFont(42);
	ttext->SetTextSize(0.032);
	ttext->Draw();
	
    // print
    c1->Print(filename.data());

    // free memory
    for (unsigned int i = 0; i < histcontainer.size(); ++i) {
        TH1D* hist = histcontainer.at(i);
        delete hist;
    }
    delete stack;
    delete c1;
//    delete hist_error_band;
    delete hist_sum;
}

