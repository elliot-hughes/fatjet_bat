// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__FATJETMODEL__H
#define __BAT__FATJETMODEL__H

//#include <BAT/BCModel.h>
#include <BAT/BCMTF.h>

#include <string>
#include <vector>
#include <map>
#include <TH1D.h>
#include <TString.h>

// This is a fatjetModel header file.
// Model source code is located in file fatjetModel/fatjetModel.cxx

// ---------------------------------------------------------
class fatjetModel : public BCMTF
{

public:

    // Constructor
    fatjetModel(const std::string& name);

    // Destructor
    ~fatjetModel();

	// Custom stuff:
	/// CDF things:
	TH1D* make_cdf(const TH1D* h, TString name="");
	double median_from_cdf(const TH1D* cdf);
	double evaluate_cdf(const TH1D* cdf, const double x);
	void differentiate_cdf(const TH1D* cdf, TH1D* out, double amplitude, double shift, double stretch);
	
	/// Containers:
	void fill_hists();
	TH1D* get_histogram(std::string channel_name, std::string process_name);
	void fill_cdfs();
	TH1D* get_cdf(std::string channel_name, std::string process_name);
	
	
	/// Morphing:
	TH1D* morph_template(std::string channel_name, std::string process_name, double amplitude, double shift, double stretch);
	
	/// Other:
	double get_parameter_value(const std::vector<double>& parameters, std::string name);
	
	// Overload LogLikelihood to implement model
	double LogLikelihood(const std::vector<double>& parameters);
	
	// Limit setting:
	TH1D* make_toy(std::string channel);
	std::vector<double> get_limit(bool observed, std::string signal_process, int n);
	std::vector<double> get_observed_limit(std::string signal_process, int n);
	std::vector<double> get_expected_limit(std::string signal_process, int n);
	
	void PrintStack(int channelindex, const std::vector<double>& parameters, const std::string& filename = "stack.pdf", const std::string& options = "e1b0stack");


private:
	std::map<std::string, TH1D*> hists;
	std::map<std::string, TH1D*> cdfs;
	std::string input_file = "theta_plots_sim_sig_sb.root";
	std::vector<std::string> processes = {"DATA", "QCD", "TTbar", "Ms100", "Ms150", "Ms200", "Ms250", "Ms300", "Ms400", "Ms500"};
	std::vector<std::string> mss = {"100", "150", "200", "250", "300", "400", "500"};
	std::vector<std::string> processes_sig = {"Ms100", "Ms150", "Ms200", "Ms250", "Ms300", "Ms400", "Ms500"};
//	std::vector<std::string> channels = {"SR", "SB"};
	std::vector<std::string> channels = {"SR"};
};
// ---------------------------------------------------------

#endif
