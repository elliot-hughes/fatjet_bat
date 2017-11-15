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
	
	
	///
	std::vector<std::string> get_process_names(std::string category="all");
	void let_signal_float(std::string signal_process="");
	
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
	void make_limit_plot(std::vector<std::vector<double>> expected_limits, std::vector<std::vector<double>> observed_limits, std::string prefix="");
	
	void PrintStack(int channelindex, const std::vector<double>& parameters, const std::string& filename = "stack.pdf", const std::string& options = "e1b0stack");


private:
	std::map<std::string, TH1D*> hists;
	std::map<std::string, TH1D*> cdfs;
	std::string input_file = "theta_plots_sim_sig_sb.root";
	std::vector<std::string> processes = {"DATA", "QCD", "TTbar", "Ms100", "Ms150", "Ms200", "Ms250", "Ms300", "Ms400", "Ms500"};
	std::vector<double> mss_double = {100, 150, 200, 250, 300, 400, 500};
	std::vector<std::string> mss = {"100", "150", "200", "250", "300", "400", "500"};
	std::vector<std::string> processes_sig = {"Ms100", "Ms150", "Ms200", "Ms250", "Ms300", "Ms400", "Ms500"};
//	std::vector<std::string> channels = {"SR", "SB"};
	std::vector<std::string> channels = {"SR"};
	std::vector<double> cross_sections = {1521.11, 574.981, 249.409, 121.416, 64.5085, 36.3818, 21.5949, 13.3231, 8.51615, 5.60471, 3.78661, 2.61162, 1.83537, 1.31169, 0.948333, 0.697075, 0.51848, 0.390303, 0.296128, 0.226118, 0.174599};
	std::vector<double> cross_section_errors = {0.154038, 0.149895, 0.147477, 0.146341, 0.144098, 0.142189, 0.140595, 0.142549, 0.139223, 0.138144, 0.136877, 0.138477, 0.136985, 0.135013, 0.134559, 0.133926, 0.133797, 0.133443, 0.132687, 0.132741, 0.132074};
};
// ---------------------------------------------------------

#endif
