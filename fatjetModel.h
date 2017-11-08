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
	
	void PrintStack(int channelindex, const std::vector<double>& parameters, const std::string& filename = "stack.pdf", const std::string& options = "e1b0stack");


private:
	std::map<std::string, TH1D*> hists;
	std::map<std::string, TH1D*> cdfs;
};
// ---------------------------------------------------------

#endif
