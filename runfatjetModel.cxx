// ***************************************************************
// This file was created using the bat-project script
// for project fatjetModel.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCMTFChannel.h>
#include <BAT/BCEngineMCMC.h>
//#include <BAT/BCIntegrate.h>

#include "fatjetModel.h"

#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TMinuit.h>

using namespace std;

TString input_file = "theta_plots_sim_sig_sb.root";
//TString input_file = "theta_plots_sim_sig_sbb.root";
vector<TString> regions = {"SR", "SB"};
vector<TString> processes = {"DATA", "QCD", "TTbar", "Ms100", "Ms150", "Ms200", "Ms250", "Ms300", "Ms400", "Ms500"};
vector<TString> processes_bkg = {"QCD", "TTbar"};
vector<string> processes_sig = {"Ms100", "Ms150", "Ms200", "Ms250", "Ms300", "Ms400", "Ms500"};

map<TString,TH1D*> fetch_input_histograms(TString fname) {
	map<TString,TH1D*> histograms;
	
	/// Read file:
	TFile* tfile = TFile::Open(fname, "READ");
	if (!tfile || !tfile->IsOpen()) {
		BCLog::OutError(Form("Could not open file %s.", fname.Data()));
		return histograms;
	}
	/// Fetch histograms:
	for (unsigned i = 0; i < regions.size(); ++i) {
		TString region = regions[i];
		for (unsigned j = 0; j < processes.size(); ++j) {
			TString process = processes[j];
			TString name = "mavg" + region + "__" + process;
			TString key = process + "_" + region;
			TH1D* h = (TH1D*) tfile->Get(name);
			if (!h) {
				BCLog::OutError(Form("Could not find histogram named %s.", name.Data()));
			}
			else {
				histograms[name] = h;
			}
		}
	}
	return histograms;
}


int main()
{
	// set nicer style for drawing than the ROOT default
	BCAux::SetStyle();
	
	// open log file
//	BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);
	BCLog::OpenLog("log.txt", BCLog::detail, BCLog::summary);
	
	// create new fatjetModel object
	fatjetModel m("fatjetModel");
//	m.fill_hists();
//	m.fill_cdfs();
	
	// set precision
	m.SetPrecision(BCEngineMCMC::kMedium);
	
	// Fetch input histograms:
//	map<TString,TH1D*> input_histograms = fetch_input_histograms(input_file);
	/// Print out information about them:
//	unsigned ninput = input_histograms.size();
//	if (ninput > 0) {
//		BCLog::OutSummary(Form("Fetched %u input histograms:", ninput));
//		for (auto const& i : input_histograms) {
////			BCLog::OutSummary(Form("  * %s: %s", i.first.Data(), i.second->GetName()));
//			BCLog::OutSummary(Form("  * %s", i.first.Data()));
//		}
//		BCLog::OutSummary("");
//	}
//	else {
//		BCLog::OutError("No input histograms were found!");
//		return 0;
//	}
	
	// Define everything:
	/// Channels:
	vector<string> channels = {"SR", "SB"};
	channels = {"SR"};
//	channels = {"SB"};
	for (unsigned i = 0; i < channels.size(); ++i) m.AddChannel(channels[i]);
	
	/// Parameter configuration:
	//// Norms:
	double normQCD_min = 0.5;
	double normQCD_max = 2.0;
	double normTTbar_min = 0.5;
	double normTTbar_max = 1.5;
	double normSig_min = 0.0;
	double normSig_max = 1.5;
	//// Shifts:
	double shiftQCD_min = -40.0;
	double shiftQCD_max = 40.0;
	double shiftTTbar_min = -40.0;
	double shiftTTbar_max = 40.0;
	//// Stretches:
	double stretchQCD_min = 0.6;
	double stretchQCD_max = 1.5;
	double stretchTTbar_min = 0.6;
	double stretchTTbar_max = 1.5;
	
	
	/// Processes:
	//// Note: this order governs the stack order.
	//// Add QCD:
	for (unsigned i = 0; i < channels.size(); ++i) m.AddProcess("normQCD" + channels[i], normQCD_min, normQCD_max);
	//// Add TTbar:
	m.AddProcess("normTTbar", normTTbar_min, normTTbar_max);
	//// Add signals:
	for (unsigned i = 0; i < processes_sig.size(); ++i) m.AddProcess("norm" + processes_sig[i], normSig_min, normSig_max);
	
	/// Templates:
	for (unsigned i = 0; i < channels.size(); ++i) {
		string channel = channels[i];
		m.SetTemplate(channel, "normTTbar", *m.get_histogram(channel, "TTbar"), 1.0, m.get_histogram(channel, "TTbar")->Integral());
		m.SetTemplate(channel, "normQCD" + channel, *m.get_histogram(channel, "QCD"), 1.0, m.get_histogram(channel, "QCD")->Integral());
		for (unsigned j = 0; j < processes_sig.size(); ++j) {
			string sig_name = processes_sig[j];
			TH1D* sig_hist = m.get_histogram(channel, sig_name);
			m.SetTemplate(channel, "norm" + sig_name, *sig_hist, 1.0, sig_hist->Integral());
		}
	}
	
	
	m.AddParameter("shiftQCD", shiftQCD_min, shiftQCD_max);
	m.AddParameter("shiftTTbar", shiftTTbar_min, shiftTTbar_max);
	m.AddParameter("stretchQCD", stretchQCD_min, stretchQCD_max);
	m.AddParameter("stretchTTbar", stretchTTbar_min, stretchTTbar_max);
	
//	m.GetParameter("normQCDSR").SetPriorConstant();
	for (unsigned i = 0; i < channels.size(); ++i) m.GetParameter("normQCD" + channels[i]).SetPriorConstant();
	m.GetParameter("normTTbar").SetPriorConstant();
	for (unsigned i = 0; i < processes_sig.size(); ++i) m.GetParameter("norm" + processes_sig[i]).SetPriorConstant();
//	m.GetParameter("normQCDSR").SetPrior(new BCGaussianPrior(1.0, 100.0));
//	m.GetParameter("normQCDSB").SetPrior(new BCGaussianPrior(1.0, 100.0));
//	m.GetParameter("normTTbar").SetPrior(new BCGaussianPrior(1.0, 0.10));
	
//	m.GetParameter("shiftQCD").SetPrior(new BCGaussianPrior(0.0, 100.0));
//	m.GetParameter("shiftTTbar").SetPrior(new BCGaussianPrior(0.0, 100.0));
//	m.GetParameter("stretchQCD").SetPrior(new BCGaussianPrior(1.0, 0.50));
	m.GetParameter("stretchTTbar").SetPrior(new BCGaussianPrior(1.0, 0.10));
	m.GetParameter("shiftQCD").SetPriorConstant();
	m.GetParameter("shiftTTbar").SetPriorConstant();
	m.GetParameter("stretchQCD").SetPriorConstant();
//	m.GetParameter("stretchTTbar").SetPriorConstant();

//	m.GetParameter("shiftQCD").Fix(-19.8);
//	m.GetParameter("shiftTTbar").Fix(4.3);
//	m.GetParameter("stretchQCD").Fix(0.98);
//	m.GetParameter("stretchTTbar").Fix(1.02);


//	cout << m.GetParameter("shiftTTbar").GetPrior(0.0) << endl;
//	cout << m.GetParameter("stretchTTbar").GetPrior(0.0) << endl;
	
	/// Data:
//	for (unsigned i = 0; i < channels.size(); ++i) m.SetData(channels[i], *input_histograms["mavg" + channels[i] + "__DATA"]);
	
	// Normalize the posterior by integrating it over the full par. space
//	m.Normalize();

//	m.SetOptimizationMethod(BCIntegrate::kOptMetropolis);
//	m.SetOptimizationMethod(BCIntegrate::kOptMinuit);
////	TMinuit* min = m.GetMinuit();
////	min->Command("SET STR 2");		// Doesn't change anything because the TMinuit instance is remade every time FindMode is called (see FindModeMinuit).
//	m.FindMode();
//	return 1;
//	m.FindMode(vector<double> {0.845, 1.15, -19.8, 4.3, 0.98, 1.03});
//	m.FindMode(vector<double> {1.0, 1.0, 0.0, 0.0, 1.0, 1.0});
//	m.FindMode(vector<double> {0.8, 1.0, -20, 5, 1.0, 1.0});
//	m.FindMode(vector<double> {0.845474, 1.14825, -19.7312, 4.2893, 0.98243, 1.02859});
//	m.FindMode(vector<double> {0.845474, 1.14825, -19.7312, 4.28935, 0.98243, 1.02859});
	
	
	
//	m.SetNIterationsPreRunMin(500000);		// Set iterations (default is 100000)
	m.SetNIterationsPreRunMax(50000);		// Set iterations (default is 100000)
	m.SetNIterationsRun(50000);		// Set iterations (default is 100000)
//	m.SetNChains(8);		// Seg. faults.
//	m.SetPrecision(BCEngineMCMC::kHigh);
	m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
//	m.SetMarginalizationMethod(BCIntegrate::kMargGrid);
//	m.SetMarginalizationMethod(BCIntegrate::kMargMonteCarlo);
//	m.SetProposeMultivariate(false);
	
	
//	m.SetInitialPositionScheme(BCEngineMCMC::kInitUserDefined);
//	m.SetInitialPositions(vector<double> {0.8, 1.0, -20.0, 5.0, 1.0, 1.0});		// ERROR: "Too few initial positions provided"
	
	// WORKS:
//	m.SetInitialPositionScheme(BCEngineMCMC::kInitUserDefined);		// Works, but it's redundant with SetInitialPositions.
//	m.SetInitialPositions(vector<vector<double>> {{0.8, 1.0, -20.0, 5.0, 1.0, 1.0}, {0.8, 1.0, -20.0, 5.0, 1.0, 1.0}, {0.8, 1.0, -20.0, 5.0, 1.0, 1.0}, {0.8, 1.0, -20.0, 5.0, 1.0, 1.0}});
//	vector<vector<double>> initial_mcmc;
//	for (unsigned i = 0; i < m.GetNChains(); ++i) initial_mcmc.push_back(m.GetBestFitParameters());
//	m.SetInitialPositions(initial_mcmc);
	
//	return 1;
//	m.MarginalizeAll();
//	m.get_histogram("SR", "DATA");
//	cout << m.GetMarginalized("normMs300").GetQuantile(0.95) << endl;
	
//	TH1* h = m.GetMarginalizedHistogram("normMs300");
//	TCanvas* tc = new TCanvas();
//	h->Draw();
//	tc->SaveAs("test.pdf");
	
//	for (int i = 0; i < 5; ++i) cout << m.make_toy("SR")->GetMean() << endl;
	
	vector<vector<double>> limits;
//	limits.push_back(m.get_expected_limit("Ms100", 4));
//	limits.push_back(m.get_observed_limit("Ms150"));
//	limits.push_back(m.get_observed_limit("Ms200"));
//	limits.push_back(m.get_observed_limit("Ms250"));
	limits.push_back(m.get_observed_limit("Ms300", 10));
//	limits.push_back(m.get_observed_limit("Ms400"));
//	limits.push_back(m.get_observed_limit("Ms500"));
	
	for (unsigned i = 0; i < limits[0].size(); ++i) cout << limits[0][i] << endl;
	
//	return 1;

//	// run mode finding; by default using Minuit
////	m.FindMode(m.GetBestFitParameters());

//	// draw all marginalized distributions into a PDF file
	m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");		// Commenting this out gives the weird error at the end.

//	// BROKEN FOR NOW:
////	for (int i = 0; i < m.GetNChannels(); ++i) {
////		BCMTFChannel* channel = m.GetChannel(i);
////		channel->PrintTemplates(channel->GetSafeName() + "_templates.pdf");
////		m.PrintStack(i, m.GetBestFitParameters(), channel->GetSafeName() + "_stack.pdf");
////	}
//	
//	
////	unsigned tote = m.GetBestFitParameterErrors().size();
////	cout << tote << endl;
////	for (unsigned i = 0; i < tote; ++i) {
////		cout << m.GetBestFitParameterErrors()[i] << endl;
////	}

//	// print summary plots
//	// m.PrintParameterPlot(m.GetSafeName() + "_parameters.pdf");
//	// m.PrintCorrelationPlot(m.GetSafeName() + "_correlation.pdf");
//	// m.PrintCorrelationMatrix(m.GetSafeName() + "_correlationMatrix.pdf");
//	// m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_update.pdf");

//	// print results of the analysis into a text file
//	m.PrintSummary();

//	// close log file
//	BCLog::OutSummary("Exiting");
//	BCLog::CloseLog();
	
	return 0;
}
