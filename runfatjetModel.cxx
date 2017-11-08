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

#include "fatjetModel.h"

#include <iostream>
#include <TFile.h>

using namespace std;

TString input_file = "theta_plots_sim_sig_sb.root";
vector<TString> regions = {"SR", "SB"};
vector<TString> processes = {"DATA", "QCD", "TTbar", "Ms100", "Ms150", "Ms200", "Ms250", "Ms300", "Ms400", "Ms500"};
vector<TString> processes_bkg = {"QCD", "TTbar"};
vector<TString> processes_sig = {"Ms100", "Ms150", "Ms200", "Ms250", "Ms300", "Ms400", "Ms500"};

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
	BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);
	
	// create new fatjetModel object
	fatjetModel m("fatjetModel");
	
	// set precision
	m.SetPrecision(BCEngineMCMC::kMedium);
	
	// Fetch input histograms:
	map<TString,TH1D*> input_histograms = fetch_input_histograms(input_file);
	/// Print out information about them:
	unsigned ninput = input_histograms.size();
	if (ninput > 0) {
		BCLog::OutSummary(Form("Fetched %u input histograms:", ninput));
		for (auto const& i : input_histograms) {
//			BCLog::OutSummary(Form("  * %s: %s", i.first.Data(), i.second->GetName()));
			BCLog::OutSummary(Form("  * %s", i.first.Data()));
		}
		BCLog::OutSummary("");
	}
	else {
		BCLog::OutError("No input histograms were found!");
		return 0;
	}
	
	// Define everything:
	/// Channels:
//	m.AddChannel("SR");
	m.AddChannel("SB");
	
	/// Processes:
	//// Note: this order governs the stack order.
//	m.AddProcess("normQCDSR", 0.0, 50.0);
	m.AddProcess("normQCDSB", 0.0, 50.0);
	m.AddProcess("normTTbar", 0.0, 50.0);
//	m->AddProcess("signal",                0., 400.);

//	m.SetTemplate("SR", "normTTbar", *input_histograms["mavgSR__TTbar"], 1.0, input_histograms["mavgSR__TTbar"]->Integral());
	m.SetTemplate("SB", "normTTbar", *input_histograms["mavgSB__TTbar"], 1.0, input_histograms["mavgSB__TTbar"]->Integral());
//	m.SetTemplate("SR", "normQCDSR", *input_histograms["mavgSR__QCD"], 1.0, input_histograms["mavgSR__QCD"]->Integral());
	m.SetTemplate("SB", "normQCDSB", *input_histograms["mavgSB__QCD"], 1.0, input_histograms["mavgSB__QCD"]->Integral());
	
	m.fill_hists();
	m.fill_cdfs();
	
	m.AddParameter("shiftQCD", -100.0, 100.0);
	m.AddParameter("shiftTTbar", -100.0, 100.0);
	m.AddParameter("stretchQCD", 0.0, 5.0);
	m.AddParameter("stretchTTbar", 0.0, 5.0);
	
	
//	m.GetParameter("normQCDSR").SetPriorConstant();
	m.GetParameter("normQCDSB").SetPriorConstant();
	m.GetParameter("normTTbar").SetPriorConstant();
//	m.GetParameter("normQCDSR").SetPrior(new BCGaussianPrior(1.0, 100.0));
//	m.GetParameter("normQCDSB").SetPrior(new BCGaussianPrior(1.0, 100.0));
//	m.GetParameter("normTTbar").SetPrior(new BCGaussianPrior(1.0, 0.20));
	
//	m.GetParameter("shiftQCD").SetPrior(new BCGaussianPrior(0.0, 100.0));
//	m.GetParameter("shiftTTbar").SetPrior(new BCGaussianPrior(0.0, 100.0));
//	m.GetParameter("stretchQCD").SetPrior(new BCGaussianPrior(1.0, 0.50));
//	m.GetParameter("stretchTTbar").SetPrior(new BCGaussianPrior(1.0, 0.50));
	m.GetParameter("shiftQCD").SetPriorConstant();
	m.GetParameter("shiftTTbar").SetPriorConstant();
	m.GetParameter("stretchQCD").SetPriorConstant();
	m.GetParameter("stretchTTbar").SetPriorConstant();

//	m.GetParameter("shiftQCD").Fix(-19.8);
//	m.GetParameter("shiftTTbar").Fix(3.4);
//	m.GetParameter("stretchQCD").Fix(0.98);
//	m.GetParameter("stretchTTbar").Fix(1.008);
	
	/// Data:
//	m.SetData("SR", *input_histograms["mavgSR__DATA"]);
	m.SetData("SB", *input_histograms["mavgSB__DATA"]);
	
	// Normalize the posterior by integrating it over the full par. space
//	m.Normalize();

	// run MCMC, marginalizing posterior
//	m.SetNIterationsPreRunMin(500000);		// Set iterations (default is 100000)
	m.SetNIterationsPreRunMax(500000);		// Set iterations (default is 100000)
	m.SetNIterationsRun(500000);		// Set iterations (default is 100000)
	
//	m.SetOptimizationMethod(BCIntegrate::kOptMetropolis);
	m.SetOptimizationMethod(BCIntegrate::kOptMinuit);
	m.FindMode();
	return 0;
	
//	m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
//	m.MarginalizeAll();

	// run mode finding; by default using Minuit
//	m.FindMode(m.GetBestFitParameters());

	// draw all marginalized distributions into a PDF file
	m.PrintAllMarginalized(m.GetSafeName() + "_plots.pdf");

	for (int i = 0; i < m.GetNChannels(); ++i) {
		BCMTFChannel* channel = m.GetChannel(i);
		channel->PrintTemplates(channel->GetSafeName() + "_templates.pdf");
		m.PrintStack(i, m.GetBestFitParameters(), channel->GetSafeName() + "_stack.pdf");
	}
	
//	unsigned tote = m.GetBestFitParameterErrors().size();
//	cout << tote << endl;
//	for (unsigned i = 0; i < tote; ++i) {
//		cout << m.GetBestFitParameterErrors()[i] << endl;
//	}

	// print summary plots
	// m.PrintParameterPlot(m.GetSafeName() + "_parameters.pdf");
	// m.PrintCorrelationPlot(m.GetSafeName() + "_correlation.pdf");
	// m.PrintCorrelationMatrix(m.GetSafeName() + "_correlationMatrix.pdf");
	// m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_update.pdf");

	// print results of the analysis into a text file
	m.PrintSummary();

	// close log file
	BCLog::OutSummary("Exiting");
	BCLog::CloseLog();
	
	return 0;
}
