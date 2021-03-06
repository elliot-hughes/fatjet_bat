#include <BAT/BCMath.h>
#include "fatjetModel.h"
#include "fatjetModel.cxx"
#include <Deracination/Straphanger/test/decortication/macros/common.cc>

//TString input_file = "theta_plots_sim_sig_sb.root";
TString input_file = "theta_plots_sim_sig_sbb.root";
TString input_file_cdf = "analysis_plots_sbb_f1.root";
//TString input_file_cdf = "analysis_plots_sig_f1.root";
int cnum = 0;

double nll(TH1D* h_data, TH1D* h_qcd, TH1D* h_ttbar) {
	double result = 0;
	int nbins = h_data->GetNbinsX();
	int nbins_flagged = 0;
	for (int ibin = 1; ibin <= nbins; ++ibin) {
		double m = h_data->GetXaxis()->GetBinUpEdge(ibin);
		if (m > 900) break;
		double obs = h_data->GetBinContent(ibin);
		double exp = h_qcd->GetBinContent(ibin) + h_ttbar->GetBinContent(ibin);
		if (exp <= 0) {
			nbins_flagged += 1;
			exp = 0.000001;
		}
//		std::cout << "obs = " << obs << ", exp = " << exp << ", L = " << BCMath::LogPoisson(obs, exp) << std::endl;
//		result += BCMath::LogPoisson(obs, exp);
		double result_component = obs*log(exp) - exp - TMath::LnGamma(obs + 1);
//		std::cout << "m = " << m << ", obs = " << obs << ", exp = " << exp << ", L = " << result_component << std::endl;		// DEBUG
		result += result_component;
	}
	return -result;
}


double check_bat(double normQCD, double shiftQCD, double stretchQCD, double normTTbar, double shiftTTbar, double stretchTTbar, bool draw=true) {
	TFile* tfile = TFile::Open(input_file);
	TFile* tfile_cdf = TFile::Open(input_file_cdf);
	fatjetModel m("fatjetModel");
	
	TString region = "SB";
	
	TH1D* h_data = (TH1D*) tfile->Get("mavg" + region + "__DATA");
	TH1D* h_qcd = (TH1D*) tfile->Get("mavg" + region + "__QCD");
	TH1D* h_ttbar = (TH1D*) tfile->Get("mavg" + region + "__TTbar");
	
//	TH1D* cdf_qcd = m.make_cdf(h_qcd);
//	TH1D* cdf_ttbar = m.make_cdf(h_ttbar);
	TH1D* cdf_qcd = (TH1D*) tfile_cdf->Get("cdf_jetht");
	TH1D* cdf_ttbar = (TH1D*) tfile_cdf->Get("cdf_ttbar");
	
	TH1D* h_qcd_morphed = (TH1D*) h_qcd->Clone(h_qcd->GetName() + TString("_morphed"));
	m.differentiate_cdf(cdf_qcd, h_qcd_morphed, normQCD, shiftQCD, stretchQCD);
	TH1D* h_ttbar_morphed = (TH1D*) h_ttbar->Clone(h_ttbar->GetName() + TString("_morphed"));
	m.differentiate_cdf(cdf_ttbar, h_ttbar_morphed, normTTbar, shiftTTbar, stretchTTbar);
	
	double result = nll(h_data, h_qcd_morphed, h_ttbar_morphed);
	cout << "L = " << result << endl;
	
//	// Draw QCD cdf:
//	TCanvas* tc_cdf = new TCanvas(TString("c") + to_string(cnum), TString("c") + to_string(cnum));
//	cnum += 1;
//	cdf_qcd->Draw();
	
	if (draw) {
		// Make stack:
		TCanvas* tc = new TCanvas(TString("c") + to_string(cnum), TString("c") + to_string(cnum));
		cnum += 1;
	
		gStyle->SetOptStat(0);
		h_data->SetTitle("");
		h_data->GetXaxis()->SetRangeUser(0, 900);
		h_qcd_morphed->SetFillColor(kBlue-9);
		h_ttbar_morphed->SetFillColor(kRed-9);
	
		THStack* hs = new THStack();
		hs->Add(h_qcd_morphed);
		hs->Add(h_ttbar_morphed);
	
		h_data->Draw("e");
		hs->Draw("hist same");
		h_data->Draw("e same");
	
		std::ostringstream oss;
		oss << "NLL = " << std::fixed << std::setprecision(2) << result;
		style_write(oss.str());
	}
	
	return result;
}



void crosscheck() {
////	check_bat(0.849546, -19.8253, 0.98038, 1.11772, 3.35824, 1.00835);
//	check_bat(0.845474, -19.7312, 0.98243, 1.14825, 4.29174, 1.0286);		// These are the params that study_analysis gives (MINUIT).
////	check_bat(1.00095, -30.1275, 0.919733, 1.84852e-11, -100, 4.88907);

////	check_bat(0.810436, -18.8426, 0.969619, 1.40995, 9.95084, 1.22927);		// These are the params that BAT gives (with MINUIT) (BEFORE CDF RES CHANGE)
//	check_bat(0.8126, -17.67, 0.9681, 1.381, 8.859, 1.2324);		// These are the params that BAT gives for MCMC.

//	// SR (BAT)
//	check_bat(0.755, -10.5, 0.987, 1.09, -13.6, 0.845);		// Float
//	check_bat(0.763, -11.5, 0.991, 1.09, -10.6, 0.882);		// Gaussian strTTbar
	
	// SBB (BAT)
	check_bat(0.501, -17.9, 1.131, 0.99, 3.2, 1.050);		// Float
	check_bat(0.501, -15.5, 1.127, 1.00, 2.8, 1.051);		// Gaussian strTTbar
	
}
