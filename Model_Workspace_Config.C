
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooDataHist.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TH1F.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
using namespace RooFit;
using namespace RooStats;

double Model_Workspace_Config()
{
	
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING); //silence fitting messages
	
	//Load Data
	TFile *F_0NuBB = new TFile("output/DBeta_Debug_0NuBB.root","READ");
	TFile *F_2NuBB = new TFile("output/DBeta_Debug_2NuBB.root","READ");
//	TFile *F_0NuBB = new TFile("/u/nobackup/lwinslow/elagin/data/db_out/Te130_0vbb_1e6.root","READ"); //Contains more data
//	TFile *F_2NuBB = new TFile("/u/nobackup/lwinslow/elagin/data/db_out/Te130_2vbb_1e6.root","READ"); //Contains more data
	
	//Load histogram from file
	TH1 * hSingle0 = (TH1F *) F_0NuBB->Get("hSingle");
	TH1 * hSingle2 = (TH1F *) F_2NuBB->Get("hSingle");
	TH1 * hSum0 = (TH1 *) F_0NuBB->Get("hSum");
	TH1 * hSum2 = (TH1 *) F_2NuBB->Get("hSum");
	TH1F * hAngle0 = (TH1F *) F_0NuBB->Get("hAngle");	
	TH1F * hAngle2 = (TH1F *) F_2NuBB->Get("hAngle");
	
	//Construct observables
	RooRealVar x_Single("x_Single","x_Single",0,3.5);
	RooRealVar x_Sum("x_Sum","x_Sum",0,3.5);
	RooRealVar x_Ang("x_Ang","x_Ang",-1.,1.);
	
	//Energy resolution for sigma of gaussian
	const double max_energy = 2.53697;
	double energy_res_10 = max_energy*0.1;
	
	//Root->Roofit histograms
	RooDataHist dhSingle0("dhSingle0","dhSingle0",x_Single, Import(*hSingle0));
	RooDataHist dhSingle2("dhSingle2","dhSingle2",x_Single, Import(*hSingle2));
	RooDataHist dhSum0("dhSum0","dhSum0",x_Sum, Import(*hSum0));
	RooDataHist dhSum2("dhSum2","dhSum2",x_Sum, Import(*hSum2));
	RooDataHist dhAngle0("dhAngle0","dhAngle0",x_Ang,Import(*hAngle0));
	RooDataHist dhAngle2("dhAngle2","dhAngle2",x_Ang,Import(*hAngle2));
	
	
	//Create PDF from histogram
	RooHistPdf pdfhSingle0("pdfhSingle0","pdfhSingle0",x_Single,dhSingle0,0);
	RooHistPdf pdfhSingle2("pdfhSingle2","pdfhSingle2",x_Single,dhSingle2,0);
	RooHistPdf pdfhSum0("pdfhSum0","pdfhSum0",x_Sum,dhSum0,0);
	RooHistPdf pdfhSum2("pdfhSum2","pdfhSum2",x_Sum,dhSum2,0);
	RooHistPdf pdfhAngle0("pdfhAngle0","pdfhAngle0",x_Ang,dhAngle0,0);
	RooHistPdf pdfhAngle2("pdfhAngle2","pdfhAngle2",x_Ang,dhAngle2,0);
	
	//Construct gaussian
	RooRealVar meang("meang","meang",0);
	RooRealVar sigmag_10("sigmag_10","sigmag_10",energy_res_10);
	RooGaussian gauss_sin_10("gauss_sin_10","gauss_sin_10",x_Single,meang,sigmag_10);
	RooGaussian gauss_sum_10("gauss_sum_10","gauss_sum_10",x_Sum,meang,sigmag_10);
	RooGaussian gauss_ang_10("gauss_ang_10","gauss_ang_10",x_Ang,meang,sigmag_10);
	
	//Construct PDF (x) Gauss	
	x_Single.setBins(10000,"cache"); 
	x_Sum.setBins(10000,"cache"); 
	x_Ang.setBins(10000,"cache");
	RooFFTConvPdf cSingle0_10("cSingle0_10","Single0 data (X) gauss_10",x_Single,pdfhSingle0,gauss_sin_10);
	RooFFTConvPdf cSingle2_10("cSingle2_10","Single2 data (X) gauss_10",x_Single,pdfhSingle2,gauss_sin_10);
	RooFFTConvPdf cSum0_10("cSum0_10","Sum0 data (X) gauss_10",x_Sum,pdfhSum0,gauss_sum_10);
	RooFFTConvPdf cSum2_10("cSum2_10","Sum2 data (X) gauss_10",x_Sum,pdfhSum2,gauss_sum_10);
	RooFFTConvPdf cAngle0_10("cAngle0_10","Angle0 data (X) gauss_ang_10",x_Ang,pdfhAngle0,gauss_ang_10);
	RooFFTConvPdf cAngle2_10("cAngle2_10","Angle2 data (X) gauss_ang_10",x_Ang,pdfhAngle2,gauss_ang_10);

	//Set Buffer Fraction: FFT assumes a cyclic variable.  Single and Sum energy are not cyclic, so a buffer zone is provided at each end of 30%.
	double buff_frac = 1;
	cSingle0_10.setBufferFraction(buff_frac);
	cSingle2_10.setBufferFraction(buff_frac);
	cSum0_10.setBufferFraction(buff_frac);
	cSum2_10.setBufferFraction(buff_frac);
	
	//Extend the PDF so that the number of events can be fit
	Double_t num_generated = 400.; //Data size generated from PDFs to be fit by PDFs
	Double_t n_0_fraction = 0.;
	Double_t n_0_generated = num_generated*n_0_fraction;
	Double_t n_2_generated = num_generated-n_0_generated;
	n_2_generated = num_generated;
	Double_t n_init = num_generated; //n_0 and n_2 initial value (not important because this is a variable in the fit.)
	RooRealVar n_0("n_0","number of 0vBB events",n_0_generated, -num_generated*10.,num_generated*10); //n_0 must have a negative lower bound in order to get an accurate q_0 distribution
	RooRealVar n_2("n_2","number of 2vBB events",n_2_generated, -num_generated*10.,num_generated*10);
	RooExtendPdf ecSingle0_10("ecSingle0_10","Extended cSingle0_10",cSingle0_10,n_0);
	RooExtendPdf ecSingle2_10("ecSingle2_10","Extended cSingle2_10",cSingle2_10,n_2);	
	RooExtendPdf ecSum0_10("ecSum0_10","Extended cSum0_10",cSum0_10,n_0);
	RooExtendPdf ecSum2_10("ecSum2_10","Extended cSum2_10",cSum2_10,n_2);
	RooExtendPdf ecAngle0_10("ecAngle0_10","Extended cAngle0_10",cAngle0_10,n_0);
	RooExtendPdf ecAngle2_10("ecAngle2_10","Extended cAngle2_10",cAngle2_10,n_2);
	
	//Create 0vBB + 2vBB Pdfs
	RooAddPdf pdf_Single("pdf_Single","Single 0vBB + 2vBB", RooArgList(ecSingle0_10,ecSingle2_10));
	RooAddPdf pdf_Sum("pdf_Sum","Sum 0vBB + 2vBB", RooArgList(ecSum0_10,ecSum2_10));
	RooAddPdf pdf_Angle("pdf_Angle","Angle 0vBB + 2vBB", RooArgList(ecAngle0_10,ecAngle2_10));
	
//	RooRandom::randomGenerator()->SetSeed(0);
	
	//Get sample events from convoluted sig+back PDFs
	n_0.setVal(n_0_generated);
	n_2.setVal(n_2_generated);
	RooDataSet * ecSingle_data_10 = pdf_Single.generate(x_Single,num_generated,Extended(kTRUE));
	RooDataSet * ecSum_data_10 = pdf_Sum.generate(x_Sum,num_generated,Extended(kTRUE));
	RooDataSet * ecAngle_data_10 = pdf_Angle.generate(x_Ang,num_generated,Extended(kTRUE));

	//Define category to distinguish observables
	RooCategory data_type("data_type","data_type");
	data_type.defineType("Single");
	data_type.defineType("Sum");
	data_type.defineType("Angle");
	
	//Construct combined dataset 
	RooDataSet data_all("data_all","all data",RooArgSet(x_Single, x_Sum, x_Ang), Index(data_type),
			Import("Single",*ecSingle_data_10),
			Import("Sum",*ecSum_data_10),
			Import("Angle",*ecAngle_data_10));
	
	//Construct a simultaneous pdf
//	n_0.setVal(n_0_generated);
//	n_2.setVal(n_2_generated);
	RooSimultaneous pdf_all("pdf_all","simultaneous pdf",data_type);
	pdf_all.addPdf(pdf_Single,"Single");
	pdf_all.addPdf(pdf_Sum,"Sum");
	pdf_all.addPdf(pdf_Angle,"Angle");		
	RooSimultaneous pdf_n2_all("pdf_n2_all","simultaneous n2 pdf",data_type);
	pdf_n2_all.addPdf(ecSingle2_10, "Single");
	pdf_n2_all.addPdf(ecSum2_10, "Sum");
	pdf_n2_all.addPdf(ecAngle2_10, "Angle");
	
//	RooDataSet *data_all = pdf_all.generate(RooArgSet(x_Single, x_Sum, x_Ang, data_type), Extended(kTRUE));

	RooWorkspace w("w");

	
	ModelConfig model("model", &w);
	model.SetPdf(pdf_Single);
	model.SetParametersOfInterest(n_0);
	model.SetObservables(RooArgSet(x_Single));
	model.SetNuisanceParameters(n_2);
    w.import(model);
//	w.import(data_type);
//	w.import(pdf_all);
	w.import(*ecSingle_data_10);
	w.Print("v");
    w.writeToFile("model_workspace_2.root",true);

//    model.Print("v");
//    w.Print("v");


};












