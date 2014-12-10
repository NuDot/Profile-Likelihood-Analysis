//This code simulatneously fits E_Single, E_Sum and Angle data usin g their expected PDFs.  It then finds the profile likelihood ratio and cooresponding test statistic (q0).
//In lines 89 and 90 you can change the number of n_2 and n_0 events to produce.
//It also plots all of the fits. (It's exactly the same as Full_Fit_and_Test_Stat.C except it produces plots too)



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
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
using namespace RooFit;
using namespace RooStats;

double Full_Fit_and_Test_Stat_wPlots()
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
	RooRealVar x("x","x",0,3.5);
	RooRealVar x_ang("x_ang","x_ang",-1.,1.);
	
	//Energy resolution for sigma of gaussian
	const double max_energy = 2.53697;
	double energy_res_10 = max_energy*0.1;
	
	//Root->Roofit histograms
	RooDataHist dhSingle0("dhSingle0","dhSingle0",x, Import(*hSingle0));
	RooDataHist dhSingle2("dhSingle2","dhSingle2",x, Import(*hSingle2));
	RooDataHist dhSum0("dhSum0","dhSum0",x, Import(*hSum0));
	RooDataHist dhSum2("dhSum2","dhSum2",x, Import(*hSum2));
	RooDataHist dhAngle0("dhAngle0","dhAngle0",x_ang,Import(*hAngle0));
	RooDataHist dhAngle2("dhAngle2","dhAngle2",x_ang,Import(*hAngle2));
	
	//Create PDF from histogram
	RooHistPdf pdfhSingle0("pdfhSingle0","pdfhSingle0",x,dhSingle0,0);
	RooHistPdf pdfhSingle2("pdfhSingle2","pdfhSingle2",x,dhSingle2,0);
	RooHistPdf pdfhSum0("pdfhSum0","pdfhSum0",x,dhSum0,0);
	RooHistPdf pdfhSum2("pdfhSum2","pdfhSum2",x,dhSum2,0);
	RooHistPdf pdfhAngle0("pdfhAngle0","pdfhAngle0",x_ang,dhAngle0,0);
	RooHistPdf pdfhAngle2("pdfhAngle2","pdfhAngle2",x_ang,dhAngle2,0);
	
	//Construct gaussian
	RooRealVar meang("meang","meang",0);
	RooRealVar sigmag_10("sigmag_10","sigmag_10",energy_res_10);
	RooGaussian gauss_10("gauss_10","gauss_10",x,meang,sigmag_10);
	RooGaussian gauss_ang_10("gauss_ang_10","gauss_ang_10",x_ang,meang,sigmag_10);
	
	//Construct PDF (x) Gauss	
	x.setBins(10000,"cache"); 
	x_ang.setBins(10000,"cache");
	RooFFTConvPdf cSingle0_10("cSingle0_10","Single0 data (X) gauss_10",x,pdfhSingle0,gauss_10);
	RooFFTConvPdf cSingle2_10("cSingle2_10","Single2 data (X) gauss_10",x,pdfhSingle2,gauss_10);
	RooFFTConvPdf cSum0_10("cSum0_10","Sum0 data (X) gauss_10",x,pdfhSum0,gauss_10);
	RooFFTConvPdf cSum2_10("cSum2_10","Sum2 data (X) gauss_10",x,pdfhSum2,gauss_10);
	RooFFTConvPdf cAngle0_10("cAngle0_10","Angle0 data (X) gauss_ang_10",x_ang,pdfhAngle0,gauss_ang_10);
	RooFFTConvPdf cAngle2_10("cAngle2_10","Angle2 data (X) gauss_ang_10",x_ang,pdfhAngle2,gauss_ang_10);
	
	//Set Buffer Fraction: FFT assumes a cyclic variable.  Single and Sum energy are not cyclic, so a buffer zone is provided at each end of 30%.
	double buff_frac = 1;
	cSingle0_10.setBufferFraction(buff_frac);
	cSingle2_10.setBufferFraction(buff_frac);
	cSum0_10.setBufferFraction(buff_frac);
	cSum2_10.setBufferFraction(buff_frac);
	
	//Extend the PDF so that the number of events can be fit
	double num_generated = 1000.; //Data size generated from PDFs to be fit by PDFs
	double n_0_fraction = 0.;
	double n_0_generated = num_generated*n_0_fraction;
	double n_2_generated = num_generated-n_0_generated;
	double n_init = num_generated; //n_0 and n_2 initial value (not important because this is a variable in the fit.)
	RooRealVar n_0("n_0","number of 0vBB events",n_0_generated, -num_generated*10,num_generated*10);
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
	
	//Get sample events from convoluted sig+back PDFs
	n_0.setVal(n_0_generated);
	n_2.setVal(n_2_generated);
	RooDataSet * ecSingle_data_10 = pdf_Single.generate(x,num_generated,Extended(kTRUE));
	RooDataSet * ecSum_data_10 = pdf_Sum.generate(x,num_generated,Extended(kTRUE));
	RooDataSet * ecAngle_data_10 = pdf_Angle.generate(x_ang,num_generated,Extended(kTRUE));
	
	//Define category to distinguish observables
	RooCategory data_type("data_type","data_type");
	data_type.defineType("Single");
	data_type.defineType("Sum");
	data_type.defineType("Angle");

	//Construct combined dataset 
	RooDataSet data_all("data_all","all data",RooArgSet(x,x_ang),Index(data_type),
			Import("Single",*ecSingle_data_10),
			Import("Sum",*ecSum_data_10),
			Import("Angle",*ecAngle_data_10));
	
	
	//Construct a simultaneous pdf
	n_0.setVal(n_0_generated);
	n_2.setVal(n_2_generated);
	RooSimultaneous pdf_all("pdf_all","simultaneous pdf",data_type);
	pdf_all.addPdf(pdf_Single,"Single");
	pdf_all.addPdf(pdf_Sum,"Sum");
	pdf_all.addPdf(pdf_Angle,"Angle");		
	
	//Null Hypothesis L: n_0 = 0, n_2 = float
	n_0.setVal(0);
	n_0.setConstant(kTRUE);
	RooFitResult * fitresults = pdf_all.fitTo(data_all, Save(), Extended(), PrintLevel(-1000));
	cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	cout<<"Fit Results  n_0=0 FIXED:"<<endl;
	cout<<"Number of 0vBB Fit Events: "<<n_0.getVal()<<" +/- "<<n_0.getError()<<endl;
	cout<<"Number of 2vBB Fit Events: "<<n_2.getVal()<<" +/- "<<n_2.getError()<<endl;
	cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	Float_t min_log_fixed = fitresults->minNll();
	
	//Plot First Fits
	TCanvas *c_Single = new TCanvas("c_Single","Single");
	RooPlot* frame_Single = x.frame(0,3);
	frame_Single->SetName("E1 PDF Fit to Sampled Convoluted Data");
	frame_Single->SetTitle("E1 PDF Fit to Sampled Convoluted Data");
	frame_Single->SetXTitle("E1 [MeV]");
	ecSingle_data_10->plotOn(frame_Single, MarkerColor(kBlack), Name("data"));
	pdf_Single.plotOn(frame_Single, LineColor(kBlue),Name("pdf"),Normalization(1.0,RooAbsReal::RelativeExpected));
	frame_Single->Draw();
	
	TCanvas *c_Sum = new TCanvas("c_Sum","Sum");
	RooPlot* frame_Sum = x.frame(0,3);
	frame_Sum->SetName("E1+E2 PDF Fit to Sampled Convoluted Data");
	frame_Sum->SetTitle("E1+E2 PDF Fit to Sampled Convoluted Data");
	frame_Sum->SetXTitle("E1+E2 [MeV]");
	ecSum_data_10->plotOn(frame_Sum, MarkerColor(kBlack), Name("data"));
	pdf_Sum.plotOn(frame_Sum, LineColor(kBlue),Name("pdf"),Normalization(1.0,RooAbsReal::RelativeExpected));
	frame_Sum->Draw();

	TCanvas *c_Angle = new TCanvas("c_Angle","Angle");
	RooPlot* frame_Angle = x_ang.frame(-1,1);
	frame_Angle->SetName("Cos(Theta) PDF Fit to Sampled Convoluted Data");
	frame_Angle->SetTitle("Cos(Theta) PDF Fit to Sampled Convoluted Data");
	frame_Angle->SetXTitle("Cos(Theta)");
	ecAngle_data_10->plotOn(frame_Angle, MarkerColor(kBlack), Name("data"));
	pdf_Angle.plotOn(frame_Angle, LineColor(kBlue),Name("pdf"),Normalization(1.0,RooAbsReal::RelativeExpected));
	frame_Angle->Draw();
	
	//Full Model L: n_0, n_2 = float
	n_0.setVal(n_0_generated);
	n_0.setConstant(kFALSE);
	RooFitResult * fitresults_n_0_0 = pdf_all.fitTo(data_all, Save(), Extended(), PrintLevel(-1000));
	cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	cout<<"Fit Results:"<<endl;
	cout<<"Number of 0vBB Fit Events: "<<n_0.getVal()<<" +/- "<<n_0.getError()<<endl;
	cout<<"Number of 2vBB Fit Events: "<<n_2.getVal()<<" +/- "<<n_2.getError()<<endl;
	cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	Float_t min_log_float = fitresults_n_0_0->minNll();
	
	//Plotting Second Fit
	pdf_Single.plotOn(frame_Single, LineColor(kGreen),Name("pdf_n0"),Normalization(1.0,RooAbsReal::RelativeExpected));
	pdf_Single.plotOn(frame_Single, Components(ecSingle0_10), LineStyle(kDashed), LineColor(kRed), Name("0vBB pdf"), Normalization(1.0,RooAbsReal::RelativeExpected));
	frame_Single->Draw();
	pdf_Sum.plotOn(frame_Sum, LineColor(kGreen),Name("pdf_n0"),Normalization(1.0,RooAbsReal::RelativeExpected));
	pdf_Sum.plotOn(frame_Sum, Components(ecSum0_10), LineStyle(kDashed), LineColor(kRed), Name("0vBB pdf"), Normalization(1.0,RooAbsReal::RelativeExpected));
	frame_Sum->Draw();
	pdf_Angle.plotOn(frame_Angle, LineColor(kGreen),Name("pdf_n0"),Normalization(1.0,RooAbsReal::RelativeExpected));
	pdf_Angle.plotOn(frame_Angle, Components(ecAngle0_10), LineStyle(kDashed), LineColor(kRed), Name("0vBB pdf"), Normalization(1.0,RooAbsReal::RelativeExpected));
	frame_Angle->Draw();
	
	//Create Legend for All Plots
	TLegend * leg = new TLegend(.6,.7,.88,.88);
	leg->SetFillColor(kWhite);
	leg->SetLineColor(kWhite);
	leg->SetTextSize(.03);
	leg->AddEntry(frame_Single->findObject("pdf"),"Total PDF, n_0=0 fit","l");
	leg->AddEntry(frame_Single->findObject("pdf_n0"),"Total PDF, n_0=float fit","l");
	leg->AddEntry(frame_Single->findObject("0vBB pdf"),"0vBB PDF, n_0=float fit","l");
	
	Float_t q0 = 2*(min_log_fixed - min_log_float);
	cout<<"q0: "<<q0<<endl;
	cout<<"Z: "<<sqrt(q0)<<endl;

	TPaveText* tbox_Single = new TPaveText(.5,.6,.88,.7,"blNDC");
	tbox_Single->SetTextAlign(31);
	tbox_Single->SetBorderSize(0);
	tbox_Single->SetFillStyle(0);
	tbox_Single->AddText(Form("Z = %.2g",sqrt(q0)));
	tbox_Single->SetTextSize(.03);
	
	//Draw Frames
	c_Single->cd();
	leg->Draw("same");
	tbox_Single->Draw("same");
	frame_Single->Draw("same");
	
	c_Sum->cd();
	leg->Draw("same");
	tbox_Single->Draw("same");
	frame_Sum->Draw("same");
	
	c_Angle->cd();
	leg->Draw("same");
	tbox_Single->Draw("same");
	frame_Angle->Draw("same");
	
	return q0;
};












