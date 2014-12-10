//////////////////////////////////////////////////////////////////////////
// 
// This code:
//	1) Loads a histogram from a ROOT file
//	2) Convolutes it with a Gaussian
//	3) Creates a PDF with the convoluted histogram
//	4) Generates data from PDF
//	5) Fits generated data to PDF
//	6) Plots results, including Chi^2 of fit
//
// 10/2014 - Chandler Schlupf
// 
/////////////////////////////////////////////////////////////////////////

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
using namespace RooFit;

void DBeta_roof_1DHistPDF_fit()
{
	//Load Data
	TFile *F_0NuBB = new TFile("output/DBeta_Debug_0NuBB.root","READ");
//	TFile *F_0NuBB = new TFile("/u/nobackup/lwinslow/elagin/data/db_out/Te130_0vbb_1e6.root","READ"); //Contains more data

	//Load Histogram from File
	TH1 * hSingle0 = (TH1F *) F_0NuBB->Get("hSingle");
	
	//Construct Observables
	RooRealVar x("x","x",-.5,3.5);
	
	//Energy Resolution for Sigma of Gaussian
	const double max_energy = 2.53697;
	double energy_res_10 = max_energy*.1;
	
	//Widen range of histograms so that the convolution isn't messed up
	hSingle0->SetBit(TH1::kCanRebin);
	hSingle0->Fill(-.1,1);	hSingle0->Fill(2.7,1);

	//Root->Roofit histograms
	RooDataHist dhSingle0("dhSingle0","dhSingle0",x, Import(*hSingle0));
	
	//Create PDF from Histogram
	RooHistPdf pdfhSingle0("pdfhSingle0","pdfhSingle0",x,dhSingle0,0);
	
	//Construct Gaussian
	RooRealVar meang("meang","meang",0);
	RooRealVar sigmag_10("sigmag_10","sigmag_10",energy_res_10);
	RooGaussian gauss_10("gauss_10","gauss_10",x,meang,sigmag_10);
	
	//Construct PDF (x) Gauss	
	x.setBins(10000,"cache"); 
	RooFFTConvPdf cSingle0_10("cSingle0_10","Single0 data (X) gauss_10",x,pdfhSingle0,gauss_10);

	Double_t total_events = 20000;
	RooRealVar nsig("nsig","number of signal events",10.,0.,total_events*2);
	RooExtendPdf ecSingle0_10("ecSingle0_10","Extended cSingle0_10",cSingle0_10,nsig);	//Extend the PDF so that the number of events can be fit
	
	//Get Sample Events from Convoluted PDF
	RooDataSet * ecSingle0_data_10 = ecSingle0_10.generate(x,total_events, Extended(kTRUE));

//	//Fit Sampled Data to Convoluted PDF
	RooFitResult * fitresults = ecSingle0_10.fitTo(*ecSingle0_data_10, Save(),Extended());
	cout<<"Number of Fit Events: "<<nsig.getVal()<<" +/- "<<nsig.getError()<<endl;;

//	//Plot
	RooPlot* frame = x.frame(0,3);
//	pcSingle0_10.plotOn(frame);
	frame->SetName("E1 0vBB PDF Fit to Sampled Convoluted Data");
	frame->SetTitle("E1 0vBB PDF Fit to Sampled Convoluted Data");
	frame->SetXTitle("E1 [MeV]");
	ecSingle0_10.plotOn(frame, LineColor(kRed),Name("pdf"),Normalization(1.0,RooAbsReal::RelativeExpected)) ;
	ecSingle0_data_10->plotOn(frame, MarkerColor(kYellow),Name("data"));
	frame->Draw();
	
	//Add Chi^2 of Fit to Plot
	double chi2 = frame->chiSquare();
	TPaveText* tbox = new TPaveText(.3,.15,.6,.3,"blNDC");
	tbox->SetBorderSize(0);
	tbox->SetFillStyle(0);
	tbox->AddText(Form("Chi Square = %g",chi2));
	tbox->AddText(Form("Event # Fit = %G +/- %G",nsig.getVal(), nsig.getError()));
	tbox->AddText(Form("Event # Generated = %G",total_events));
	tbox->SetTextSize(.04);
	frame->addObject(tbox) ;
	
	//Add Legend to Plot
	TLegend * leg = new TLegend(.75,.7,.88,.88);
	leg->SetFillColor(kWhite);
	leg->SetLineColor(kWhite);
	leg->AddEntry(frame->findObject("data"),"data","p");
	leg->AddEntry(frame->findObject("pdf"),"PDF","l");
	leg->Draw("same");
	frame->Draw("same");

};


