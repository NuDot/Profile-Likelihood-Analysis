//This file is an altered version of the RooStat tutorial "SimpleHypoTest.C".  
//It takes the signal/background models from files created in "Model_Workspace_Config.C".
//Still unsure about how many events it is using when testing the test statistic.  
//At the end I added a section that does the following: 
	//Finds median of distribution : integrates null distribution from median of alt distribution to infinity to determine p value, then uses p->Z calc to find expected significance. 


#include "TMath.h"
using namespace RooStats;
using namespace RooFit;

void SimpleHypoTest( const char* infile =  "model_workspace_2.root", 
                     const char* workspaceName = "w",
                     const char* modelConfigName = "model",
                     const char* dataName = "pdf_SingleData" )
{

  /////////////////////////////////////////////////////////////
  // First part is just to access the workspace file 
  ////////////////////////////////////////////////////////////
	
RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING); //silence fitting messages
	
  // open input file 
  TFile *file = TFile::Open(infile);
  if (!file) return;

  // get the workspace out of the file
  RooWorkspace* w = (RooWorkspace*) file->Get(workspaceName);

  // get the data  out of the file
  RooAbsData* data = w->data(dataName);
  
  // get the modelConfig (S+B) out of the file
  // and create the B model from the S+B model

  double n_0_events = 10;
  
  ModelConfig*  sbModel = (RooStats::ModelConfig*) w->obj(modelConfigName);
  sbModel->SetName("S+B Model");      
  RooRealVar* poi = (RooRealVar*) sbModel->GetParametersOfInterest()->first(); //n_0
  poi->setVal(n_0_events);  // set POI snapshot in S+B model for expected significance
  sbModel->SetSnapshot(RooArgSet(*poi));
  ModelConfig * bModel = (ModelConfig*) sbModel->Clone();
  bModel->SetName("B Model");  
  poi->setVal(0); //Background only, so n_0 = 0
  bModel->SetSnapshot(RooArgSet(*poi));
  sbModel->Print("v");
  bModel->Print("v");

  FrequentistCalculator   fc(*data, *sbModel, *bModel);
  const int ntoymc_int = 100;
  fc.SetToys(ntoymc_int,ntoymc_int);    //null (B), alt (S+B) 

  // create the test statistics
  ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
  // use one-sided profile likelihood
  profll.SetOneSidedDiscovery(true);
  
//  Double_t eval = profll.Evaluate(*data, *poi); //If you want to evaluate the test statistic of only your input data, uncomment these 3 lines.
//  cout<<"eval : "<<eval<<endl;
//  return;

  // configure  ToyMCSampler and set the test statistics
  ToyMCSampler *toymcs = (ToyMCSampler*)fc.GetTestStatSampler();
  toymcs->SetTestStatistic(&profll);
//  toymcs->SetNEventsPerToy(2000); //I haven't figured out what this function does.  My guess was that would be the number of events generated from the PDF, but results are telling me that might not be right.
  
  // run the test
  HypoTestResult * fqResult = fc.GetHypoTest();
  fqResult->Print();

  // plot test statistic distributions
  TCanvas *c1 = new TCanvas();
  HypoTestPlot * plot = new HypoTestPlot(*fqResult);
  plot->SetLogYaxis(true);
  plot->Draw();
  SamplingDistribution* null_dist = fqResult->GetNullDistribution();
  SamplingDistribution* alt_dist = fqResult->GetAltDistribution();
  
  
  SamplingDistPlot sampdist_plot;
  sampdist_plot.AddSamplingDistribution(null_dist);
  	

  	//Find median of distribution - integrate null distribution from median of alt distribution to infinity to determine p value, then use p->Z calc to find expected significance. 
  	Long64_t ntoymc = ntoymc_int;
  	const vector<Double_t> null_dist_vals = null_dist->GetSamplingDistribution();
  	Double_t null_dist_vals_d[ntoymc_int];
  	const vector<Double_t> alt_dist_vals = alt_dist->GetSamplingDistribution();
  	Double_t alt_dist_vals_d[ntoymc_int];
  	
  	for(int i=0;i<ntoymc_int;i++)
  	{
  		null_dist_vals_d[i] = null_dist_vals[i];
  		alt_dist_vals_d[i] = alt_dist_vals[i];
  	}
  	double q_median_null = TMath::Median(ntoymc, null_dist_vals_d);
  	double q_median_alt = TMath::Median(ntoymc, alt_dist_vals_d);
  	
  	cout<<"q_med_null : "<<q_median_null<<endl;
  	cout<<"q_med_alt : "<<q_median_alt<<endl;
  	
//	int bin_med_null = plot->GetTH1F(null_dist)->GetXaxis()->FindBin(q_median_null);
	int bin_med_alt = plot->GetTH1F(null_dist)->GetXaxis()->FindBin(q_median_alt);
	double bin_max_null = plot->GetTH1F(null_dist)->GetMaximumBin();
//	double q_med_null_P = plot->GetTH1F(null_dist)->Integral(bin_med_null, bin_max_null,"width");
	double q_med_alt_P = plot->GetTH1F(null_dist)->Integral(bin_med_alt, bin_max_null, "width");
//	double med_null_Z_roostat = PValueToSignificance(q_med_null_P);
	double med_alt_Z_roostat = PValueToSignificance(q_med_alt_P);
//	double null_Z_roostat = PValueToSignificance(q_null_P);
	cout<<"Integral from q_med_alt->inf : "<< q_med_alt_P<<endl;
	cout<<"Z significance from P value : "<<med_alt_Z_roostat<<endl;
  
  
  
  
  
  


}