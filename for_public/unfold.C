void unfold( ){

  TProfile::Approximate();
  gStyle->SetOptStat(0);
  TH1D *meas, *ntrk_gen, *nch_gen;
  TH2D *resp;
  TLegend* leg = 0x0;
  
  TCanvas* c = new TCanvas("unfolding", "N_{ch} unfolding", 1000, 800);
  TCanvas* c2 = new TCanvas("distributions", "N_{ch}, Ntrk distributions", 1000, 800);
  c2->SetLogy();
  
  
  
  
  // load measured ntrk distribution
  
  
  TFile fm("data.root");
  TList* lm = (TList*) fm.Get("jpsi2eeHistos");
  TList* llm = (TList*) lm->FindObject("Event_AfterCuts");
  meas = (TH1D*) llm->FindObject("Trk_corr");
  meas->SetDirectory(0);
  meas->SetTitle( "N_{trk, corr.} distribution" );
  
  
  meas->Scale(1./meas->GetSumOfWeights() );
  fm.Close();
  
  
  
  // load response matrix from Monte Carlo.
  // expects a TH2D with number of tracklets on x-axis
  // and Nch on y axis
  // x axis should have same binning as ntrk distribution in data
  
  TFile fr("mc.root");
  TList* lr = (TList*) fr.Get("jpsi2eeHistos");
  TList* llr = (TList*) lr->FindObject("Event_AfterCuts");
  resp = (TH2D*) llr->FindObject( "NchEta1_trk");
  resp->SetDirectory(0);
  fr.Close();
  
  
  
  // get the ntrk and nch distributions in MC from projections of the response matrix
  
  ntrk_gen = (TH1D*) resp->ProjectionX();
  nch_gen = (TH1D*) resp->ProjectionY();
  
  ntrk_gen->Scale( 1./ntrk_gen->GetSumOfWeights() );
  nch_gen->Scale( 1./nch_gen->GetSumOfWeights() );
  
  
  
  
  // ---------------------------------------------
  // here somes the unfolding part
  
  RooUnfoldResponse* rooResp = new RooUnfoldResponse ( 0, 0 , resp);
  RooUnfoldBayes   unfold (rooResp, meas, 4);
  
  TH1D* nch_unfolded = (TH1D*) unfold.Hreco();
  nch_unfolded->Scale(1./nch_unfolded->GetSumOfWeights() );

  
  
  // unfolding done
  // -----------------------------------------------------
  
  
  
  
  nch_gen->SetLineColor(kRed);
  nch_gen->SetMarkerColor(kRed);
  nch_gen->SetMarkerStyle(21);
  
  nch_unfolded->SetLineColor(kGreen+2);
  nch_unfolded->SetMarkerColor(kGreen+2);
  nch_unfolded->SetMarkerStyle(22);
    
    
  c2->cd();
    
  leg = new TLegend(.6,.7,.88,.88);
  leg->SetBorderSize(0);  
  leg->AddEntry(nch_gen, "original N_{ch} distr.", "PL");
  leg->AddEntry(nch_unfolded, "unfolded N_{ch} distr.", "PL");
  nch_gen->GetYaxis()->SetRangeUser(  1e-6   ,.1)   ;
  nch_gen->GetXaxis()->SetRangeUser(   0., 120.   )   ;
  nch_gen->GetYaxis()->SetTitle( "Probability")   ;
  
  nch_gen->Draw();
  nch_unfolded->Draw("same");
  leg->Draw();
 
 c2->SaveAs( "nch_distributions.pdf");
 
 
 
 nch_unfolded->Sumw2();
 nch_gen->Sumw2();
 
 
 TH1D* ratio_nch = (TH1D*)nch_unfolded->Clone("ratio_nch");
 ratio_nch->Divide(nch_gen);
 
 ratio_nch->SetTitle(  "Ratio of unfolded and generated N_{ch} distribution"   );
 
 
 c->cd();
 ratio_nch->GetYaxis()->SetRangeUser(0.,10.);
 ratio_nch->GetYaxis()->SetTitle(  "Ratio unfolded / generated N_{ch}"   );
 ratio_nch->Draw();
 
 
 c->SaveAs( "ratio_nch.pdf" );

// scale response matrix with unfolded nch distribution
  c->SetLogy(0);

  TH2D* resp_corr = (TH2D*)resp->Clone("corrected");
  for(int ix = 0; ix< resp->GetNbinsX(); ++ix){
    for(int iy =0; iy< resp->GetNbinsY(); ++iy){
      if(nch_gen->GetBinContent(iy)  ){
        resp_corr->SetBinContent( ix, iy, resp_corr->GetBinContent(ix, iy) * nch_unfolded->GetBinContent(iy) / nch_gen->GetBinContent(iy) );
      }
    }
  }
  resp_corr->Draw("colz");
  
  resp_corr->SaveAs( "detector_response_mid_corr.root" );
  
  
  // profiles of  mean n_ch as a function of n_trk
  // prof_orig : using nch distribution as in pythia
  // prof_corr: using unfolded nch distribution
  // prof_corr_fit: using unfolded nch,and pol2 fit of mean nch as a function of ntrk at high ntrk 
  // (nch is fitted between 40< ntrk<80, and result of fit used in the range 80<ntrk<120)
  
  
  // prof_corr_fit is the one I used in my analysis
  
  
  
  TProfile* prof_orig = resp->ProfileX("prof_orig" );
  TProfile* prof_corr = resp_corr->ProfileX("prof_corr");
  TH1D* prof_corr_fit = new TH1D("prof_corr_fit", "prof_corr_fit",
                              prof_corr->GetNbinsX() ,  
                              prof_corr->GetXaxis()->GetBinLowEdge(1) ,
                              prof_corr->GetXaxis()->GetBinUpEdge( prof_corr->GetNbinsX() ) )  ; 
  
  TF1 * pol2 = new TF1("pol2","pol2",40.,80.);
  
  prof_corr->Fit(pol2,"","",40.,80.);
  
  
  TF1 * pol22 = new TF1("pol22","pol2",80.,120.);
  
  pol22->SetParameter(0, pol2->GetParameter(0) );
  pol22->SetParameter(1, pol2->GetParameter(1) );
  pol22->SetParameter(2, pol2->GetParameter(2) );
  
  
  
  for( int j= prof_corr_fit->FindBin(0.); j <= prof_corr_fit->FindBin(80.); ++j  ){
    prof_corr_fit->SetBinContent( j,  prof_corr->GetBinContent(j) );
  }
  for( int j= prof_corr_fit->FindBin(80.); j <= prof_corr_fit->FindBin(120.); ++j  ){
    double y = pol22->Eval( prof_corr_fit->GetBinCenter(j) );
    prof_corr_fit->SetBinContent( j,  y );
  }
  
  
  meas->SetLineColor(kBlue);
  meas->SetMarkerColor(kBlue);
  meas->SetMarkerStyle(20);
  meas->GetXaxis()->SetRangeUser(0., 120.);
  meas->GetXaxis()->SetTitle("N_{trk, corr.}");
  meas->GetYaxis()->SetTitle("Probability");
  
  
  ntrk_gen->SetLineColor(kRed);
  ntrk_gen->SetMarkerColor(kRed);
  ntrk_gen->SetMarkerStyle(21);
  
  
  
  
  TH1D* ntrk_unfolded = resp_corr->ProjectionX();
  ntrk_unfolded->Scale(1./ntrk_unfolded->GetSumOfWeights() );
  ntrk_unfolded->SetLineColor(kGreen+2);
  ntrk_unfolded->SetMarkerColor(kGreen+2);
  ntrk_unfolded->SetMarkerStyle(22);
  
  c2->cd();
  
  meas->Draw();
  ntrk_gen->Draw("same");
  ntrk_unfolded->Draw("same");
  leg = new TLegend(.6,.7,.88,.88);
  leg->SetBorderSize(0);
  leg->AddEntry(meas, "Data", "PL");
  leg->AddEntry(ntrk_gen, "MC (original N_{ch})", "PL");
  leg->AddEntry(ntrk_unfolded, "MC (unfolded N_{ch})", "PL");
  
  leg->Draw();
  c2->SaveAs( "ntrk_distributions.pdf");
  
    
    
        
  int nx = resp->GetNbinsX();
  double minx = resp->GetXaxis()->GetBinLowEdge( 1 );
  double maxx = resp->GetXaxis()->GetBinUpEdge( resp->GetNbinsX() );
    
  // histograms of alpha factors, i.e. mean n_ch/n_trk as a function of n_trk
  // alpha_orig : using nch distribution as in pythia
  // alpha_corr: using unfolded nch distribution
  // alpha_corr_fit: using unfolded nch,and pol2 fit of mean nch as a function of ntrk at high ntrk 
  
    
  TH1D* alpha_orig = new TH1D("alpha_orig", "Ratio between N_{ch} and N_{trk, corr.} from MC simulations",nx, minx, maxx );
  TH1D* alpha_corr = new TH1D("alpha_corr", "alpha_corr",nx, minx, maxx );
  TH1D* alpha_corr_fit = new TH1D("alpha_corr", "alpha_corr",nx, minx, maxx );
  alpha_orig->Sumw2();
  alpha_corr->Sumw2();
  alpha_corr_fit->Sumw2();
    
    
    
  for( int i=0; i< prof_orig->GetNbinsX(); ++i){
    
    double x = prof_orig->GetBinCenter(i+1);
    double y_orig = prof_orig->GetBinContent(i+1);
    double y_corr = prof_corr->GetBinContent(i+1);
    double y_corr2 = prof_corr_fit->GetBinContent(i+1);
    
    double ye_orig = prof_orig->GetBinError(i+1);
    double ye_corr = prof_corr->GetBinError(i+1);
    double ye_corr2 = prof_corr_fit->GetBinError(i+1);
    if(x){
      alpha_orig->SetBinContent( i+1 , y_orig/x );
      alpha_corr->SetBinContent( i+1 , y_corr/x );
      alpha_corr_fit->SetBinContent( i+1 , y_corr2/x );
      alpha_orig->SetBinError( i+1 , ye_orig/x );
      alpha_corr->SetBinError( i+1 , ye_corr/x );
      alpha_corr_fit->SetBinError( i+1 , 0.001);
    }
  }
    
    
    
  alpha_corr->SetLineColor(kBlue);
  alpha_corr->SetMarkerColor(kBlue);
  alpha_corr->SetMarkerStyle(20);
  
  alpha_corr_fit->SetLineColor(kBlack);
  alpha_corr_fit->SetMarkerColor(kBlack);
  alpha_corr_fit->SetMarkerStyle(24);
  
  alpha_orig->SetLineColor(kRed);
  alpha_orig->SetMarkerColor(kRed);
  alpha_orig->SetMarkerStyle(21);
  
  alpha_orig->GetXaxis()->SetRangeUser( 1.,120. );
  alpha_orig->GetYaxis()->SetRangeUser( 1.,1.5 );
  alpha_orig->GetXaxis()->SetTitle("N_{trk, corr}");
  alpha_orig->GetYaxis()->SetTitle("#LTN_{ch}#GT / N_{trk, corr}");
    
  
  c->cd();
  alpha_orig->Draw();
  c->SaveAs(  "alpha_orig.pdf"  );
  
  alpha_corr->Draw("same");  
  alpha_corr_fit->Draw("same");  
  leg = new TLegend(.4,.7,.88,.88);
  leg->SetBorderSize(0);
  leg->AddEntry(alpha_orig, "original", "PL");
  leg->AddEntry(alpha_corr, "unfolded", "PL");
  leg->AddEntry(alpha_corr_fit, "unfolded, fitted at high N_{trk, corr.}", "PL");
  leg->Draw();
  c->SaveAs( "alpha.pdf" );
  

    
  alpha_corr->SaveAs( "alpha_mid_nofit.root" );
  alpha_corr_fit->SaveAs( "alpha_mid.root" );
  prof_corr_fit->SaveAs( "prof_mid.root" );

}
