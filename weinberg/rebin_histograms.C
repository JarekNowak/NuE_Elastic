#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

void rebin_histograms(const char* infile = "weinberg/inputs.root",
                      const char* outfile = "weinberg/outputs_reduced.root") {
    // 1. Open the source file
    TFile *fIn = TFile::Open(infile, "READ");
    if (!fIn || fIn->IsZombie()) {
        printf("Error: Could not open %s\n", infile);
        return;
    }



   int rebin =1;

    // 2. Retrieve the histograms
    TH1D *h_data = (TH1D*)fIn->Get("h_data");
    TH1D *h_eff  = (TH1D*)fIn->Get("h_eff");
    TH1D *h_bkg  = (TH1D*)fIn->Get("h_bkg");

    
    
    TH2D *h_smear = (TH2D*)fIn->Get("h_smear");

    if (!h_data || !h_smear) {
        printf("Error: One or more histograms not found!\n");
        return;
    }

    // 3. Create the output file
    TFile *fOut = TFile::Open(outfile, "RECREATE");

    // 4. Rebin h_data (1D)
    // Rebin(ngroup) merges 'ngroup' adjacent bins into one.
    // To reduce the number of bins by half, we use 2.
//    TH1D *h_data_new = (TH1D*)h_data->Rebin(2, "h_data_rebin");
    TH1D *h_data_new = (TH1D*)h_data->Rebin(rebin, "h_data");
    TH1D *h_eff_new  = (TH1D*)h_eff ->Rebin(rebin, "h_eff");
    TH1D *h_bkg_new  = (TH1D*)h_bkg ->Rebin(rebin, "h_bkg");


  h_eff_new->Scale(double(1.0/rebin));




    // 5. Rebin h_smear (2D)
    // Rebin2D(nxgroup, nygroup) merges bins in both dimensions.
    TH2D *h_smear_new = (TH2D*)h_smear->Rebin2D(rebin, rebin, "h_smear");


    h_data_new->Sumw2(kFALSE); 
    h_data_new->Sumw2(kTRUE);  // Recalculates based on current bin content

    h_eff_new->Sumw2(kFALSE);
  //  h_eff_new->Sumw2(kTRUE);  // Recalculates based on current bin content

    h_bkg_new->Sumw2(kFALSE);
    h_bkg_new->Sumw2(kTRUE);  // Recalculates based on current bin content




    h_smear_new->Sumw2(kFALSE);
    h_smear_new->Sumw2(kTRUE);


    // 6. Save and Close
    fOut->Write();
    fOut->Close();
    fIn->Close();

    printf("Success: Histograms rebinned and saved to %s\n", outfile);
}
