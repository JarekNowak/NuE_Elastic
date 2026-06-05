#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>

void get_sbnd_flux() {
    // 1. Setup Binning: 0 to 3.0 GeV with 50 MeV (0.05 GeV) steps
    const double xmin = 0.0;
    const double xmax = 3.0;
    const int nbins = 60; // (3.0 / 0.05)

    // 2. Digitized y-axis values for the Muon Neutrino (Blue Curve)
    // Units of raw_vals: [10^6 POT * m^2 * 50 MeV]^-1
    std::vector<double> raw_vals = {
        0.5, 2.5, 4.8, 6.2, 7.5, 8.2, 8.5, 8.4, 8.2, 7.8, // 0.0 - 0.5 GeV
        7.2, 6.5, 5.8, 5.0, 4.4, 3.8, 3.3, 2.8, 2.4, 2.1, // 0.5 - 1.0 GeV
        1.8, 1.5, 1.3, 1.1, 0.95, 0.82, 0.70, 0.60, 0.52, 0.45, // 1.0 - 1.5 GeV
        0.38, 0.32, 0.28, 0.24, 0.21, 0.18, 0.16, 0.14, 0.12, 0.11, // 1.5 - 2.0 GeV
        0.10, 0.09, 0.08, 0.07, 0.065, 0.06, 0.055, 0.05, 0.045, 0.04, // 2.0 - 2.5 GeV
        0.038, 0.035, 0.032, 0.03, 0.028, 0.026, 0.024, 0.022, 0.02, 0.018 // 2.5 - 3.0 GeV
    };

    // 3. Apply the conversion factor: 2e-9
    double factor = 2.0e-9;

    // 4. Create ROOT file and Histogram
    TFile *f = new TFile("sbnd_flux_converted.root", "RECREATE");
    TH1D *h_flux = new TH1D("h_flux", "SBND #nu_{#mu} Flux;E_{#nu} [GeV];#Phi [#nu/cm^{2}/POT/GeV]", nbins, xmin, xmax);

    for (int i = 0; i < nbins; ++i) {
        h_flux->SetBinContent(i + 1, raw_vals[i] * factor);
        // Setting a small dummy error for visualization (e.g., 2% relative)
        h_flux->SetBinError(i + 1, (raw_vals[i] * factor) * 0.02);
    }

    // 5. Save and Display
    h_flux->SetLineColor(kBlue+1);
    h_flux->SetLineWidth(2);
    h_flux->Write();
    
    std::cout << "ROOT file 'sbnd_flux_converted.root' created." << std::endl;
    std::cout << "Peak Flux: " << h_flux->GetMaximum() << " [nu/cm^2/POT/GeV]" << std::endl;

    // Draw result
    TCanvas *c1 = new TCanvas("c1", "SBND Flux", 800, 600);
    h_flux->Draw("HIST E");
}
