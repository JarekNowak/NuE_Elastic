import ROOT
import array

# 1. Setup the binning based on the plot (0 to 2.5 GeV, 50 MeV bins)
xmin = 0.0
xmax = 2.5        # digitized values cover 0–2.5 GeV only
bin_width = 0.05  # 50 MeV
nbins = int((xmax - xmin) / bin_width)

# 2. Approximate digitized values for the muon neutrino flux (blue line)
# These represent the 'Flux' values at the center of each bin
# Estimated from the log scale in Fig 1.
v_mu_flux_values = [
    2.5, 3.8, 5.2, 6.5, 7.5, 8.0, 8.2, 8.1, 7.8, 7.2,  # 0.0 - 0.5 GeV
    6.5, 5.8, 5.2, 4.5, 3.9, 3.4, 3.0, 2.6, 2.2, 1.8, # 0.5 - 1.0 GeV
    1.5, 1.3, 1.1, 0.9, 0.75, 0.62, 0.52, 0.44, 0.38, 0.32, # 1.0 - 1.5 GeV
    0.27, 0.23, 0.20, 0.17, 0.15, 0.13, 0.11, 0.10, 0.09, 0.08, # 1.5 - 2.0 GeV
    0.07, 0.065, 0.06, 0.058, 0.055, 0.052, 0.05, 0.048, 0.046, 0.045 # 2.0 - 2.5 GeV
]

# 3. Create the ROOT file and Histogram
file = ROOT.TFile("sbnd_flux_data.root", "RECREATE")
h_flux = ROOT.TH1D("h_flux", "Muon Neutrino Flux at SBND;E (GeV);Flux [(10^{6} POT m^{2} 50 MeV)^{-1}]", nbins, xmin, xmax)

# Fill the histogram
for i, val in enumerate(v_mu_flux_values):
    h_flux.SetBinContent(i + 1, val)

# 4. Save and Close
h_flux.Write()
file.Close()

print("File 'sbnd_flux_data.root' has been created with histogram 'h_flux'.")
