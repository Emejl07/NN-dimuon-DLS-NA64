import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def RatioPlot(npy_file, histogramNameSim, filename, title, xmin, xmax):
    fig, axes = plt.subplots(2, 1, figsize=(6, 6), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    
    # Load the .npy file data for experimental data
    dimuon_features = np.load(npy_file)
    
    # Calculate the histogram of the npy file data
    histExp, binsExp = np.histogram(dimuon_features[:,3], bins=65, range=(xmin, xmax), density=False)
    
    # Open the ROOT file and access the histogram for simulation data
    file_sim = uproot.open("/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Output/job008747_out.root") # Use experimental data
    #file_sim = uproot.open("./../../DIMUONS/data/25042024_MM_real_seperated_beamtime_out.root") # Use simulated data
    histSim = file_sim[histogramNameSim].to_numpy()
    
    # Rebin the simulation histogram to match the experimental histogram
    histSim_rebinned, _ = np.histogram(histSim[1][:-1], bins=binsExp, weights=histSim[0], density=False)
    
    bin_centers = (binsExp[:-1] + binsExp[1:]) / 2
    
    # Calculate the errors
    errorExp = np.sqrt(histExp) / histExp
    errorSim = np.sqrt(histSim_rebinned) / histSim_rebinned
    
    # Plot the experimental (npy) and simulation histograms
    axes[0].step(bin_centers, histExp, 'r', where="mid", linewidth=1, label='NN selection')
    axes[0].errorbar(bin_centers, histExp, errorExp, lw=1, fmt="", linestyle='', color="k", capsize=0)
    
    axes[0].step(bin_centers, histSim_rebinned, 'b', where="mid", linewidth=1, label='Traditional selection')
    axes[0].errorbar(bin_centers, histSim_rebinned, errorSim, lw=1, fmt="", linestyle='', color="k", capsize=0)

    #axes[0].hist(bin_centers, weights=histExp, bins=len(bin_centers), color='r', alpha=1, label='NN selected data')

    # Plot the histogram for traditional selection
    #axes[0].hist(bin_centers, weights=histSim_rebinned, bins=len(bin_centers), color='b', alpha=1, label='Traditional selection')
    
    # Set x-axis range
    #axes[0].set_yscale('log')
    
    # Add labels and title
    axes[0].set_ylabel('Frequency', fontsize=20)
    axes[0].set_title(title, fontsize=21)
    leg = axes[0].legend(fontsize=14)
    leg.get_frame().set_edgecolor('none')
    leg.get_frame().set_facecolor('none')
    axes[0].text(0.7, 0.6, r"$\mathbf{e^+ e^- \rightarrow \mu^+ \mu^-}$", horizontalalignment='center', verticalalignment='center', transform=axes[0].transAxes, fontsize=18, fontweight='bold')
    
    # Plot the ratio plot
    #ratio = (histSim_rebinned / normSim) / (histExp / normExp)
    ratio = histExp / histSim_rebinned
    #ratio_error = ratio * np.sqrt((errorSim / (histSim_rebinned / normSim))**2 + (errorExp / (histExp / normExp))**2) * 10
   
    ratio_error = errorSim * errorExp
    
    axes[1].axhline(y=1.2, color='k', linestyle='--')
    axes[1].axhline(y=1.0, color='k', linestyle='-')
    axes[1].axhline(y=0.8, color='k', linestyle='--')
   
    # Error bar plot on the second axis (axes[1])
    axes[1].errorbar(
        bin_centers, 
        ratio, 
        yerr=ratio_error, 
        xerr=(xmax - xmin) / (2 * len(bin_centers)), 
        lw=0,  # Line width for the marker lines (not used here)
        fmt='s',  # Marker style: square
        color='black',  # Marker color
        ecolor='black',  # Error bar color
        elinewidth=1,  # Error bar line width
        markersize=3,
        capsize=0  # No caps on the error bars
    ) 
    #axes[1].errorbar(bin_centers, ratio, ratio_error, xerr=(xmax - xmin) / (2 * len(bin_centers)), lw=2, fmt='', linestyle='', color='g', capsize=0)
    
    axes[1].set_xlabel('Energy [GeV]', fontsize=20)
    axes[1].set_ylabel('Ratio', fontsize=20)
    axes[1].set_ylim(0.1, 10)
    axes[1].set_xlim(xmin, xmax)
    axes[1].set_ylim(0, 2)

    axes[1].xaxis.set_major_locator(MaxNLocator(integer=True))

    # Set the y-axis to display only integer values
    axes[0].yaxis.set_major_locator(MaxNLocator(integer=True))
    axes[1].yaxis.set_major_locator(MaxNLocator(integer=True))
    
    # Show plot
    axes[0].tick_params(axis='y', labelsize=15)
    plt.tight_layout(pad=2.0)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0, wspace=0)
    plt.savefig(filename)
    plt.show()

# Example usage:
RatioPlot("/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Data/Input/NN_Selected_experiment.npy", "Dimuons/HCAL2_11_T_18", "/eos/user/e/ezaya/simulation_output/NN/NN_dimuons/Plots/HCAL2_RatioPlot_Exp.pdf", r"Central cell of HCAL$_2$", 0, 20)
