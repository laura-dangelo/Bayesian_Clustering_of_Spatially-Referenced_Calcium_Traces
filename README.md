# Bayesian_Clustering_of_Spatially-Referenced_Calcium_Traces

This repository contains the code for replication of the results in the paper "Decoding Neuronal Ensembles from Spatially-Referenced Calcium Traces: A Bayesian Semiparametric Approach"  by D'Angelo, Denti, Canale, and Guindani (available on [arXiv](https://arxiv.org/abs/2508.09576)). It also contains the R package `bSCDCsampler`, which implements the Gibbs sampler used for posterior inference (see [Installation](#installation)).

## Structure of the repository
The `Code` folder contains the R project and should be used as the main directory to run the code.
The structure of the repository is designed to facilitate the streamlined execution of all scripts required to replicate the analyses. Please do not change the folders’ names or paths to avoid errors. Each folder should be run following the order given by the numbers in the name; similarly, within each folder, scripts should be executed following the numbering.

Key structure:
<ul>
  <li> Code 
    <ul>
      <li> 01_data_preprocessing : preparation of the data from the original .mat file. </li>
      <li> 02_data_analysis : contains all the scripts to replicate the analyses on the real data.</li>
      <li> 03_simulation_study : contains the scripts to replicate the simulation studies reported in the Supplementary Material.</li>
      <li> bSCDCsampler : R package implementing the Gibbs sampler.</li>
    </ul>
  </li>
  <li> Data </li>
  <li> Manuscript </li>
</ul>


## Data
All analyses can be replicated starting from the original .mat file (M3424F_data_togo_neuron_behav_multiTrials_072621.mat), which is freely downloadable from the [Mendeley data repository](https://doi.org/10.17632/tnbwfw2pg2.2). However, some computations can be time-consuming, and we do not recommend starting from scratch.
To this end, we provide intermediate .RDS outputs to ease replication. Large files are available in the Google Drive [folder](https://drive.google.com/drive/folders/1-1xf57mZBc1usA-iCZGkp4KPF8oX_5mV?usp=sharing). The name of each folder corresponds to the path where the files should be copied.

- `Data` contains the data extracted from the original Matlab file. They can be extracted from the .mat file using the `Code > 01_data_preprocessing` scripts. 
- `*/results` contain the output of the Gibbs sampler algorithm.



## Code
#### Code > 02_data_analysis
The key elements of this folder are the following:

├── **01_bSCDC_individual_trials** &ensp;# single-window analyses <br/>
│   ├── 01_plot_data.R <br/>
│   ├── 02_run_gibbs.R <br/>
│   ├── 03_analyze_individual_runs.R <br/>
│   └── 04_PSM_neurons.R <br/>
├── **02_bSCDC_neuronal_response_to_position** &ensp;# analyses on the overall neuronal activity <br/>
│   ├── 00_auxiliary_functions_DONT_RUN.R <br/>
│   ├── 01_plot_cluster_complexity.R <br/>
│   └── 02_plot_activation_maps.R <br/>
├── **03_bSCDC_sensitivity_study** &ensp;# sensitivity studies on the prior parameters <br/>
│   ├── A_sensitivity_alow_01_run_gibbs.R <br/>
│   ├── A_sensitivity_alow_02_analyze_results.R <br/>
│   ├── A_sensitivity_alow_03_compare_results.R <br/>
│   ├── B_sensitivity_PSBP_01_run_gibbs.R <br/>
│   ├── B_sensitivity_PSBP_02_analyze_results.R <br/>
│   ├── B_sensitivity_PSBP_03_compare_results.R <br/>
│   └── B_sensitivity_PSBP_04_toy_example.R <br/>
├── **04_comparison_two-step** &ensp;# comparison with the two-step deconvolution and clustering <br/>
│   ├── 01_spike_detection_JW.R <br/>
│   ├── 02_run_kmeans.R <br/>
│   └── 03_plot_results.R <br/>
└── **05_comparison_GPFA** &ensp;# comparison with Bayesian Gaussian process factor analysis <br/>




#### Code > 03_simulation_study
The key elements of this folder are the following:

├── 01_sensitivity_study&ensp;# sensitivity studies on synthetic data <br/>
│   ├── 00_auxiliary_functions_DONT_RUN.R <br/>
│   ├── 01_simulate_data.R <br/>
│   ├── 02_run_gibbs.R <br/>
│   ├── 03_extract_results.R <br/>
│   └── 04_plot_results.R <br/>
├── 02_comparison_two-step_synthetic_data&ensp;# comparison with the two-step method on synthetic data <br/>
│   ├── 01_spike_detection_JW.R <br/>
│   ├── 02_run_kmeans_compute_ARI.R <br/>
│   ├── 03_extract_results.R <br/>
│   └── 04_plot_results.R <br/>
└── 03_computational_cost&ensp;# simulation study to evaluate the computational cost <br/>
&ensp;&ensp;&ensp;├── 01_simulate_data.R <br/>
&ensp;&ensp;&ensp;├── 02_run_gibbs.R <br/>
&ensp;&ensp;&ensp;└── 03_extract_results.R <br/>



## Installation
The Code folder contains the binary file `bSCDCsampler_0.0.1.tar.gz`. You can install the package on R using
``` install.packages(bSCDCsampler_0.0.1.tar.gz, repos = NULL, type="source") ```

