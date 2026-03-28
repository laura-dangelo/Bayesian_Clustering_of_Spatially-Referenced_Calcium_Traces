# Bayesian_Clustering_of_Spatially-Referenced_Calcium_Traces

This repository contains the code to replicate the results in the paper "_Decoding Neuronal Ensembles from Spatially-Referenced Calcium Traces: A Bayesian Semiparametric Approach_" by D'Angelo, Denti, Canale, and Guindani (available on [arXiv](https://arxiv.org/abs/2508.09576)). It also contains the R package `bSCDCsampler`, which implements the Gibbs sampler for posterior inference (see [Installation](#installation)).

## Structure of the Repository

Main structure:
<ul>
  <li> Code 
    <ul>
      <li> 01_data_preprocessing : preparation of the data from the original .mat file. </li>
      <li> 02_data_analysis : contains the scripts to replicate the real data analyses.</li>
      <li> 03_simulation_study : contains the scripts to replicate the simulation studies reported in the Supplementary Material.</li>
      <li> bSCDCsampler : R package implementing the Gibbs sampler.</li>
    </ul>
  </li>
  <li> Data </li>
</ul>

The `Code` folder contains the R project and should be used as the main directory for running the code.
The structure of the repository is designed to facilitate the streamlined execution of all scripts required to replicate the analyses. Please do not change the folder names or paths to avoid errors. 

The `Data` folder contains the data used in the analyses. 

Large files (typically, outputs of the Gibbs sampler algorithm or the deconvolution algorithm) are excluded from the repository and are available in the Google Drive [folder](https://drive.google.com/drive/folders/1-1xf57mZBc1usA-iCZGkp4KPF8oX_5mV?usp=sharing). The name of each subfolder corresponds to the path of the files in the repository. Note that these files are not necessary to replicate the analyses.  




## Using the Repo
All analyses can be replicated starting from the original `.mat` file (`M3424F_data_togo_neuron_behav_multiTrials_072621.mat`), which is freely available for download from the [Mendeley data repository](https://doi.org/10.17632/tnbwfw2pg2.2). However, some computations can be time-consuming, and we do not recommend starting from scratch.<br/>


**We provide RDS files containing precomputed outputs of the inference procedure to facilitate replication of the analyses.** <br/>
**To replicate the analyses and plots of the paper starting from these precomputed quantities, you only need to execute the scripts marked by "RUN".**



## Structure of the `Code` folder
### Code > 02_data_analysis
The key elements of this folder are the following:

├── **01_bSCDC_individual_trials** &ensp;&ensp;&ensp;&ensp;# single-window analyses <br/>
│   ├── 01_plot_data.R <br/>
│   ├── 02_gibbs_sampler.R <br/>
│   ├── 03_RUN_analyze_individual_runs.R <br/>
│   └── 04_RUN_PSM_neurons.R <br/>
|<br/>
├── **02_bSCDC_neuronal_response_to_position** &ensp;&ensp;&ensp;&ensp;# analyses on the overall neuronal activity <br/>
│   ├── 00_auxiliary_functions_DONT_RUN.R <br/>
│   ├── 01_RUN_plot_cluster_complexity.R <br/>
│   └── 02_RUN_plot_activation_maps.R <br/>
|<br/>
├── **03_bSCDC_sensitivity_study** &ensp;&ensp;&ensp;&ensp;# sensitivity studies on the prior parameters <br/>
│   ├── A_sensitivity_alow_01_gibbs_sampler.R <br/>
│   ├── A_sensitivity_alow_02_extract_results.R <br/>
│   ├── A_sensitivity_alow_03_RUN_compare_results.R <br/>
│   ├── B_sensitivity_PSBP_01_gibbs_sampler.R <br/>
│   ├── B_sensitivity_PSBP_02_extract_results.R <br/>
│   ├── B_sensitivity_PSBP_03_RUN_compare_results.R <br/>
│   └── B_sensitivity_PSBP_04_toy_example.R <br/>
|<br/>
├── **04_comparison_two-step** &ensp;&ensp;&ensp;&ensp;# comparison with the two-step deconvolution and clustering <br/>
│   ├── 01_spike_detection_JW.R <br/>
│   ├── 02_kmeans.R <br/>
│   └── 03_RUN_plot_results.R <br/>
|<br/>
└── **05_comparison_GPFA** &ensp;&ensp;&ensp;&ensp;# comparison with Bayesian Gaussian process factor analysis <br/>




### Code > 03_simulation_study
The key elements of this folder are the following:

├── **01_sensitivity_study** &ensp;&ensp;&ensp;&ensp;# sensitivity studies on synthetic data <br/>
│   ├── 00_auxiliary_functions_DONT_RUN.R <br/>
│   ├── 01_simulate_data.R <br/>
│   ├── 02_gibbs_sampler.R <br/>
│   ├── 03_extract_results.R <br/>
│   └── 04_RUN_plot_results.R <br/>
|<br/>
├── **02_comparison_two-step_synthetic_data** &ensp;&ensp;&ensp;&ensp;# comparison with the two-step method on synthetic data <br/>
│   ├── 01_spike_detection_JW.R <br/>
│   ├── 02_kmeans_compute_ARI.R <br/>
│   ├── 03_extract_results.R <br/>
│   └── 04_RUN_plot_results.R <br/>
|<br/>
└── **03_computational_cost** &ensp;&ensp;&ensp;&ensp;# simulation study to evaluate the computational cost <br/>
&ensp;&ensp;&ensp;├── 01_simulate_data.R <br/>
&ensp;&ensp;&ensp;├── 02_gibbs_sampler.R <br/>
&ensp;&ensp;&ensp;└── 03_RUN_extract_results.R <br/>



## Installation
The Code folder contains the binary file `bSCDCsampler_0.0.1.tar.gz`. You can install the package on R using
``` install.packages("bSCDCsampler_0.0.1.tar.gz", repos = NULL, type="source") ```

