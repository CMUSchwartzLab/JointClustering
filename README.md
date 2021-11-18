# JointClustering
This repository aims to generate joint clusters for single cell sequencing copy number data and miFISH data (multiplex interphase fluorescence in situ hybridization). The code is based on python3.7. The method is introduced in Xuecong Fu, Haoyun Lei, Yifeng Tao, Kerstin Heselmeyer-haddad, Irianna Torres, Michael Dean, Thomas Ried, and Russell Schwartz.Journal of Computational Biology.Nov 2021.1035-1051.http://doi.org/10.1089/cmb.2021.0255

## Pre-clustering
The method breaks into two steps. First separate clustering needs to be conducted for both types of data. The input includes the csv files for both datasets (-scs and -fish) and the cluster numbers for both datasets(-k_scs and -k_fish). The following code is an example for pre-clustering. You can also change the method for clustering. The default method for both datasets are Kmeans but Gaussian mixture model is also an option.
```
python clustering.py 
        -scs ${csv of scSeq data} -fish ${csv of scFISH data} 
        -k_scs ${number of clusters for scSeq} -k_fish ${number of clusters for scFISH} 
        -method_scs ${clustering method for scSeq} -method_fish ${clustering method for scFISH} 
```
The output would be two csv files including the cluster centers for both datasets ("means_fish.csv" and "means_scs.csv") and two csv files including the estimated frequencies of the clusters for both datasets ("freqs_fish.csv" and "freqs_scs.csv"). You can also generate these four files with your own clustering method.

## MCMC Sampling
After the pre-clustering, we apply the MCMC sampling for searching the joint matching with the maximum likelihood. The inputs include the output files from pre-clustering, the output folder, iterations and probe_idx.pickle which is a list showing the genome indices of scSeq data corresponding to the FISH probes. You can also change the scaling factor for the gaussian error distribution which helps estimate the likelihood and the bound iterations and likelihood which are used in early stop determination. You can run the following example code if you decide to set everything as default.

```
python mcmc_newll.py
    -sp / --scseq_profiles ${cluster centers of scSeq}
    -fp / --scfish_profiles ${cluster centers of scFISH}
    -sf / --scseq_freqs ${cluster frequencies of scSeq}
    -ff / --scfish_freqs ${cluster frequencies of scFISH}
    -o / --output_folder ${output folder}
    -sigma / --sigma ${The scaling factor for gaussian error distribution}
    -trial / --trial ${The number of trial for training, used for naming, default is 1}
    -i / --iterations ${max iterations}
    -bound_ll / --bound_ll ${Bound log likelihood for early stop, default is -10^7}
    -bound_iter / --bound_iter ${Bound iteration for early stop, default is 0}
    -dataset / --dataset ${the dataset name, used for naming}
    -probe_idx / --probe_index_list ${the pickle file contains a list which shows the genome indices of scSeq data corresponding to the FISH probes, default is probe_idx.pickle}
```
You can also try ``` python mcmc_newll.py ``` if you set all options to default.

The output would be the whole genome unnormalized cluster centers profiles (means_scs_ploidy.csv) and the corresponding cluster frequencies (weights_scs_ploidy.csv), along with the trained model saved in a pickle file.

## Simulation
We also provide codes for simulation of both types of data. You can customize several parameters such as the level of perturbation of single cells from the ground truth cluster centers, the ground truth cluster number, the mutation rate etc.
```
python simulation_data.py
    -perturb_scs / --perturb_scs ${Perturbation rate of scSeq from the cluster centers when method is int, default=0.05}
    -perturb_fish / --perturb_fish ${Perturbation rate of scFISH from the cluster centers when method is int, default=0.05}
    -scs / --scseq_data ${template scSeq data if any}
    -chromo_info / --chromo_info ${A csv file with chromosome info of each feature if any}
    -probe_idx / --probe_index_list ${a pickle file which saves a list which shows the genome indices of scSeq data corresponding to the FISH probes, default is probe_idx.pickle}
    -k / --num_leaves ${number of leaf nodes}
    -lam / --lambda ${The mutation rate}
    -sigma2 / --sigma2 ${The noise between the ground truth scFISH marker copy numbers and corresponding scSeq features, default is 0.2}
    -n / --num_features ${number of features}
    -wgd / --wgd ${Boolean indicating whether or not to incorporate whole genome duplication}
    -perturb_method / --perturb_method ${Perturbation method, either int or gaussian, default is int}
    -perturb_sigma2 / --perturb_sigma2 ${Perturbation sigma2 when method is gaussian, default is 0.3}
    -m_scs / --num_scs ${total number of scSeq cells, default is 200}
    -m_fish / --num_fish ${total number of scFISH cells, default is 450}
```
You can also try ``` python simulation_data.py ``` if you set all options to default. A csv file with scSeq data (scs_simu.csv) and a csv file with miFISH data (fish_simu.csv) will be generated.

