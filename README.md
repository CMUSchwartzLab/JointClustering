# JointClustering
This repository aims to generate joint clusters for single cell sequencing copy number data and single cell miFISH data (multiplex interphase fluorescence in situ hybridization). The code is based on python3.7. The method is introduced in [1].

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

## Phylogenetic Reconstruction

We used a heuristic method which combines the tree reconstruction method FISHTrees and a tree summary method ASTRAL to reconstruct the phylogenetic tree from the joint clusters.

### Requirements
* FISHTrees[2,3] (ftp://ftp.ncbi.nlm.nih.gov/pub/FISHtrees)
* ASTRAL[4] (https://github.com/smirarab/ASTRAL/)

```
python run_fishtree_mcmc.py 
    -d / --depth ${Sequencing depth for phylogenetic reconstruction}
    -p / --num_probe ${Number of probes per subset for each run of FISHTree, maximum=10 for FISHTree software}
    -cpu / --cpu_parellel ${number of cpu for running FISHTRee in parallel}
    -chromo_info / --chromo_info ${A csv file with chromosome info of each feature, default is chromo_list.csv}
    -o / --output_folder ${output folder}
    -m_scs / --num_scs ${total number of scSeq cells}
    -m_fish / --num_fish ${total number of scFISH cells}
    -astral / --astral_directory ${The directory where ASTRAL is installed}
    -fishtree / --fishtree_directory ${The directory where FISHTree is installed}
```

## Simulation
We also provide codes for simulation of both types of data.
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
You can also try ``` python mcmc_newll.py ``` if you set all options to default.

## Reference
[1] X. Fu et al. Joint clustering of single cell sequencing and fluorescencein situ hybridization data to infer tumor copy number phylogenies (in preparation)
[2] S. A. Chowdhury, S. E. Shackney, K. Heselmeyer-Haddad, T. Ried, A. A. Sch ̈affer, and R. Schwartz. Phylogeneticanalysis of multiprobe fluorescence in situ hybridization data from tumor cell populations.  InBioinformatics,volume 29, 2013.
[3] E. M. Gertz, S. A. Chowdhury, W.-J. Lee, D. Wangsa, K. Heselmeyer-Haddad, T. Ried, R. Schwartz, and A. A.Sch ̈affer.  FISHtrees 3.0: Tumor Phylogenetics Using a Ploidy Probe.PLOS ONE, 11(6):e0158569, 6 2016.
[4] C. Zhang, M. Rabiee, E. Sayyari, and S. Mirarab. ASTRAL-III: Polynomial time species tree reconstruction from partially resolved gene trees.BMC Bioinformatics,2018.
