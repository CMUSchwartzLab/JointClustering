#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 13:17:18 2019

Simulate the scSeq and scFISH data

@author: cindyfu
"""
import numpy as np
from collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import pickle
import sys

# simulate toy example
def simulate_toy(alpha=10000,beta=10000,k=3,lam=5,chr_num=10,chr_len=500, seg_len_max=200,m_scs=100,m_fish=200):
    p_wgd = {2: 0.3, 1: 0.2, 3: 0.2, 4:0.2, 5: 0.1, 6: 0.1, 7: 0.01, 8: 0.01}
    n = chr_num * chr_len
    info = np.zeros((n))
    for ch in range(chr_num):
        info[ch*chr_len: (ch+1)*chr_len] = ch+1
    # generate the binary tree
    root = [2] * n 

    seg_size_range = np.arange(1,seg_len_max+1)
    lam = lam
    b_rvs = np.random.beta(alpha, beta, size=k)
    u_rvs = np.random.uniform(0, 1, size=k)
    p_cd = 0.2
    tree = {}
    node = 1
    label = {0:[0, 1]}
    mutation = {}  # the corresponding mutated regions and copy number changes in each branch
    leaf_nodes = [0]
    real_scs = {0: np.array(root)}
    ploid_scs = {0: 2}
    for i in range(k-1):
        for j in range(i+1): # traverse all current leaf nodes j
            current_node = leaf_nodes[-j-1]
            left_ = label[current_node][0]
            right_ = label[current_node][1]
            if left_ <= u_rvs[i] <= right_:
                current_scs = real_scs[current_node]
                current_ploid = ploid_scs[current_node]
                # generate left child of the binary tree
                leaf_nodes.pop(-j-1)
                tree.setdefault(node, current_node)
                label.setdefault(node, [left_, left_ + (right_ - left_) * b_rvs[i]])
                leaf_nodes.append(node)
                # mutate from parent
                wgd = False
                if current_ploid in p_wgd.keys():
                    if np.random.uniform(0,1) <= p_wgd[current_ploid]:
                        ploid_scs.setdefault(node, current_ploid * 2)
                        wgd = True
                        mutation.setdefault((current_node, node),(1, ["wgd", 1]))
                    else:
                        ploid_scs.setdefault(node, current_ploid)
                else:
                    ploid_scs.setdefault(node, current_ploid)
                    
                if not wgd:
                    cd = False
                    if np.random.uniform(0,1) <= p_cd:
                        cd = True
                        chromo = np.random.choice(chr_num, 1)[0]
                        cd_cn = np.random.choice([-1,1], 1)[0]
                        mutation.setdefault((current_node, node), (1, ["cd", chromo+1, cd_cn, 1]))
                        left_scs = current_scs.copy().astype(float)
                        left_scs[chromo*chr_len: (chromo+1)*chr_len] += cd_cn / current_ploid * 2
                        real_scs.setdefault(node, left_scs)
                    if not cd:
                        size = 1 #np.random.poisson(lam)
                        seg_size = np.random.choice(seg_size_range, size=size)
                        left_scs = current_scs.copy()
                        
                        loc = np.random.choice(n, size=size)
                        cn = np.random.choice([-1,1], size=size)
                        left_scs = left_scs.astype(float)
                        temp_mutation = []
                        for q, (loci, copynum, segment_size) in enumerate(zip(loc, cn, seg_size)):
                            ch = info[loci]
                            ch_end = np.max(np.where(info == ch)[0])+1
                            mutation_end = np.min(np.array([loci + segment_size, ch_end]))
                            temp_mutation.append(((loci, mutation_end), copynum))
                            left_scs[loci: mutation_end] += copynum / current_ploid * 2
                        left_scs[left_scs <= 0] = 0
                        mutation.setdefault((current_node, node),(size, temp_mutation))
                        #left_scs = np.round(left_scs).astype(int)
                        real_scs.setdefault(node, left_scs)
                else:
                    real_scs.setdefault(node, current_scs)
                    wgd = False
                node += 1
                # generate right child of binary tree
                tree.setdefault(node, current_node)
                label.setdefault(node, [left_ + (right_ - left_) * b_rvs[i], right_])
                leaf_nodes.append(node)
                # mutate from parent
                wgd = False
                if current_ploid in p_wgd.keys():                 
                    if np.random.uniform(0,1) <= p_wgd[current_ploid]:
                        ploid_scs.setdefault(node, current_ploid * 2)
                        wgd = True
                        mutation.setdefault((current_node, node),(1, ["wgd", 1]))
                    else:
                        ploid_scs.setdefault(node, current_ploid)
                else:
                    ploid_scs.setdefault(node, current_ploid)
                    
                if not wgd:
                    cd = False
                    if np.random.uniform(0,1) <= p_cd:
                        cd = True
                        chromo = np.random.choice(chr_num, 1)[0]
                        cd_cn = np.random.choice([-1,1], 1)[0]
                        mutation.setdefault((current_node, node), (1, ["cd", chromo+1, cd_cn, 1]))
                        right_scs = current_scs.copy().astype(float)
                        right_scs[chromo*chr_len: (chromo+1)*chr_len] += cd_cn / current_ploid * 2
                        real_scs.setdefault(node, right_scs)
                    if not cd:
                
                        size = 1 #np.random.poisson(lam)
                        seg_size = np.random.choice(seg_size_range, size=size)
                        right_scs = current_scs.copy()
                        loc = np.random.choice(n, size=size)
                        cn = np.random.choice([-1,1], size=size)
                        right_scs = right_scs.astype(float)
                        temp_mutation = []
                        for q, (loci, copynum, segment_size) in enumerate(zip(loc, cn, seg_size)):
                            ch = info[loci]
                            ch_end = np.max(np.where(info == ch)[0])+1
                            mutation_end = np.min(np.array([loci + segment_size, ch_end]))
                            temp_mutation.append(((loci, mutation_end), copynum))
                            right_scs[loci: mutation_end] += copynum / current_ploid * 2
                        right_scs[right_scs <= 0] = 0
                        mutation.setdefault((current_node, node),(size, temp_mutation))
                        #right_scs = np.round(right_scs).astype(int)
                        real_scs.setdefault(node, right_scs)
                else:
                    real_scs.setdefault(node, current_scs)
                    wgd = False
                node += 1
                break
    print(tree)
    print(ploid_scs)
    simulated_clones = []
    for node, scs in real_scs.items():
        simulated_clones.append(scs)
    simulated_clones = np.array(simulated_clones)
    
    # generate the single cell numbers for Seq and FISH
    #phi_tru = np.random.beta(a=2, b=2, size=2*k-1)
    phi_tru = np.array([1] * (2*k-1))
    while True:
        phi_obs = np.random.dirichlet(phi_tru) # size=1
        cell_scs = np.round(phi_obs * m_scs).astype(int) #np.random.multinomial(m_scs, phi_obs)
        cell_fish = np.round(phi_obs * m_fish).astype(int) #np.random.multinomial(m_fish, phi_obs)
        if np.min(cell_scs) != 0 and np.min(cell_fish) != 0:
            break
    print(cell_scs, cell_fish)
    
    ploid_list = np.zeros((2*k-1))
    ploid_count = 0
    for clone, ploid in ploid_scs.items():
        ploid_list[ploid_count] = ploid
        ploid_count += 1
    print(ploid_list)
    real_clones=(simulated_clones.T * ploid_list/2).T
    mutation_plus1 = {}
    for key in mutation.keys():
        mutation_plus1[tuple(np.array(key) + 1)] = mutation[key]
    print(mutation_plus1)
    return simulated_clones, ploid_list, cell_scs/m_scs, cell_fish/m_fish, tree, real_clones, mutation_plus1, info


# main function for simulation
def simulate(perturb_scs, perturb_fish, sc_07_npmatrix, probe_idx, chromo_list, alpha = 10000, beta = 10000, k=9, lam=10,
             cn_rand=False, n=None, wgd_simu=True,sigma2=0.3, perturb_method='int', perturb_sigma=0.3, m_fish=450, m_scs=None):
    # simulate data from tree
    # phi_tru is the same, phi_obs is the observed freqs
    if sc_07_npmatrix is not None:
        sc_07_npmatrix = sc_07_npmatrix.astype(int)
        if m_scs is None:
            m_scs = sc_07_npmatrix.shape[0] # cell number
        if n is None:
            n = sc_07_npmatrix.shape[1] # dimension
        rate = Counter(sc_07_npmatrix.reshape(sc_07_npmatrix.shape[0]*sc_07_npmatrix.shape[1]).tolist())
        p = (sc_07_npmatrix != 2).sum(axis=0) / (sc_07_npmatrix != 2).sum()
    else:
        if m_scs is None:
            raise Exception("Number of scSeq cells should not be None!")
        if n is None:
            n = 9934
        if chromo_list is None:
            chromo_list_index = [0, 791, 1649, 2359, 3037, 3674, 4282, 4820, 5336, 5730, 6193, 6660, 7127, 7478,
                           7793, 8073, 8336, 8599, 8872, 9045, 9259, 9382, 9495, 9900, 9934] # obtained from real GBM data
            chromo_list = []
            for chromo in range(1, 25):
                for i in range(chromo_list_index[chromo-1], chromo_list_index[chromo]):
                    chromo_list.append(chromo)
            chromo_list = np.array(chromo_list)
            np.savetxt("chromo_list.csv", chromo_list, delimiter=',')
        rate = {2: 1823873, 1: 160262, 3: 60571, 4: 14429, 0: 3164,5.0: 3045,6.0: 835,7.0: 62, 8: 21, 9: 5, 10: 5}  # calculated from real GBM data
        p = np.loadtxt("loc_p.csv", delimiter='\t')              # calculated from real GBM data

    rate_p = np.zeros((11))
    for cn, freq in rate.items():
        rate_p[cn] = freq
    rate_p = np.delete(rate_p,2)
    rate_p /= np.sum(rate_p) 
    if wgd_simu:
        p_wgd = {2: 0.2, 4:0.2, 8: 0.01}
    else:
        p_wgd = {2: 0.0}
    p_cd = 0.3

    mutation = {}
    chr_num = np.max(chromo_list)
    
    
    # generate the binary tree
    root = [2] * n 
    
    cn_range = [0,1,3,4,5,6,7,8,9,10]
    seg_size_range = np.arange(10,110)
    seg_size_range = np.append(seg_size_range, np.arange(250,350))
    seg_size_range = np.append(seg_size_range, np.arange(580,590))
    lam = lam
    b_rvs = np.random.beta(alpha, beta, size=k)
    u_rvs = np.random.uniform(0, 1, size=k)
    tree = {}
    node = 1
    label = {0:[0, 1]}
    leaf_nodes = [0]
    real_scs = {0: np.array(root)}
    ploid_scs = {0: 2}
    
    for i in range(k-1):
        for j in range(i+1): # traverse all current leaf nodes j
            current_node = leaf_nodes[-j-1]
            left_ = label[current_node][0]
            right_ = label[current_node][1]
            if left_ <= u_rvs[i] <= right_:
                current_scs = real_scs[current_node]
                current_ploid = ploid_scs[current_node]
                # generate left child of the binary tree
                leaf_nodes.pop(-j-1)
                tree.setdefault(node, current_node)
                label.setdefault(node, [left_, left_ + (right_ - left_) * b_rvs[i]])
                leaf_nodes.append(node)
                # mutate from parent
                wgd = False
                if current_ploid in p_wgd.keys():
                    if np.random.uniform(0,1) <= p_wgd[current_ploid]:
                        ploid_scs.setdefault(node, current_ploid * 2)
                        mutation.setdefault((current_node, node),(1, [("wgd", 1)]))
                        wgd = True
                    else:
                        ploid_scs.setdefault(node, current_ploid)
                else:
                    ploid_scs.setdefault(node, current_ploid)
                if not wgd:
                    temp_mutation = []
                    cd = False
                    if np.random.uniform(0,1) <= p_cd:
                        cd = True
                        chromo = np.random.choice(np.arange(1, chr_num+1), 1)[0]
                        cd_cn = np.random.choice([-1,1], 1)[0]
                        ch_end = np.max(np.where(chromo_list == chromo)[0])+1
                        ch_start = np.min(np.where(chromo_list == chromo)[0])
                        temp_mutation.append(("cd", chromo, cd_cn, ch_end-ch_start))
                        left_scs = current_scs.copy().astype(float)
                        left_scs[ch_start: ch_end] += cd_cn
                        real_scs.setdefault(node, left_scs)
                    size = np.random.poisson(lam)
                    seg_size = np.random.choice(seg_size_range, size=size)
                    left_scs = current_scs.copy()
                    if cn_rand:
                        loc = np.random.choice(n, size=size, p=p)
                        cn = np.random.choice(cn_range, size=size, p=rate_p)
                        for q, (loci, copynum, segment_size) in enumerate(zip(loc, cn, seg_size)):
                            left_scs[loci: loci + segment_size] = copynum
                    else:   
                        loc = np.random.choice(n, size=size)
                        cn = np.random.choice([-1,1], size=size)
                        left_scs = left_scs.astype(float)
                        temp_mutation = []
                        for q, (loci, copynum, segment_size) in enumerate(zip(loc, cn, seg_size)):
                            ch = chromo_list[loci]
                            ch_end = np.max(np.where(chromo_list == ch)[0])+1
                            mutation_end = np.min(np.array([loci + segment_size, ch_end]))
                            temp_mutation.append(((loci, mutation_end), copynum, mutation_end - loci))
                            left_scs[loci: mutation_end] += copynum
                        left_scs[left_scs <= 0] = 0
                        mutation.setdefault((current_node, node),(size+cd, temp_mutation))
                    real_scs.setdefault(node, left_scs)
                else:
                    real_scs.setdefault(node, current_scs)
                    wgd = False
                node += 1
                # generate right child of binary tree
                tree.setdefault(node, current_node)
                label.setdefault(node, [left_ + (right_ - left_) * b_rvs[i], right_])
                leaf_nodes.append(node)
                # mutate from parent
                wgd = False
                if current_ploid in p_wgd.keys():                 
                    if np.random.uniform(0,1) <= p_wgd[current_ploid]:
                        ploid_scs.setdefault(node, current_ploid * 2)
                        wgd = True
                    else:
                        ploid_scs.setdefault(node, current_ploid)
                else:
                    ploid_scs.setdefault(node, current_ploid)
                if not wgd:
                    cd = False
                    temp_mutation = []
                    if np.random.uniform(0,1) <= p_cd:
                        cd = True
                        chromo = np.random.choice(np.arange(1, chr_num+1), 1)[0]
                        cd_cn = np.random.choice([-1,1], 1)[0]
                        left_scs = current_scs.copy().astype(float)
                        ch_end = np.max(np.where(chromo_list == chromo)[0])+1
                        ch_start = np.min(np.where(chromo_list == chromo)[0])
                        temp_mutation.append(("cd", chromo, cd_cn, ch_end-ch_start))
                        left_scs[ch_start: ch_end] += cd_cn
                        real_scs.setdefault(node, left_scs)
                    size = np.random.poisson(lam)
                    seg_size = np.random.choice(seg_size_range, size=size)
                    right_scs = current_scs.copy()
                    if cn_rand:
                        loc = np.random.choice(n, size=size, p=p)
                        cn = np.random.choice(cn_range, size=size, p=rate_p)
                        for q, (loci, copynum, segment_size) in enumerate(zip(loc, cn, seg_size)):
                            right_scs[loci: loci + segment_size] = copynum
                    else:
                        loc = np.random.choice(n, size=size)
                        cn = np.random.choice([-1,1], size=size)
                        right_scs = right_scs.astype(float)
                        temp_mutation = []
                        for q, (loci, copynum, segment_size) in enumerate(zip(loc, cn, seg_size)):
                            ch = chromo_list[loci]
                            ch_end = np.max(np.where(chromo_list == ch)[0])+1
                            mutation_end = np.min(np.array([loci + segment_size, ch_end]))
                            temp_mutation.append(((loci, mutation_end), copynum, mutation_end - loci))
                            right_scs[loci: mutation_end] += copynum
                        right_scs[right_scs <= 0] = 0
                        mutation.setdefault((current_node, node),(size + cd, temp_mutation))
                    real_scs.setdefault(node, right_scs)
                else:
                    real_scs.setdefault(node, current_scs)
                    wgd = False
                node += 1
                break
    simulated_clones = []
    for node, scs in real_scs.items():
        simulated_clones.append(scs)
    simulated_clones = np.array(simulated_clones)
    
    # generate the single cell numbers for Seq and FISH
    #phi_tru = np.random.beta(a=2, b=2, size=2*k-1)
    phi_tru = np.array([1] * (2*k-1))
    while True:
        phi_obs = np.random.dirichlet(phi_tru) # size=1
        cell_scs = np.round(phi_obs * m_scs).astype(int)
        cell_fish = np.round(phi_obs * m_fish).astype(int)
        if np.min(cell_scs) != 0 and np.min(cell_fish) != 0:
            break
    count_actual_cluster = 0
    
    # generate single cell Sequencing from the clones
    for cluster in range(len(cell_scs)):
        count_actual_cluster += 1
        current_clone_scs = np.tile(simulated_clones[cluster], (cell_scs[cluster],1))
        if count_actual_cluster > 1:
            simulated_scs = np.append(simulated_scs, current_clone_scs, axis=0)
        else:
            simulated_scs = current_clone_scs.copy()  # first clone condition
    # extract FISH probes
    simulated_clones_fish = np.zeros((simulated_clones.shape[0], len(probe_idx)))
    for idx in range(len(probe_idx)):
        if isinstance(probe_idx[idx][1], int):
            simulated_clones_fish[:, idx] = simulated_clones[:,probe_idx[idx][1]]
        else:
            simulated_clones_fish[:, idx] = (simulated_clones[:, probe_idx[idx][1][0]]
            + simulated_clones[:, probe_idx[idx][1][1]])/2

    ploid_list = np.zeros((len(cell_scs)))
    ploid_count = 0
    for clone, ploid in ploid_scs.items():
        ploid_list[ploid_count] = ploid
        ploid_count += 1
    print(ploid_list)
    for clone_num in range(len(ploid_list)):
        simulated_clones_fish[clone_num, :] *= ploid_list[clone_num]/2
    count_fish = 0
    for cluster in range(len(cell_fish)):
        count_fish += 1
        current_clone_fish = np.tile(simulated_clones_fish[cluster], (cell_fish[cluster],1))
        if count_fish > 1:
            simulated_fish = np.append(simulated_fish, current_clone_fish, axis=0)
        else:
            simulated_fish = current_clone_fish.copy()
            
    # create gt_mat (ground truth matrix A) and gt_01_mat (ground truth matrix Z)
    gt_mat = []
    gt_01_mat = []
    num_clones_fish = list(range(simulated_clones_fish.shape[0]))
    for i in range(simulated_clones_fish.shape[0]):
        if i not in num_clones_fish:
            continue
        add = np.zeros((simulated_clones.shape[0]))        
        add2 = np.zeros((simulated_clones.shape[0]+1))

        add[i] = phi_obs[i]
        add2[-1] += cell_fish[i]/m_fish
        add2[i] = 1
        for j in range(i+1,simulated_clones_fish.shape[0]):
            if np.array_equal(simulated_clones_fish[i], simulated_clones_fish[j]):
                add[j] = phi_obs[j]
                add2[-1] += cell_fish[j]/m_fish
                add2[j] = 1
                num_clones_fish.remove(j)
        gt_mat.append(add)
        gt_01_mat.append(add2)
    gt_01_mat.append(np.append(cell_scs/m_scs, np.array([0])))
    gt_mat = np.array(gt_mat)
    gt_01_mat = np.array(gt_01_mat)
    simulated_clones_fish_compress = simulated_clones_fish[num_clones_fish]     

    # add fish system noise ~ N(0,sigma2)
    simulated_clones_fish_compress += np.round(np.random.normal(0, sigma2, size=simulated_clones_fish_compress.shape))
    
    remove_idx = []
    for i in range(gt_01_mat.shape[1]-1):
        for j in range(i+1, gt_01_mat.shape[1]-1):
            if np.array_equal(simulated_clones[i], simulated_clones[j]):
                gt_mat[:,i] += gt_mat[:,j]
                gt_01_mat[:,i] += gt_01_mat[:,j]
                remove_idx.append(j)
                print(j)
    simulated_clones_compress = np.copy(simulated_clones)
    simulated_clones_compress = np.delete(simulated_clones_compress, remove_idx, axis=0)
    gt_mat = np.delete(gt_mat, remove_idx, axis=1)
    gt_01_mat = np.delete(gt_01_mat, remove_idx, axis=1)
    gt_mat_freq = np.c_[gt_mat, np.tile(gt_01_mat[:-1, -1],[1,1]).T]
    gt_mat_freq = np.r_[gt_mat_freq, np.tile(gt_01_mat[-1, :],[1,1])]
    plt.figure()
    sns.heatmap(gt_mat_freq, annot=True, cmap="Reds")
    plt.figure()
    sns.heatmap(gt_01_mat, annot=True, cmap="Reds")

    # perturb scSeq    
    if perturb_scs is not None:
        perturb_frac = perturb_scs
        for scs in range(np.sum(cell_scs)):
            perturb_vec = np.zeros((n))
            if perturb_method == 'int':
                perturb_size = np.random.poisson(perturb_frac * n)
                perturb_loc = np.random.choice(n, size=perturb_size)
                perturb_vec[perturb_loc] = np.random.choice([-1,1], size=perturb_size)
            elif perturb_method == 'gaussian':
                perturb_vec = np.random.normal(0, perturb_sigma, size=sc_07_npmatrix.shape[1])
            simulated_scs[scs] += np.round(perturb_vec).astype(int)

    if perturb_fish is not None:
        # perturb scFISH
        perturb_fish_num = perturb_fish
        for fish in range(np.sum(cell_fish)):
            perturb_vec = np.zeros((len(probe_idx)))
            if perturb_method == 'int':
                perturb_size_fish = np.random.poisson(perturb_fish_num)
                perturb_probe = np.random.choice(len(probe_idx), size=perturb_size_fish)
                perturb_vec[perturb_probe] = np.random.choice([-1,1], size=perturb_size_fish)
            elif perturb_method == 'gaussian':
                perturb_vec = np.random.normal(0, perturb_sigma, size=len(probe_idx))
            simulated_fish[fish] += np.round(perturb_vec).astype(int)
    
    simulated_fish[simulated_fish < 0] = 0
    simulated_scs[simulated_scs < 0] = 0
    simulated_fish = np.round(simulated_fish).astype(int)
    simulated_scs = np.round(simulated_scs).astype(int)
    mutation_plus1 = {}
    if not cn_rand: 
        for key in mutation.keys():
            mutation_plus1[tuple(np.array(key) + 1)] = mutation[key]
        print(mutation_plus1)
    return simulated_clones, ploid_list, simulated_scs, simulated_fish,\
cell_scs/m_scs, cell_fish/m_fish, tree, gt_mat, gt_01_mat,simulated_clones_compress, \
    simulated_clones_fish_compress,phi_obs, mutation_plus1


def main(argv):
    args = get_args(argv)
    alpha = 10 ** np.random.uniform(-4, 4)
    beta = alpha
    perturb_scs = args["perturb_scs"]
    perturb_fish = args["perturb_fish"]
    if args["scs"] is None:
        sc_07_npmatrix=None
    else:
        sc_07_npmatrix=np.loadtxt(args["scs"], delimiter="\t")
    with open(args["probe_idx"], 'rb') as f:
        probe_idx = pickle.load(f)
    chromo_list = args["chromo_info"]
    k = args["num_leaves"]
    lam = args["lambda"]
    sigma2 = args["sigma2"]
    n = args["num_features"]
    wgd = args["wgd"]
    perturb_method = args["perturb_method"]
    if perturb_method == 'gaussian':
        perturb_sigma2 = args["perturb_sigma2"]
    else:
        perturb_sigma2 = None
    m_fish = args["m_fish"]
    m_scs = args["m_scs"]
    simulated_clones, ploid_list, simulated_scs, simulated_fish, freqs_scs, freqs_fish, tree, gt_mat, gt_01_mat, \
    simulated_clones_compress, simulated_clones_fish_compress, phi_obs, mutation_plus1 = simulate(perturb_scs,
    perturb_fish, sc_07_npmatrix, probe_idx, chromo_list, alpha, beta, k, lam,
             False, n, wgd, sigma2, perturb_method, perturb_sigma2, m_fish, m_scs)
    np.savetxt("scs_simu.csv", simulated_scs, delimiter=',')
    np.savetxt("fish_simu.csv", simulated_fish, delimiter=',')


def get_args(argv):
    parser = argparse.ArgumentParser(prog='simulate_data.py', description="Generate simulation data")
    parser.add_argument('-perturb_scs', '--perturb_scs', type=float, dest="perturb_scs", default=0.05, help="Perturbation rate of scSeq from the cluster centers when method is int")
    parser.add_argument('-perturb_fish', '--perturb_fish', type=float, dest="perturb_fish", default=0.05, help="Perturbation rate of scFISH from the cluster centers when method is int")
    parser.add_argument('-scs', '--scseq_data', type=str, dest="scs", default=None)
    parser.add_argument('-chromo_info', '--chromo_info', type=str, dest="chromo_info", default=None, help="A csv file with chromosome info of each feature")
    parser.add_argument('-probe_idx', '--probe_index_list', type=str, dest="probe_idx", default="probe_idx.pickle", help="A list which shows the genome indices of scSeq data corresponding to the FISH probes")
    parser.add_argument('-k', '--num_leaves', type=int, dest="num_leaves", default=5)
    parser.add_argument('-lam', '--lambda', type=int, dest="lambda", default=10, help="The mutation rate")
    parser.add_argument('-sigma2', '--sigma2', type=float, dest="sigma2", default=0.2, help="The noise between the ground truth scFISH marker copy numbers and corresponding scSeq features")
    parser.add_argument('-n', '--num_features', type=int, dest="num_features", default=None)
    parser.add_argument('-wgd', '--wgd', action="store_true", dest="wgd", default=True)
    parser.add_argument('-perturb_method', '--perturb_method', choices=['int', 'gaussian'], dest="perturb_method", default='int', help="Perturbation method")
    parser.add_argument('-perturb_sigma2', '--perturb_sigma2', dest="perturb_sigma2", default=0.3, help="Perturbation sigma2 when method is gaussian")
    parser.add_argument('-m_scs', '--num_scs', type=int, dest="m_scs", default=200)
    parser.add_argument('-m_fish', '--num_fish', type=int, dest="m_fish", default=450)
    return vars(parser.parse_args(argv))


if __name__ == '__main__':
    main(sys.argv[1:])
