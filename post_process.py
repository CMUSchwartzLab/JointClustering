#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 16:10:31 2020

post-process the data after MCMC sampling

@author: cindyfu
"""
import numpy as np


def scs_probes(scs, probe_idx):
    if not isinstance(probe_idx[0], list):
        scs_probes = scs[:, probe_idx]
    else:
        scs_probes = np.zeros((scs.shape[0], len(probe_idx)))
        for idx in range(len(probe_idx)):
            if isinstance(probe_idx[idx][1], int):
                scs_probes[:, idx] = scs[:,probe_idx[idx][1]]
            else:
                scs_probes[:, idx] = (scs[:, probe_idx[idx][1][0]]\
            + scs[:, probe_idx[idx][1][1]])/2
    return scs_probes


def generate_ploidy_matrix(means_scs, means_fish, matrix_freq, probe_idx):
    ploidy_matrix = np.zeros((matrix_freq.shape))
    means_scs_probe = scs_probes(means_scs, probe_idx)
    for i in range(matrix_freq.shape[0]):
        for j in range(matrix_freq.shape[1]):
            if matrix_freq[i, j] >= 0.001:
                ploidy_matrix[i, j] = np.median(2*means_fish[i,:][means_scs_probe[j,:]!=0]
                /means_scs_probe[j,:][means_scs_probe[j,:]!=0])
    return ploidy_matrix


def generate_means_ploidy(means_scs, means_fish, ploidy_matrix, matrix_freq, fish_num, scs_num):
    means_scs_ploidy = []
    weights_scs_ploidy = []
    for scs in range(ploidy_matrix.shape[1]):
        mean_list = []
        for p in range(ploidy_matrix.shape[0]):
            if ploidy_matrix[p, scs] != 0:
                mean_scs_ploidy = np.round(means_scs[scs] * ploidy_matrix[p, scs] / 2).astype(int)
                freq = matrix_freq[p, scs]
                if freq * (fish_num + scs_num) < 2:
                    print("too little weights",p,scs)
                    continue
                same = False
                for m in mean_list:
                    if np.array_equal(mean_scs_ploidy,m):
                        same = True
                        print("same")
                        break
                if not same:
                    means_scs_ploidy.append(mean_scs_ploidy)
                    mean_list.append(mean_scs_ploidy)
                    weights_scs_ploidy.append(freq)
    means_scs_ploidy = np.array(means_scs_ploidy)
    weights_scs_ploidy = np.array(weights_scs_ploidy)
    return means_scs_ploidy, weights_scs_ploidy