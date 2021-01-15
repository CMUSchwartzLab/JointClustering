#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 14:36:51 2020

Some useful functions

@author: cindyfu
"""

import numpy as np
from pyclustering.cluster.kmedians import kmedians
from scipy.optimize import linear_sum_assignment
from scipy.stats import multivariate_normal
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
import sys
import os
import argparse

def scs_probes(scs, probe_idx, cov_scs=None):
    """
    :param scs: scSeq cluster profiles
    :param probe_idx: probe list for mapping the probe to genome index
    :param cov_scs: covariance of the scSeq cluster profiles, default is None
    :return: the mapped scSeq profiles at FISH probe positions,
            if cov_scs is not None, also return the mapped scSeq profile covariance matrix at FISH probe position
    """
    if not isinstance(probe_idx[0], list):
        scs_probes = scs[:, probe_idx]
    else:
        scs_probes = np.zeros((scs.shape[0], len(probe_idx)))
        for idx in range(len(probe_idx)):
            if isinstance(probe_idx[idx][1], int):
                scs_probes[:, idx] = scs[:,probe_idx[idx][1]]
            else:
                scs_probes[:, idx] = (scs[:, probe_idx[idx][1][0]] + scs[:, probe_idx[idx][1][1]])/2
    if cov_scs is None:
        return scs_probes
    else:
        if not isinstance(probe_idx[0], list):
            cov_scs_probes = cov_scs[:, probe_idx]
        else:
            cov_scs_probes = np.zeros((cov_scs.shape[0], len(probe_idx)))
            for idx in range(len(probe_idx)):
                if isinstance(probe_idx[idx][1], int):
                    cov_scs_probes[:, idx] = cov_scs[:,probe_idx[idx][1]]
                else:
                    cov_scs_probes[:, idx] = (cov_scs[:, probe_idx[idx][1][0]]\
                + cov_scs[:, probe_idx[idx][1][1]])/4
        return scs_probes, cov_scs_probes


########################
# clustering functions #
########################

def kmeans_fish(fish_p1_matrix, k_fish):
    """
    :param fish_p1_matrix: FISH profiles
    :param k_fish: number of clusters
    :return: KMeans cluster centers, cluster labels and frequencies (counts) of each cluster
    """
    cluster = KMeans(n_clusters=k_fish).fit(fish_p1_matrix)
    frac_fish = np.zeros((k_fish))
    for c in cluster.labels_:
        frac_fish[c] += 1
    #frac_fish /= np.sum(frac_fish)
    return cluster.cluster_centers_, cluster.labels_, frac_fish


def kmeans_scs(sc_07_npmatrix, k_scs):
    """
    :param sc_07_npmatrix: scSeq profiles
    :param k_scs: number of clusters
    :return: KMeans cluster centers, cluster labels and frequencies (counts) of each cluster
    """
    cluster = KMeans(n_clusters=k_scs).fit(sc_07_npmatrix)
    frac_scs = np.zeros((k_scs))
    for c in cluster.labels_:
        frac_scs[c] += 1
    #frac_scs /= np.sum(frac_scs)
    return cluster.cluster_centers_, cluster.labels_, frac_scs

def gmm_scs(sc_07_npmatrix, k_scs, cov_type="diag"):
    """
    :param sc_07_npmatrix: scSeq profiles
    :param k_scs: number of clusters
    :param cov_type: covariance type for GMM clustering
    :return: Gaussian mixture cluster centers, cluster labels, frequencies (counts) of each cluster,
    the total loglikelihood and BIC score
    """
    em_scs = GaussianMixture(n_components=k_scs, covariance_type=cov_type, init_params='kmeans')
    em_scs.fit(sc_07_npmatrix)
    mean_ll = em_scs.score(sc_07_npmatrix)/sc_07_npmatrix.shape[1]
    freq_scs = sc_07_npmatrix.shape[0] * em_scs.weights_
    bic_scs = em_scs.bic(sc_07_npmatrix)
    return em_scs.means_, em_scs.covariances_, freq_scs, mean_ll * sc_07_npmatrix.shape[0], bic_scs

def gmm_fish(fish_p1_npmatrix, k_fish, cov_type="diag"):
    """
    :param sc_07_npmatrix: scSeq profiles
    :param k_scs: number of clusters
    :param cov_type: covariance type for GMM clustering
    :return: Gaussian mixture cluster centers, cluster labels, frequencies (counts) of each cluster,
    the total loglikelihood and BIC score
    """
    em_fish = GaussianMixture(n_components=k_fish, covariance_type=cov_type, init_params='kmeans')
    em_fish.fit(fish_p1_npmatrix)
    mean_ll = em_fish.score(fish_p1_npmatrix)/fish_p1_npmatrix.shape[1]
    freq_fish = fish_p1_npmatrix.shape[0] * em_fish.weights_
    bic_fish = em_fish.bic(fish_p1_npmatrix)
    return em_fish.means_, em_fish.covariances_, freq_fish, mean_ll * fish_p1_npmatrix.shape[0], bic_fish



def kmedians_fish(fish_p1_matrix, k_fish):
    """
    KMedians clustering for scFISH
    """
    initial_idx_fish = np.random.choice(fish_p1_matrix.shape[0], size=k_fish)
    initial_median_fish = fish_p1_matrix[initial_idx_fish]
    kmedians_instance_fish = kmedians(fish_p1_matrix, initial_median_fish)

    kmedians_instance_fish.process()
    clusters_fish = kmedians_instance_fish.get_clusters()
    medians_fish = kmedians_instance_fish.get_medians()

    freq_fish = np.zeros((len(clusters_fish)))
    for c in range(len(clusters_fish)):
        freq_fish[c] = len(clusters_fish[c])
    return medians_fish, clusters_fish, freq_fish


def kmedians_scs(sc_07_npmatrix, k_scs):
    """
    KMedians clustering for scSeq
    """
    initial_idx = np.random.choice(sc_07_npmatrix.shape[0], size=k_scs)
    initial_median = sc_07_npmatrix[initial_idx]
    kmedians_instance = kmedians(sc_07_npmatrix, initial_median)

    kmedians_instance.process()
    clusters = kmedians_instance.get_clusters()
    medians = kmedians_instance.get_medians()

    freq = np.zeros((len(clusters)))
    for c in range(len(clusters)):
        freq = len(clusters[c])

    # frac_scs = np.sum(freq, axis=1)/np.sum(freq)
    return medians, clusters, freq

def arrange_components(mat1, mat2):
    assert mat1.shape[0] == mat2.shape[0]
    num_components = mat1.shape[0]
    distmat = np.zeros((num_components, num_components))
    for i in range(num_components):
        for j in range(num_components):
            distmat[i, j] = np.sum(np.abs(mat2[i] - mat1[j])) # the true clones should be in mat2
    _, order = linear_sum_assignment(distmat)
    return order

def main(argv):
    args = get_args(argv)
    sc_npmatrix = np.loadtxt(args["scs"], delimiter=',')
    fish_npmatrix = np.loadtxt(args["fish"], delimiter=',')
    k_scs = args["k_scs"]
    k_fish = args["k_fish"]
    method_scs = args["method_scs"]
    method_fish = args["method_fish"]
    if method_scs == "kmeans":
        means_scs, labels_scs, freqs_scs = kmeans_scs(sc_npmatrix, k_scs)
    elif method_scs == "gmm":
        means_scs, labels_scs, freqs_scs, _, __ = gmm_scs(sc_npmatrix, k_scs)
    else:
        raise Exception("The input method is not available.")

    if method_fish == "kmeans":
        means_fish, labels_fish, freqs_fish = kmeans_fish(fish_npmatrix, k_fish)
    elif method_fish == "gmm":
        means_fish, labels_fish, freqs_fish, _, __ = gmm_fish(fish_npmatrix, k_fish)
    else:
        raise Exception("The input method is not available.")

    # save the outputs to csv files
    np.savetxt("means_scs.csv", means_scs, delimiter=',')
    np.savetxt("means_fish.csv", means_fish, delimiter=',')
    np.savetxt("freqs_scs.csv", freqs_scs, delimiter=',')
    np.savetxt("freqs_fish.csv", freqs_fish, delimiter=',')


def get_args(argv):
    parser = argparse.ArgumentParser(prog='clustering.py', description="Run separate clustering for scSeq and scFISH")
    parser.add_argument('-scs', '--scseq_data', type=str, dest="scs", required=True)
    parser.add_argument('-fish', '--scfish_data', type=str, dest="fish", required=True)
    parser.add_argument('-k_scs', '--clusternum_scseq', type=int, dest="k_scs", required=True)
    parser.add_argument('-k_fish', '--clusternum_fish', type=int, dest="k_fish", required=True)
    parser.add_argument('-method_scs', '--method_scseq', type=str, dest="method_scs", choices=["kmeans", "gmm"], default="kmeans", help="The pre-clustering method for scSeq")
    parser.add_argument('-method_fish', '--method_fish', type=str, dest="method_fish", choices=["kmeans", "gmm"], default="kmeans", help="The pre-clustering method for scFISH")
    return vars(parser.parse_args(argv))


if __name__ == "__main__":
    main(sys.argv[1:])