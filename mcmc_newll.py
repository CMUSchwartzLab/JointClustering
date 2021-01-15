#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 16:41:12 2020

Choose nonzero matrix by MCMC samplers based on the likelihood l(A,Z),

calculate A matrix by optimizating the likelihood l(A)

@author: cindyfu

"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.special import erf, erfc
import argparse
import sys
import os

from clustering import kmeans_scs, kmeans_fish, gmm_scs, scs_probes
from scipy.optimize import minimize, differential_evolution, Bounds, LinearConstraint
from post_process import generate_means_ploidy, generate_ploidy_matrix
import pickle


class MCMC_sampler2:
    def __init__(self, means_fish, means_scs, frac_fish, frac_scs, probe_idx, sigma=1):
        """
            The MCMC sampler for optimizing the likelihood of the matching of scSeq and scFISH clusters

         :param means_fish - preclustered scFISH cluster profiles
         :param means_scs - preclustered scSeq cluster profiles
         :param frac_fish - preclustered scFISH cluster frequencies (counts)
         :param frac_scs - preclustered scSeq cluster frequencies (counts)
         :param probe_idx - a list which shows the genome indices of scSeq data corresponding to the FISH probes,
                  i.e. [[0,1000],[1, (2000,2001)]] means the scFISH probe with index 0 corresponds to the scSeq profile at the index of 1000,
                  the scFISH probe with index 1 corresponds to the scSeq profile at the indices of both 2000 and 2001

        """
        self.means_fish = np.array(means_fish)
        self.means_scs = np.array(means_scs)
        self.frac_fish = frac_fish
        self.frac_scs = frac_scs
        self.probe_idx = probe_idx
        self.m = len(means_fish)
        self.n = len(means_scs)
        self.error_matrix = self.gaussian_error(sigma)
        self.proposal = self.normalize(self.error_matrix)
        self.ll_list = []
        self.ll_ind_list = []
        self.matrix01_list = []
        self.matrix_list = []
        self.reject_rate_list = []
        self.count_valid_iter = 0
        self.reject_count = 0
        self.accept_count = 0
    
    def initialize(self, matrix01_init=None):
        if matrix01_init is not None:
            self.matrix01 = matrix01_init
        else:
            self.matrix01 = np.zeros((self.m, self.n))
            self.matrix01[np.where(np.random.uniform(0, 10^(-9),size=(self.m, self.n)) < self.proposal)] = 1
        self.matrix01_ = self.matrix01
        res = self.solve_ll()
        self.matrix_freq = self.reshape_res(res.x)
        self.ll = self.calculate_ll(self.matrix_freq)
        self.ll_ind =  self.calculate_indicator_prob(self.matrix01)
        return self.ll + self.ll_ind            

    def load_model(self, path_model):
        load_dict = pickle.load(open(path_model, 'rb'))
        self.ll_list = load_dict['ll_list']
        self.ll_ind_list = load_dict['ll_ind_list']
        self.matrix_list = load_dict['matrix_list']
        self.count_valid_iter = load_dict['count_valid_iter']
        self.reject_rate_list = load_dict['reject_rate_list']
        self.matrix01_list = load_dict['matrix01_list']
        self.matrix01 = load_dict['matrix01']
        self.matrix_freq = load_dict['matrix_freq']
        self.reject_count = load_dict['reject_count']
        self.accept_count = load_dict['accept_count']
        self.ll = self.ll_list[-1]
        self.ll_ind = self.ll_ind_list[-1]

    def sampling(self, min_iter, max_iter, reject_threshold=0.93):
        """

        :param min_iter: the minimum number of iteration
        :param max_iter: the maximum number of iteration
        :param reject_threshold: if the iteration is between min_iter and max_iter, if the rejection rate
        is larger than reject_threshold for 5 consecutive record of rejection rate, early stop the training
        :return: the maximum log likelihood l(A,Z)
        """
        while True:
            posterior_proposal = np.abs(self.matrix01 - self.proposal)

            # normalize the probability of matching for the proposal
            move_matrix = np.reshape(np.random.multinomial(1, posterior_proposal.flatten()/\
                                    np.sum(posterior_proposal)), (self.m, self.n))
            move_spot = np.where(move_matrix == 1)
            self.matrix01_ = np.abs(self.matrix01 - move_matrix)

            # make sure there is at least one matching for every column and row, otherwise skip and redo the random move
            if np.min(np.sum(self.matrix01_, 0)) == 0 or np.min(np.sum(self.matrix01_, 1)) == 0:
                continue

            res = self.solve_ll()
            self.matrix_freq_ = self.reshape_res(res.x)
            self.ll_ = self.calculate_ll(self.matrix_freq_)
            self.ll_ind_ = self.calculate_indicator_prob(self.matrix01_)
        
            self.count_valid_iter += 1
            u = np.random.uniform(0, 1)
            posterior_proposal_ = np.abs(self.matrix01_ - self.proposal)
            
            prob_post = posterior_proposal[move_spot][0]/np.sum(posterior_proposal)
            prob_post_ = posterior_proposal_[move_spot][0]/np.sum(posterior_proposal_)
            ratio = np.exp(self.ll_ - self.ll + self.ll_ind_ - self.ll_ind) * prob_post_ / prob_post
            if  u < min(1, ratio):
                self.matrix01 = self.matrix01_
                self.matrix_freq = self.matrix_freq_
                self.ll = self.ll_
                self.ll_ind = self.ll_ind_
                self.matrix01_list.append(self.matrix01_)
                self.matrix_list.append(self.matrix_freq_)
                self.ll_list.append(self.ll_)
                self.ll_ind_list.append(self.ll_ + self.ll_ind_)
                self.accept_count += 1
            else:
                self.ll_list.append(self.ll)
                self.ll_ind_list.append(self.ll + self.ll_ind)
                self.matrix01_list.append(self.matrix01)
                self.matrix_list.append(self.matrix_freq)
                self.reject_count += 1
            rr = self.reject_count / (self.reject_count + self.accept_count)
            if self.count_valid_iter % 500 == 0:
                self.reject_rate_list.append(rr)
                idx = np.argmax(self.ll_ind_list)
                print("iteration: ", self.count_valid_iter, ", current best likelihood: l(A):", self.ll_list[idx], ", l(A,Z):", self.ll_ind_list[idx])
                if self.count_valid_iter >= min_iter:
                    if np.sum(np.array(self.reject_rate_list[-5:-1]) >= reject_threshold) == 5:
                        break
            if self.count_valid_iter >= max_iter:
                break
        return np.max(self.ll_ind_list)
            
    def save_model(self, save_model_dir):
        pickle_out = open(save_model_dir, 'wb')
        index = np.argmax(np.array(self.ll_ind_list))
        self.best_matrix_freq = self.matrix_list[index]
        self.best_ll_ind = self.ll_ind_list[index]
        save_dict = {'ll_list': self.ll_list, 'll_ind_list': self.ll_ind_list,\
                     'matrix_list': self.matrix_list, 'best_matrix_freq': self.best_matrix_freq ,\
                         'best_ll_ind': self.best_ll_ind, 'count_valid_iter': self.count_valid_iter,\
                        'reject_rate_list': self.reject_rate_list, 'matrix01_list': self.matrix01_list,\
                            'matrix01': self.matrix01, 'matrix_freq': self.matrix_freq,\
                            'reject_count': self.reject_count, 'accept_count': self.accept_count}
        pickle.dump(save_dict, pickle_out)
        pickle_out.close()
    
    def reshape_res(self, res_x):
        matrix_freq = np.zeros((self.m, self.n))
        matrix_freq[self.index_x, self.index_y] = res_x
        return matrix_freq

    def solve_ll(self):
        """
        the function for solving the optimization in the M-H method
        :return: the resolution of the matrix A given Z
        """
        self.index_x, self.index_y = np.where(self.matrix01_ == 1)
        num_var = len(self.index_x)
        bounds = Bounds(np.zeros((num_var)), np.ones(num_var))
        linear_constraint = LinearConstraint(np.ones(num_var)[np.newaxis,:], [1], [1])
        x0 = np.abs(np.random.normal(0,1,num_var))
        x0 /= np.sum(x0)
        neg_likelihood = Neg_likelihood(self.frac_fish, self.frac_scs, self.index_x, self.index_y)
        res = minimize(neg_likelihood.neg_ll, x0, constraints=[linear_constraint], bounds=bounds)
        return res

    def calculate_ll(self, matrix_freq):
        ll = np.sum(self.frac_fish * np.log(np.sum(matrix_freq, axis=1)))\
    + np.sum(self.frac_scs * np.log(np.sum(matrix_freq, axis=0)))
        return ll
    
    def calculate_indicator_prob(self, matrix01):
        log_ind_prob = np.sum(np.log(self.error_matrix[np.where(matrix01 == 1)]))
        return log_ind_prob

    def calculate_l(self, matrix_freq):
        l = np.prod((np.sum(matrix_freq, axis=1))**self.frac_fish)\
    * np.prod((np.sum(matrix_freq, axis=0))**self.frac_scs)
        return l
    
    def gaussian_error(self, sigma):
        self.means_scs_probe = scs_probes(self.means_scs, self.probe_idx)
        ploidy_range = list(range(1, 9))
        error_matrix = np.zeros((self.m, self.n))
        for j in range(self.means_scs.shape[0]):
            for i in range(self.means_fish.shape[0]):
                for ploidy in ploidy_range:
                    p = np.prod(erfc((np.abs(self.means_fish[i,:] / ploidy * 2 - self.means_scs_probe[j,:]))/sigma))
                    if p > error_matrix[i,j]:
                        error_matrix[i,j] = p
        return error_matrix

    @staticmethod
    def normalize(error_matrix):  
        vertical_prob = error_matrix / np.sum(error_matrix, axis=0)[np.newaxis, :]
        horizontal_prob = error_matrix / np.sum(error_matrix, axis=1)[:, np.newaxis]
        prob_matrix_normalized = np.sqrt(vertical_prob * horizontal_prob)
        return prob_matrix_normalized


class Neg_likelihood:
    def __init__(self, frac_fish, frac_scs, index_x, index_y):
        """
            Negative likelihood class, used for optimization
        """
        self.frac_fish = frac_fish
        self.frac_scs = frac_scs
        self.m = len(frac_fish)
        self.n = len(frac_scs)
        self.index_x = index_x
        self.index_y = index_y
        
    def neg_l(self, x):
        # for rows
        l = 1
        for i in range(self.m):
            l *= np.sum(x[self.index_x == i])**self.frac_fish[i]
        # for columns
        for j in range(self.n):
            l *= np.sum(x[self.index_y == j])**self.frac_scs[j]
        return -l
    
    def neg_ll(self, x):
        # for rows
        ll = 0
        for i in range(self.m):
            ll += self.frac_fish[i] * np.log(np.sum(x[self.index_x == i]))
        # for columns
        for j in range(self.n):
            ll += self.frac_scs[j] * np.log(np.sum(x[self.index_y == j]))
        return -ll

def plot_results_real(matrix_list, ll_list, ll_ind_list, freq_fish, freq_scs):
    """
    :param matrix_list: the resulting matrix list from the sampling
    :param ll_list: the resulting likelihood l(A) list from the sampling
    :param ll_ind_list: the resulting likelihood l(Z) list from the sampling
    :param freq_fish: frequencies of scFISH from preclustering
    :param freq_scs: frequencies of scSeq from preclustering
    :return: visualization of resulting matrix using seaborn package
    """
    index = np.argmax(np.array(ll_ind_list))
    print(index, ll_list[index], ll_ind_list[index])
    matrix_freq = matrix_list[index]
    matrix_freq = np.c_[matrix_freq, np.tile(freq_fish/np.sum(freq_fish),[1,1]).T]
    matrix_freq = np.r_[matrix_freq, np.tile(np.append(freq_scs/np.sum(freq_scs),np.array([0])),[1,1])]
    plt.figure()
    sns.heatmap(matrix_freq, annot=True,cmap='Blues')
    plt.title("Estimated on real data")
    plt.xlabel("scSeq cluster index")
    plt.ylabel("scFISH cluster index")


def run_mcmc(means_fish, means_scs, freqs_fish, freqs_scs, probe_idx, iterations,
             sigma, trial, dataset, bound_iter, bound_ll, directory,save=True):
    """
        the function for running the algorithm
    ```
    :param means_scs: preclustered scSeq cluster profiles
    :param means_fish: preclustered scFISH cluster profiles
    :param freqs_fish: preclustered scFISH cluster frequencies
    :param freqs_scs: preclustered scSeq cluster frequencies
    :param probe_idx: a list which shows the genome indices of scSeq data corresponding to the FISH probes
    :param iterations: the number of iterations
    :param sigma: scaling factor of gaussian error function
    :param trial: the number of trial
    :param dataset: the dataset name, used for saving
    :param bound_iter: the number of iteration used for bound step
    :param bound_ll: the likelihood threshold for bound step
    :param directory: the directory for saved model
    :param save: the boolean, True if you want to save the model
    :return: the MCMC_sampler2 instance
    """
    while True:    
        mcmc = MCMC_sampler2(means_fish, means_scs,freqs_fish, freqs_scs, probe_idx, sigma)
        mcmc.initialize()
        ll_ind_max = mcmc.sampling(min_iter=bound_iter,max_iter=bound_iter)
        if ll_ind_max > bound_ll:
            break
    ll_ind_max = mcmc.sampling(min_iter=iterations,max_iter=iterations)
    if save:
        if not os.path.exists(directory):
            os.mkdir(directory)
        mcmc.save_model(save_model_dir= directory + '/' + dataset + str(iterations) + '_' + str(sigma) + '_' + str(trial) + '.pickle')
    return mcmc


def main(argv):
    args = get_args(argv)
    means_fish = np.loadtxt(args["means_fish"], delimiter=',')
    means_scs = np.loadtxt(args["means_scs"], delimiter=',')
    freqs_fish = np.loadtxt(args["freqs_fish"], delimiter=',')
    freqs_scs = np.loadtxt(args["freqs_scs"], delimiter=',')
    output_folder = args["output_folder"]
    sigma = args["sigma"]
    trial = args["trial"]
    iterations = args["max_iter"]
    bound_ll = args["bound_ll"]
    bound_iter = args["bound_iter"]
    dataset = args["dataset"]
    with open(args["probe_idx"], 'rb') as f:
        probe_idx = pickle.load(f)

    mcmc_result = run_mcmc(means_fish, means_scs, freqs_fish, freqs_scs, probe_idx, iterations,
             sigma, trial, dataset, bound_iter, bound_ll, directory=output_folder, save=True)

    ploidy_matrix = generate_ploidy_matrix(means_scs, means_fish, mcmc_result.best_matrix_freq, probe_idx)
    means_scs_ploidy, weights_scs_ploidy = generate_means_ploidy(means_scs, means_fish, ploidy_matrix,
                                                mcmc_result.best_matrix_freq, np.sum(freqs_scs), np.sum(freqs_fish))
    np.savetxt(output_folder + "/means_scs_ploidy.csv",means_scs_ploidy, delimiter=',')
    np.savetxt(output_folder + "/weights_scs_ploidy.csv", weights_scs_ploidy, delimiter=',')


def get_args(argv):
    parser = argparse.ArgumentParser(prog='mcmc_newll.py', description="Run joint matching of scFISH and scSeq clusters based on MCMC sampling")
    parser.add_argument('-sp', '--scseq_profiles', type=str, dest="means_scs", default="means_scs.csv")
    parser.add_argument('-fp', '--scfish_profiles', type=str, dest="means_fish", default="means_fish.csv")
    parser.add_argument('-sf', '--scseq_freqs', type=str, dest="freqs_scs", default="freqs_scs.csv")
    parser.add_argument('-ff', '--scfish_freqs', type=str, dest="freqs_fish", default="freqs_fish.csv")
    parser.add_argument('-o', '--output_folder', type=str, dest="output_folder", default="output")
    parser.add_argument('-sigma', '--sigma', type=float, dest="sigma", default=1.0, help="The scaling factor for gaussian error distribution")
    parser.add_argument('-trial', '--trial', type=int, dest="trial", default=1, help="The number of trial for training, used for naming")
    parser.add_argument('-i', '--iterations', type=int, dest="max_iter", default=50000)
    parser.add_argument('-bound_ll', '--bound_ll', type=float, dest="bound_ll", default=-10**7, help="Bound log likelihood for branching process in the beginning of training")
    parser.add_argument('-bound_iter', '--bound_iter', type=int, dest="bound_iter", default=0, help="Bound iteration for branching process in the beginning of training")
    parser.add_argument('-dataset', '--dataset', type=str, dest="dataset", default="simulation", help="Used for naming.")
    parser.add_argument('-probe_idx', '--probe_index_list', type=str, dest="probe_idx", default="probe_idx.pickle", help="A list which shows the genome indices of scSeq data corresponding to the FISH probes")
    return vars(parser.parse_args(argv))


if __name__ == '__main__':
    main(sys.argv[1:])
