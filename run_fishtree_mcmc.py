#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 16:58:03 2020

run FISHTree for MCMC results

@author: cindyfu
"""
import numpy as np
from ete3 import Tree
import random
import os
from astral_input2 import Astral
import re
import glob
from collections import deque
import copy
import time
import argparse
import sys

def internal2leaf(parent_child_table, outfile=None):
    """
    Transform the tumor evolutionary tree (with intermediate nodes) to phylogenetic tree (all leaf nodes)

    :param parent_child_table: the parent child table list of the result of FISHTree
     where each item is [parent node, child node, branch length]
    :param outfile: the filename for saving, the newick format of phylogenetic tree will be saved with name outfile
    :return: the string of the phylogenetic tree
    """
    parent_child_dict_nwk = {}
    parent_child_table = np.array(parent_child_table)
    output_dict = {}
    reverse_output = {}
    extra_node_start = 5000
    steiner_node_start = 100
    extra_node = extra_node_start
    pct = np.array(parent_child_table)

    root = 0
    for i in pct[:,0]:
        if isinstance(i, tuple):
            if 0 in i:
                root = i

    node_stack = deque([root])

    while True:
        current_node = node_stack.popleft()
        #print("current node ",current_node)
        children = 0
        for i in range(parent_child_table.shape[0]):
            if parent_child_table[i, 0] == current_node:
                children += 1
                if children == 1 and current_node == root:
                    extra_node += 1
                    current_extra = extra_node
                    output_dict.setdefault(current_extra, [current_node])
                    reverse_output.setdefault(current_node, current_extra)
                    output_dict[current_extra].append(parent_child_table[i, 1])
                    reverse_output.setdefault(parent_child_table[i,1], current_extra)
                    parent_child_dict_nwk.setdefault((current_extra, current_node), 0)
                    parent_child_dict_nwk.setdefault((current_extra, parent_child_table[i,1]), parent_child_table[i,2])

                elif children == 1 and current_node != root:
                    extra_node += 1
                    current_extra = extra_node
                    previous_extra = reverse_output[current_node]
                    output_dict[previous_extra].remove(current_node)
                    output_dict[previous_extra].append(current_extra)
                    reverse_output[current_extra] = previous_extra
                    reverse_output[current_node] = current_extra
                    output_dict.setdefault(current_extra,[current_node, parent_child_table[i, 1]])
                    reverse_output.setdefault(parent_child_table[i,1], current_extra)
                    internal_branch = parent_child_dict_nwk.pop((previous_extra, current_node), None)
                    parent_child_dict_nwk.setdefault((previous_extra, current_extra), internal_branch)
                    parent_child_dict_nwk.setdefault((current_extra, current_node), 0)
                    parent_child_dict_nwk.setdefault((current_extra, parent_child_table[i, 1]), parent_child_table[i, 2])

                elif children >= 2:
                    output_dict[current_extra].append(parent_child_table[i, 1])
                    reverse_output.setdefault(parent_child_table[i,1], current_extra)
                    parent_child_dict_nwk.setdefault((current_extra, parent_child_table[i, 1]), parent_child_table[i, 2])

                node_stack.append(parent_child_table[i, 1])
        if len(node_stack) == 0:
            break

    # generate parent child table of phylogeny tree and output newick format
    parent_child_table_phy = []
    #print("output",output_dict)

    # exclude those nodes with steiner nodes in FISHTree output
    for parent, children in output_dict.items():
        for child in children:
            if not isinstance(child, tuple):
                if int(child) <= extra_node_start and int(child) > steiner_node_start:
                    output_dict[parent].remove(child)

    #split the collapsed nodes two nodes with branch length 0
    output_dict2 = copy.deepcopy(output_dict)
    for parent, children in output_dict2.items():
        for child in children:
            branch_length = parent_child_dict_nwk[(parent, child)]
            if isinstance(child, tuple):
                if branch_length == 0:
                    output_dict[parent].remove(child)
                    del parent_child_dict_nwk[(parent, child)]
                    for c in child:
                        output_dict[parent].append(c)
                        parent_child_dict_nwk.setdefault((parent, c), 0)
                else:
                    extra_node += 1
                    current_extra = extra_node
                    output_dict[parent].remove(child)
                    del parent_child_dict_nwk[(parent, child)]
                    output_dict[parent].append(current_extra)
                    parent_child_dict_nwk.setdefault((parent, current_extra), branch_length)
                    output_dict.setdefault(current_extra, [])
                    for c in child:
                        output_dict[current_extra].append(c)
                        parent_child_dict_nwk.setdefault((current_extra, c), 0)

    # merge the branches after exclusion of steiner nodes, seems useless because it's impossible
    '''
    output_dict2 = output_dict.copy()
    for parent, children in output_dict2.items():
        if len(children) == 1:
            if isinstance(children[0], int) and children[0] > extra_node_start:
                output_dict[parent]=output_dict[children[0]]
                del output_dict[children[0]]
    output_dict2 = output_dict.copy()
    
    for parent, children in output_dict2.items():
        if len(children) == 1:
            if isinstance(children[0], int) and children[0] not in output_dict2.keys():
                for parent2, children2 in output_dict2.items():
                    if parent in children2:
                        output_dict[parent2].append(children[0])
                        output_dict[parent2].remove(parent)
                        del output_dict[parent]
    '''

    # go through final output_dict to obtain parent child table for phylogeny
    for parent, children in output_dict.items():
        for child in children:
            parent_child_table_phy.append((parent, child, parent_child_dict_nwk[(parent, child)]))
    if len(parent_child_table_phy) >= 6:
        #print("pcphy", parent_child_table_phy)
        phylo_tree_nwk = Tree.from_parent_child_table(parent_child_table_phy)
        if outfile is not None:
            phylo_tree_nwk.write(format=9, outfile=outfile)
        print(phylo_tree_nwk.write(format=9))
    return phylo_tree_nwk.write()


class FISHTreeSCS:
    def __init__(self, m_scs, m_fish, scs_chromo, scs_centers, scs_weights, num_probe=8, parallel_num=10):
        """
        :param m_scs: number of scSeq clusters
        :param m_fish: number of scFISH clusters
        :param scs_chromo: chromosomal information of scSeq
        :param scs_centers: scSeq cluster profiles
        :param scs_weights: scSeq cluster weights, the sum of weights equal to 1
        :param num_probe: the number of probes for each subsets after the scSeq genome features are splitted as the input of FISHTree
        :param parallel_num: number of parallel cores for running the FISHTree
        """
        self.scs_centers = scs_centers
        self.scs_weights = scs_weights
        self.scs_fishfmt = (np.round(self.scs_centers)).astype(int) # turn to integers
        self.scs_chromo = scs_chromo
        self.num_sample = m_scs #+ fish_p1_matrix.shape[0]
        self.scs_dim = self.scs_chromo.shape[0]
        self.num_scs_centers, self.total_dim = self.scs_centers.shape 
        # at the beginning and reduced SCS data followed.
        self.num_probe = num_probe
        self.parallel_num = parallel_num
        # unnormalize the weights to get the frequencies (counts)
        self.freq = np.round(self.num_sample * self.scs_weights).astype(int)

        # remove the clusters with frequency of 0
        self.scs_fishfmt = self.scs_fishfmt[np.where(self.freq != 0)[0]]
        self.freq = self.freq[self.freq != 0]
        if self.scs_dim % self.num_probe != 0:
            self.probe_list = [self.num_probe] * (self.scs_dim // self.num_probe) + [self.scs_dim % self.num_probe]
        else:
            self.probe_list = [self.num_probe] * (self.scs_dim // self.num_probe)

    def split_data(self, random_list=None, depth=1, fishtree_dir="~/Documents/softwares/FISHTrees3.1"):
        """
        split the genome features into different subsets randomly if random_list is not given
        create the bash file for running the FISHTree
        :param fishtree_dir: the directory where FISHTree is installed
        :param random_list: the random permutation of genome indices
        :param depth:
        :return: the random_list
        """
        self.rand_list = []
        self.probe_list = self.probe_list * depth
        if random_list is None:
            for d in range(depth):
                rand_list = list(range(self.scs_dim))
                random.shuffle(rand_list)
                self.rand_list += rand_list
        else:
            self.rand_list = random_list

        total_num = np.sum(self.freq)
        root = np.array([[2]] * self.scs_fishfmt.shape[0])

        for p in range(self.parallel_num):
            name='g'+str(p)
            locals()[name] = open("shell"+str(p)+".sh", 'w')
        for k in range(len(self.probe_list)):
            file_g = 'g'+str(k % self.parallel_num)
            scs_fishfmt_split = np.c_[root, self.scs_fishfmt[:, self.rand_list[sum(self.probe_list[:k]): sum(self.probe_list[:k+1])]], self.freq]
            direct = 'FishtreeSCS_' + str(k)

            # print FISHTrees command line
            print("timeout 20m " + fishtree_dir + "/fish " + direct + "/", end=' ', file=locals()[file_g])

            # print splitted SCS reduced centers with FISHTree format
            if not os.path.exists(direct):
                os.makedirs(direct)
            f = open(direct + "/scs_split" + str(k) + ".txt", 'w')
            print("1    diploid " + str(self.probe_list[k]), end='\t', file=f)
            for i in range(self.probe_list[k]):
                print(str(self.rand_list[i + sum(self.probe_list[:k])]), end='\t', file=f)
                print("diploid " + str(self.rand_list[i + sum(self.probe_list[:k])]), end=' ', file=locals()[file_g])
            print("-3" + "\n", end='', file=locals()[file_g])
            print("\n", end='', file=f)
            print(str(self.scs_fishfmt.shape[0]), end='\t' + str(total_num) + '\t', file=f)
            for i in range(self.probe_list[k]):
                print(int(self.scs_chromo[self.rand_list[i + sum(self.probe_list[:k])]]), end='\t', file=f)
            print('\n', end='', file=f)
            for i in range(scs_fishfmt_split.shape[0]):
                if scs_fishfmt_split[i, -1] == 0:
                    continue
                else:
                    for j in range(scs_fishfmt_split.shape[1]):
                        print(scs_fishfmt_split[i, j], end='\t', file=f)
                    print('\n', end='', file=f)
            f.close()
            print("dot -Tpng scs_split" + str(k) + ".ploidyless.dot -O", file=locals()[file_g])
        for p in range(self.parallel_num):
            name='g'+str(p)
            locals()[name].close()
        r = open("random_list", 'w')
        print(self.rand_list, file=r)
        r.close()
        return self.rand_list

    def match_cluster(self, direc="scs_split*.ploidyless.dot", random_list=None):
        """
        Match the cluster profiles from FISHTree output trees back to the whole genome sequencing profiles,
        use indices to perform the matching

        :param direc: the directory name of the FISHTree outputs, inherited from split_data()
        :param random_list: if self.random_list is lost, input the random_list file saved in split_data()
        :return: None
        """
        self.directory = glob.glob(direc)
        if random_list != None:
            self.rand_list = random_list
        steiner = 101
        for dir in range(len(self.directory)):
            match_dict = {}
            dirname = self.directory[dir].split('/')[-1]
            scs_num = int(re.search(r"\d+", dirname).group())
            print("number of scs: ", scs_num)
            astral = Astral(self.directory[dir], self.probe_list, scs_num)
            scs_fishfmt_capped = self.scs_fishfmt.copy()
            scs_fishfmt_capped[self.scs_fishfmt >= 9] = 9
            root_exist = False
            for i in range(astral.fish_centers.shape[0]):
                match = False

                # search for the root in the splitted profiles from FISHTree, generate the mapping dictionary for index before and after
                if np.array_equal(np.array(astral.fish_centers[i, 1:]), np.array([2] * (astral.fish_centers.shape[1]-1))):
                    
                    if i not in match_dict.keys():
                        match_dict.setdefault(i, [0])
                    else:
                        match_dict[i].append(0)
                    
                    match = True
                    root_exist=True

                # match the profiles from FISHTree with the centers, generate the mapping dictionary for index before and after
                for j in range(scs_fishfmt_capped.shape[0]):
                    if np.array_equal(np.array(astral.fish_centers[i, 1:]), scs_fishfmt_capped[
                        j, self.rand_list[sum(self.probe_list[:scs_num]) : sum(self.probe_list[:scs_num]) + self.probe_list[scs_num]]]):
                        if i not in match_dict.keys():
                            match_dict.setdefault(i, [j+1])
                        else:
                            match_dict[i].append(j+1)
                        match = True
                if not match:
                    match_dict.setdefault(i, [steiner])
                    steiner += 1
            '''
            if not root_exist:
                continue
            '''
            try:
                pc_table_after = astral.change_table(match_dict)
                astral.substitute(match_dict, self.directory[dir])
                internal2leaf(pc_table_after, outfile=self.directory[dir] + "_output.nwk")
            except UnboundLocalError:
                continue
            if steiner >= 5000: 
                try:
                    raise RuntimeError('Steiner node index reaches upper bound!')
                except RuntimeError as e:
                    print(e)


def check_if_finish():
    ifFinish = "ps -ef | grep fish | grep -v \'grep\'"
    if os.system(ifFinish) == 256:
        print("FISHTree finished!")
        return True
    else:
        return False


def main(argv):
    args = get_args(argv)
    depth = args["depth"]
    parallel_num = args["cpu_parallel"]
    probe_num = args["num_probe"]
    output_folder = args["output_folder"]
    m_scs = args["m_scs"]
    m_fish = args["m_fish"]
    chromo_list = np.loadtxt(args["chromo_info"], delimiter=',')
    means_scs_ploidy = np.loadtxt(output_folder + "/means_scs_ploidy.csv", delimiter=',')
    weights_scs_ploidy = np.loadtxt(output_folder + "/weights_scs_ploidy.csv", delimiter=',')
    astral_directory = args["astral_directory"]
    fishtree_directory = args["fishtree_directory"]

    directory = output_folder + "/mcmc_fishtree/FishtreeSCS_mcmc"
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.chdir(directory)
    fishtreescs = FISHTreeSCS(m_scs, m_fish, chromo_list,means_scs_ploidy, weights_scs_ploidy / np.sum(weights_scs_ploidy),
                              probe_num, parallel_num)
    fishtreescs.split_data(depth=depth, fishtree_dir=fishtree_directory)
    rand_list = fishtreescs.rand_list
    os.system("chmod u+x shell*.sh")
    for p in range(parallel_num):
        os.system("nohup bash shell" + str(p) + ".sh &")

    while True:
        if check_if_finish():
            break
        else:
            time.sleep(10)

    fishtreescs.match_cluster(direc='scs_split*.ploidyless.dot')

    os.system("sed 's/$/&/g' scs_split*_output.nwk > scs_mcmc_nwks_total.nwk")
    name = "fishtree_astral.nwk"
    os.system(
        "java -jar " + astral_directory + " -i scs_mcmc_nwks_total.nwk -o " + name + " 2>out.log")

    tree_nwk_result = Tree(name)
    for node in tree_nwk_result.iter_descendants():
        if node.name == '0':
            node.delete()
    tree_nwk_result.write(outfile='fishtree_astral_remove0.nwk') # output a tree without the default root node
    print(tree_nwk_result)

def get_args(argv):
    parser = argparse.ArgumentParser(prog='run_fishtree_mcmc.py', description="Run FISHTree for results from joint clustering.")
    parser.add_argument('-d', '--depth', type=int, dest="depth", default=1, help="Sequencing depth for phylogenetic reconstruction")
    parser.add_argument('-p', '--num_probe', type=int, dest="num_probe", default=8, help="Number of probes per subset for each run of FISHTree, maximum=10 for FISHTree software")
    parser.add_argument('-cpu', '--cpu_parellel', type=int, dest="cpu_parallel", default=12)
    parser.add_argument('-chromo_info', '--chromo_info', type=str, dest="chromo_info", default="chromo_list.csv", help="A csv file with chromosome info of each feature")
    parser.add_argument('-o', '--output_folder', type=str, dest="output_folder", default="output")
    parser.add_argument('-m_scs', '--num_scs', type=int, dest="m_scs", default=200)
    parser.add_argument('-m_fish', '--num_fish', type=int, dest="m_fish", default=450)
    parser.add_argument('-astral', '--astral_directory', type=str, dest="astral_directory", default=None, help="The directory where ASTRAL is installed")
    parser.add_argument('-fishtree', '--fishtree_directory', type=str, dest="fishtree_directory", default=None, help="The directory where FISHTree is installed")
    return vars(parser.parse_args(argv))


if __name__ == "__main__":
    main(sys.argv[1:])