#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 15:53:28 2020

Class Astral, used for applying ASTRAL for FISHTree output.

@author: cindyfu
"""


from graphviz import Source,render
import os
import sys
import io as inout
import re
import numpy as np
from ete3 import Tree
from collections import deque
import glob
import random

#sys.stdout = inout.TextIOWrapper(sys.stdout.buffer,encoding='utf8')  # change the default encode of standard io

class Astral:
    def __init__(self, directory, probe_list, scs_num):
        self.probe_list = probe_list
        self.scs_num = scs_num
        self.fish_centers, self.pc_table_origin = self.read_in(directory)

    def read_in(self, directory):
        file = open(directory, encoding="utf8")  # transform each fishtree to utf8 encoded
        line = file.readline()
        line = line.strip()
        node_dict = {}
        fish_centers = []   # a matrix with node index at 0 and fish copy numbers at 1-8
        parent_child_table = []
        while len(line) > 1:
            line = file.readline()
            line = line.strip()
            match = re.search(r'\d+\s\[label', line)
            if match != None:
                node = re.search(r'\d+', match.group())
                if self.probe_list[self.scs_num] == 8:
                    fish = re.match(r'(\d+).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\.\d+', line).groups()
                elif self.probe_list[self.scs_num] == 6:
                    fish = re.match(r'(\d+).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\.\d+', line).groups()
                elif self.probe_list[self.scs_num] == 10:
                    fish = re.match(r'(\d+).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\.\d+', line).groups()
                elif self.probe_list[self.scs_num] == 12:
                    fish = re.match(r'(\d+).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\((\d).*\.\d+', line).groups()
                elif self.probe_list[self.scs_num] == 4:
                    fish = re.match(r'(\d+).*\((\d).*\((\d).*\((\d).*\((\d).*\.\d+', line).groups()
                elif self.probe_list[self.scs_num] == 2:
                    fish = re.match(r'(\d+).*\((\d).*\((\d).*\.\d+', line).groups()
                # fish = re.match(r'(\d+\s\[label=\))\"\{PDGFA\((\d)\(\\nAPC\((\d)\)\\nEGFR\((\d)\)\\nMET\((\d)\)\\nMYC\((\d)\)\\nCCND1\((\d)\)\\nCHEK1\((\d)\)\\nERG\((\d)\)\|(.*)\}\",style=(.*)', line).group()
                fish_centers.append(np.array(fish).astype(int))
                node_dict.setdefault(node, )
            else:
                match = re.search(r'\d+\s->\s\d+.*', line)
                if match != None:
                    pc_pair = re.match(r'(\d+)\s->\s(\d+).*(\d)\.',line).groups()
                    parent_child_table.append(np.array(list(pc_pair)).astype(int))
        fish_centers = np.array(fish_centers)
        parent_child_table = np.array(parent_child_table).astype(int)
        file.close()
        return fish_centers, parent_child_table

    def change_table(self, map_dict):
        """
        Change the index of parent children table to real index of centers
        :param map_dict: the mapping dictionary of the indices of nodes before and after applying FISHTree
        :return: the changed parent children table
        """
        parent_child_table_after = []
        for i in range(self.pc_table_origin.shape[0]):
            if len(map_dict[self.pc_table_origin[i, 0]]) > 1:
                a = tuple(map_dict[self.pc_table_origin[i, 0]])
            else:
                a = map_dict[self.pc_table_origin[i, 0]][0]
            if len(map_dict[self.pc_table_origin[i, 1]]) > 1:
                b = tuple(map_dict[self.pc_table_origin[i, 1]])
            else:
                b = map_dict[self.pc_table_origin[i, 1]][0]
            parent_child_table_after.append([a, b, self.pc_table_origin[i, 2]])
        parent_child_table_after = np.array(parent_child_table_after)
        return parent_child_table_after
    
    def substitute(self, map_dict, directory):
        """
        Generate the dot file and the corresponding tree png images after replacing the FISHTree profiles with the real node indices
        :param map_dict: mapping dictionary
        :param directory: the directory for saving
        :return: None
        """
        file = open(directory, encoding="utf8")  # transform each fishtree to utf8 encoded
        line = file.readline()
        line = line.strip()
        new_directory = directory.split('.')[0] + '_sub.' + directory.split('.')[-1]
        g = open(new_directory, 'w')
        print(line, file=g, end='\n')
        while len(line) > 1:
            line = file.readline()
            line = line.strip()
            match = re.search(r'\d+\s\[label', line)
            if match != None:
                node = re.search(r'\d+', match.group()).group()
                new_line = re.sub(r'{', r'{'+str(map_dict[int(node)])+'\\n', line)
                print(new_line, file=g, end='\n')
            else:
                print(line, file=g, end='\n')
        file.close()
        g.close()
        os.system("dot -Tpng " + new_directory + " -O")


# if __name__ == "__main__":
    # Astral.internal2leaf(np.array([[0,1,1],[0,2,1],[1,3,1],[1,4,1]]), outfile=None)
