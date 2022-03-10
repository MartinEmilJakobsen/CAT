# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 18:06:18 2021

@author: Martin Emil Jakobsen
"""
import sys
import numpy as np
import networkx as nx
import math
from networkx.algorithms.tree import branchings

# Input Format: Filename QueryScheme p LowerAdjacency UpperAdjacency Hypothesis
# weighted adjacency matrices are printed with rowwise space separated elements
# Hypothesis space separated restrictions of the form "from to restriction", where
# restriction is 1 for present edge, -1 for missing edge, r for root (here from=to=root).

#sys.argv = ["filename","exact", "3", "1","1","1","1","1","1","1","1","1", "2","2","2","2","2","2","2","2","2","1","2","1","1","2","-1","1","1","r"]

#Query scheme
queryscheme = sys.argv[1]
#SystemSize
p = int(sys.argv[2])
#Adjacency matrices
LowerAdjacency = np.array([[float(sys.argv[3+i+p*j]) for i in range(p)] for j in range(p)])
UpperAdjacency = np.array([[float(sys.argv[3+i+p*j+p*p]) for i in range(p)] for j in range(p)])
#Hypothesis / substructure restriction to test
Hyp = sys.argv[3+ 2*p*p:]
nRestrictions = int(len(Hyp)/3)
Hyp = [[Hyp[i+j*3] for i in range(3)] for j in range(nRestrictions) ]


if queryscheme == "asymptotic":
    HypProductMatrix_Low = np.array([[1 for i in range(p)] for j in range(p)])
    HypProductMatrix_High = np.array([[0 for i in range(p)] for j in range(p)])
    for h in Hyp: 
        fr = int(h[0])-1
        to = int(h[1])-1
        hyp = h[2]
        if hyp == "1":
            HypProductMatrix_Low[:,to] = 0
            HypProductMatrix_Low[to,fr] = 0
            HypProductMatrix_Low[fr,to] = 1
            HypProductMatrix_High[:,to] = 1
            HypProductMatrix_High[to,fr] = 1
            HypProductMatrix_High[fr,to] = 0
        elif hyp == "-1":
            HypProductMatrix_Low[fr,to] = 0
            HypProductMatrix_High[fr,to] = 1
        elif hyp == "r":
            HypProductMatrix_Low[:,to] = 0
            HypProductMatrix_High[:,to] = 1
    #calculate the weighted adjacency   
    WeightedAdjacency = LowerAdjacency * HypProductMatrix_Low + UpperAdjacency*HypProductMatrix_High
    #creaty fully connected graph with edge weights corresponding to the weightedadjacency    
    G = nx.from_numpy_matrix(np.array(WeightedAdjacency),create_using=nx.DiGraph)
    #find MWDST of the above graph via Chu-Liu-Edmonds algorithm
    ed = branchings.minimum_spanning_arborescence(G)
    #pull the edges of the MWDST
    edges = ed.edges
    
    output = 1
    
    for h in Hyp:
        fr = int(h[0])-1
        to = int(h[1])-1
        hyp = h[2]
        if hyp == "1":
            if (fr,to) in edges:
                continue
            else:
                break
        elif hyp == "-1":
            if (fr,to) in edges:
                break
            else:
                continue
        elif hyp == "r":
            if to in [to1 for (fr1,to1) in edges]:
                break
            else:
                continue
    else:
        output = 0
    
    print(output)
    
elif queryscheme == "exact":
    #fully connected graph with upper weights
    G_up = nx.from_numpy_matrix(np.array(UpperAdjacency),create_using=nx.DiGraph)
    #MWDST of upper weights via Chu-Liu-Edmonds algorithm
    ed_up = branchings.minimum_spanning_arborescence(G_up)
    #score of upper weight MWDST
    score_up = branchings.branching_weight(ed_up)
    
    HypProductMatrix = np.array([[1 for i in range(p)] for j in range(p)])
    for h in Hyp:
        fr = int(h[0])-1
        to = int(h[1])-1
        hyp = h[2]
        if hyp == "1":
            HypProductMatrix[:,to] = 0
            HypProductMatrix[fr,to] = 1
            HypProductMatrix[to,fr] = 0
        elif hyp == "-1":
            HypProductMatrix[fr,to] = 0
        elif hyp == "r":
            HypProductMatrix[:,to] = 0
    #The graph with lower score and restricted to satisfy hyp    
    LowerAdjacencyMod = LowerAdjacency * HypProductMatrix    
    G_low = nx.from_numpy_matrix(np.array(LowerAdjacencyMod),create_using=nx.DiGraph)
    #find MWDST of lower restricted graph via Chu-Liu-Edmonds algorithm
    ed_low = branchings.minimum_spanning_arborescence(G_low)
    #score of lower restricted
    score_low = branchings.branching_weight(ed_low)
    
    
    if score_low <= score_up:
        output = 0
    else:
        output = 1
    print(output)
else:
    print("Error queryscheme must be set to either asymptotic or exact")