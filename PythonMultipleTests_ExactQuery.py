"""
Created on Thu Dec 30 18:06:18 2021

@author: Martin Emil Jakobsen
"""
import sys
import numpy as np
import networkx as nx
import math
from networkx.algorithms.tree import branchings

#sys.argv = ["filename", "3", "1","1","1","1","1","1","1","1","1", "2","2","2","2","2","2","2","2","2","1","2","1","1","2","-1","1","1","r"]

nSysArg = len(sys.argv)

p = int(sys.argv[1])

LowerAdjacency = np.array([[float(sys.argv[2+i+p*j]) for i in range(p)] for j in range(p)])
UpperAdjacency = np.array([[float(sys.argv[2+i+p*j+p*p]) for i in range(p)] for j in range(p)])


G_up = nx.from_numpy_matrix(np.array(UpperAdjacency),create_using=nx.DiGraph)

#find minimum spanning arborescence via Chu-Liu-Edmonds algorithm
ed_up = branchings.minimum_spanning_arborescence(G_up)
score_up = branchings.branching_weight(ed_up)


Hyps = sys.argv[2+ 2*p*p:]
nHyp = int(len(Hyps)/3)
Hyps = [[Hyps[i+j*3] for i in range(3)] for j in range(nHyp) ]
resList = []

for h in Hyps:
    HypProductMatrix = np.array([[1 for i in range(p)] for j in range(p)])
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
    
    LowerAdjacencyMod = LowerAdjacency * HypProductMatrix
    
    G_low = nx.from_numpy_matrix(np.array(LowerAdjacencyMod),create_using=nx.DiGraph)

    #find minimum spanning arborescence via Chu-Liu-Edmonds algorithm
    ed_low = branchings.minimum_spanning_arborescence(G_low)
    score_low = branchings.branching_weight(ed_low)

    if score_low <= score_up:
        resList.append("0")
    else:
        resList.append("1")


res = ' '.join(resList)

print(res)

