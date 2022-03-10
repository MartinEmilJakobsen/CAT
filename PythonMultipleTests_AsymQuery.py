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

Hyps = sys.argv[2+ 2*p*p:]
nHyp = int(len(Hyps)/3)
Hyps = [[Hyps[i+j*3] for i in range(3)] for j in range(nHyp) ]
resList = []

#h = Hyps[0]

for h in Hyps: 
    HypProductMatrix_Low = np.array([[1 for i in range(p)] for j in range(p)])
    HypProductMatrix_High = np.array([[0 for i in range(p)] for j in range(p)])
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
    
    WeightedAdjacency = LowerAdjacency * HypProductMatrix_Low + UpperAdjacency*HypProductMatrix_High
    
    G = nx.from_numpy_matrix(np.array(WeightedAdjacency),create_using=nx.DiGraph)

    #find minimum spanning arborescence via Chu-Liu-Edmonds algorithm
    ed = branchings.minimum_spanning_arborescence(G)
    edges = ed.edges
    
    if hyp == "1":
        if (fr,to) in edges:
            resList.append("0")
        else:
            resList.append("1")
    elif hyp == "-1":
        if (fr,to) in edges:
            resList.append("1")
        else:
            resList.append("0")
    elif hyp == "r":
        if to in [to1 for (fr1,to1) in edges]:
            resList.append("1")
        else:
            resList.append("0")
    
res = ' '.join(resList)

print(res)

