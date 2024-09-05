from scipy.optimize import milp
from scipy.optimize import LinearConstraint
from scipy.optimize import Bounds
from scipy.sparse import bsr_matrix
from scipy.sparse import bsr_array

import numpy as np
import sys
import os
import copy
# import time

def GetEquivalentH(atom, bond):
    N = len(atom)
    visit = []
    parent = []
    for i in range(N):
        visit.append(0)
        parent.append(0)
    child = {}
    for a, b in bond:
        if atom[a] == 'H':
            visit[a] += 1
            parent[a] = b
        if atom[b] == 'H':
            visit[b] += 1
            parent[b] = a
    for i in range(N):
        if visit[i] == 1 and (atom[parent[i]] != 'H' or parent[i] > i):
            if not (parent[i] in child): child[parent[i]] = []
            child[parent[i]].append(i)
    for a in child:
        if atom[a] == 'H': child[a].append(a)
    return child

def ReactionMap(R_atom, P_atom, R_bond, P_bond, maxSol=1, printRes=False):
    # start_time = time.perf_counter()
    # printRes=True to print full results
    num_R_atom = len(R_atom)
    num_P_atom = len(P_atom)
    num_R_bond = len(R_bond)
    num_P_bond = len(P_bond)
    num_P_p_atom = num_R_atom - num_P_atom

    # P_p denotes P', made by connecting any remaining atoms
    # from R (not included in P)
    P_p_atom = list(copy.deepcopy(R_atom))
    for atom in P_atom:
        P_p_atom.remove(atom)

    # P_p_bond is fully connected, such that the reaction map results in 
    # minimal R bond breaking 
    P_p_bond = list()
    for i in range(num_P_atom, num_R_atom):
        for j in range(i+1, num_R_atom):
            P_p_bond.append((i, j))

    P_atom_all = P_atom + P_p_atom
    P_bond_all = P_bond + P_p_bond
    num_P_atom_all= len(P_atom_all)
    num_P_bond_all= len(P_bond_all)
            
    # Number of variables to optimize
    N = num_R_atom*num_R_atom + num_R_bond*num_P_bond_all
        
    P_graph = np.zeros((num_P_atom, num_P_atom))
    for a, b in P_bond:
        P_graph[a, b] = P_graph[b, a] = 1

    A = []
    b_lb = []
    b_ub = []

    # mapping between atoms of 
    # R and P + P' should be one-to-one
    for i in range(num_R_atom):
        tmp = np.zeros(N)
        for j in range(num_R_atom):
            tmp[i*num_R_atom+j] = 1
        # Constraint 2
        A.append(tmp)
        b_lb.append(1)
        b_ub.append(1)

    for i in range(num_R_atom):
        tmp = np.zeros(N)
        for j in range(num_R_atom):
            tmp[j*num_R_atom+i] = 1
        # Constraint 3
        A.append(tmp)
        b_lb.append(1)
        b_ub.append(1)

    tmp = np.zeros(N)
    for i in range(num_R_atom):
        for j in range(num_P_atom_all):
            if R_atom[i] != P_atom_all[j]: tmp[i*num_R_atom+j] = 1
    # Constraint 4
    A.append(tmp)
    b_lb.append(0)
    b_ub.append(0)

    for a in range(num_R_bond):
        for b in range(num_P_bond_all):
            i, j = R_bond[a]
            k, l = P_bond_all[b]
            # Constraint 5
            tmp = np.zeros(N)
            tmp[i*num_R_atom+k] = tmp[i*num_R_atom+l] = -1
            tmp[num_R_atom*num_R_atom+a*num_P_bond_all+b] = 1
            A.append(tmp)
            b_lb.append(-1)
            b_ub.append(0)
            # Constraint 6
            tmp = np.zeros(N)
            tmp[j*num_R_atom+k] = tmp[j*num_R_atom+l] = -1
            tmp[num_R_atom*num_R_atom+a*num_P_bond_all+b] = 1
            A.append(tmp)
            b_lb.append(-1)
            b_ub.append(0)
            # Constraint 11
            tmp = np.zeros(N)
            tmp[i*num_R_atom+k] = tmp[j*num_R_atom+k] = -1
            tmp[num_R_atom*num_R_atom+a*num_P_bond_all+b] = 1
            A.append(tmp)
            b_lb.append(-1)
            b_ub.append(0)
            # Constraint 12
            tmp = np.zeros(N)
            tmp[i*num_R_atom+l] = tmp[j*num_R_atom+l] = -1
            tmp[num_R_atom*num_R_atom+a*num_P_bond_all+b] = 1
            A.append(tmp)
            b_lb.append(-1)
            b_ub.append(0)
   
     
    R_adj = [[] for _ in range(num_R_atom)]
    P_adj = [[] for _ in range(num_R_atom)]
    for bond in R_bond:
        t0, t1 = R_atom[bond[0]], R_atom[bond[1]]
        if not (t1 in R_adj[bond[0]]): R_adj[bond[0]].append(t1)
        if not (t0 in R_adj[bond[1]]): R_adj[bond[1]].append(t0) 
    for bond in P_bond_all:
        t0, t1 = P_atom_all[bond[0]], P_atom_all[bond[1]]                  
        if not (t1 in P_adj[bond[0]]): P_adj[bond[0]].append(t1)
        if not (t0 in P_adj[bond[1]]): P_adj[bond[1]].append(t0)
    for a in range(num_R_bond):
        i, j = R_bond[a]
        tmp1 = np.zeros(N)
        tmp2 = np.zeros(N)
        for b in range(num_P_bond_all):
            tmp1[num_R_atom*num_R_atom+a*num_P_bond_all+b] = tmp2[num_R_atom*num_R_atom+a*num_P_bond_all+b] = 1
        for k in range(num_P_atom_all):
            if not (R_atom[j] in P_adj[k]): tmp1[i*num_R_atom+k] = 1
            if not (R_atom[i] in P_adj[k]): tmp2[j*num_R_atom+k] = 1
        #Constraint 13
        A.append(tmp1)
        b_lb.append(0)
        b_ub.append(1)
        #Constraint 14
        A.append(tmp2)
        b_lb.append(0)
        b_ub.append(1)
    for a in range(num_P_bond):
        k, l = P_bond[a]
        tmp1 = np.zeros(N)
        tmp2 = np.zeros(N)
        for b in range(num_R_bond):
            tmp1[num_R_atom*num_R_atom+b*num_P_bond_all+a] = tmp2[num_R_atom*num_R_atom+b*num_P_bond_all+a] = 1
        for i in range(num_R_atom):
            if not (P_atom[l] in R_adj[i]): tmp1[i*num_R_atom+k] = 1
            if not (P_atom[k] in R_adj[i]): tmp2[i*num_R_atom+l] = 1
        #Constraint 15
        A.append(tmp1)
        b_lb.append(0)
        b_ub.append(1)
        #Constraint 16
        A.append(tmp2)
        b_lb.append(0)
        b_ub.append(1)
     
    #Constraint 17
    eq_R = GetEquivalentH(R_atom, R_bond)
    eq_P = GetEquivalentH(P_atom_all, P_bond_all)
    for a in eq_R:
        for b in range(len(eq_R[a])):
            i = eq_R[a][b]
            for c in range(b+1, len(eq_R[a])):
                j = eq_R[a][c]
                for d in eq_P:
                    for e in range(len(eq_P[d])):
                        k = eq_P[d][e]
                        for f in range(e+1, len(eq_P[d])):
                            l = eq_P[d][f]
                            tmp = np.zeros(N)
                            tmp[i*num_R_atom+l] = tmp[j*num_R_atom+k] = 1
                            A.append(tmp)
                            b_lb.append(0)
                            b_ub.append(1)

    # print("constraints", time.perf_counter() - start_time)
    # start_time = time.perf_counter()
     
    # Optimization
    eta = np.zeros(N) # Objective function
    for i in range(num_R_bond):
        for j in range(num_P_bond):
            eta[num_R_atom*num_R_atom+i*num_P_bond_all+j] = -2
        # Fully connected P'
        for j in range(num_P_bond, num_P_bond_all):
            eta[num_R_atom*num_R_atom+i*num_P_bond_all+j] = -1
    integrality = np.ones(N)
    A = bsr_matrix(A)
    #A = np.array(A)
    b_lb = np.array(b_lb)
    b_ub = np.array(b_ub)
    constraints = LinearConstraint(A,b_lb,b_ub)
    bounds = Bounds(np.zeros(N), np.ones(N))
    result = milp(c=eta, constraints=constraints, bounds=bounds, integrality=integrality, options=dict(disp=False))
    if printRes:
        print(result)
        # for i in range(num_R_bond):
        #     for j in range(num_P_bond):
        #         if np.isclose(result.x[num_R_atom*num_R_atom+i*num_P_bond_all+j], 1):
        #             print(str(R_bond[i])+' >> '+str(P_bond[j]))
    mapping = np.full(num_R_atom, -1)
    for i in range(num_R_atom):
        for j in range(num_P_atom):
            if np.isclose(result.x[i*num_R_atom+j], 1): mapping[i] = j

    # print("milp", time.perf_counter() - start_time)

    return mapping



