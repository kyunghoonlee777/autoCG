# -*- coding: utf-8 -*-
'''
--ReactionMapper.py--
For solving linear programming(LP)
For mapping atomic connectivity matrices of the reactant and the product
'''

# JCIM 2012, 52, 84.
from gurobipy import *
import numpy as np
import os

def GetEquivalentH(typ, bondlist):
    N = len(typ)
    #************constraint 17'
    visit = np.zeros((N),dtype=int)
    parentH = np.zeros((N),dtype=int)
    child = {}
    for bond in bondlist:
        if typ[bond[0]]==1:
            visit[bond[0]]+=1
            parentH[bond[0]]=bond[1]
        if typ[bond[1]]==1:
            visit[bond[1]]+=1
            parentH[bond[1]]=bond[0]
    for i in range(N):
        if visit[i]==1 and (typ[parentH[i]]!=1 or parentH[i]>i):
            if not(parentH[i] in child): child[parentH[i]] = []
            child[parentH[i]].append(i)
    for a in child:
        if typ[a]==1: child[a].append(a)
    return child


def ReactionMap(type_R,type_P,bondlist_R,bondlist_P,maxSolutions=1,timelimit=100.0):
    """ It finds mapping of atoms between reactant and product.

    type_R, type_P - list of atom elements of reactant and product  (e.g. ['C','O','H','H'])
    bondlist_R, bondlist_P - list of bonds of reactant and product. Each bond is given as a tuple of two atom indices.
    maxSolutions - the number of mapping solutions to be obtained by LP
    timelimit - time limit of LP calculation
    
    returns: indices of permutations (e.g. [1,3,0,2] : Reactant atom 0->Product atom 1, Reactant atom 1->Product atom 3, ...) 
    """
    m = Model("AtomMapping")
    m.setParam(GRB.Param.OutputFlag, 0)
    N = len(type_R)
    if len(type_R)!=len(type_P):
        print ('diff # of atoms!')
        return -1
    y = m.addVars(N, N, name = "y",vtype=GRB.BINARY)
    alpha = m.addVars(len(bondlist_R), len(bondlist_P), name = "alpha",vtype=GRB.BINARY)

    for i in range(N):
        tmp1 = LinExpr()
        tmp2 = LinExpr()
        for j in range(N):
            tmp1.add(y[i,j])
            tmp2.add(y[j,i])
        m.addConstr(tmp1==1) # constraint 2
        m.addConstr(tmp2==1) # constraint 3

    for i in range(N):
        for j in range(N):
            if type_R[i]!=type_P[j]: m.addConstr(y[i,j]==0) # constraint 4

    for a in range(len(bondlist_R)):
        for b in range(len(bondlist_P)):
            i = bondlist_R[a][0]
            j = bondlist_R[a][1]
            k = bondlist_P[b][0]
            l = bondlist_P[b][1]

            m.addConstr(alpha[a,b] <= y[i,k] + y[i,l]) # constraint 5
            m.addConstr(alpha[a,b] <= y[j,k] + y[j,l]) # constraint 6
            m.addConstr(alpha[a,b] <= y[i,k] + y[j,k]) # constraint 11
            m.addConstr(alpha[a,b] <= y[i,l] + y[j,l]) # constraint 12

    #constraint 13~16
    neighbor_R = [[] for i in range(N)]
    neighbor_P = [[] for i in range(N)]

    for bond in bondlist_R:
        t0 = type_R[bond[0]]
        t1 = type_R[bond[1]]
        if not(t1 in neighbor_R[bond[0]]): neighbor_R[bond[0]].append(t1)
        if not(t0 in neighbor_R[bond[1]]): neighbor_R[bond[1]].append(t0)
    for bond in bondlist_P:
        t0 = type_P[bond[0]]
        t1 = type_P[bond[1]]
        if not (t1 in neighbor_P[bond[0]]): neighbor_P[bond[0]].append(t1)
        if not (t0 in neighbor_P[bond[1]]): neighbor_P[bond[1]].append(t0)

    for a in range(len(bondlist_R)):
        i = bondlist_R[a][0]
        j = bondlist_R[a][1]
        tmp = LinExpr()
        for b in range(len(bondlist_P)):
            tmp.add(alpha[a,b])
        tmp1 = LinExpr()
        tmp2 = LinExpr()
        for k in range(N):
            if not(type_R[j] in neighbor_P[k]): tmp1.add(y[i,k])
            if not(type_R[i] in neighbor_P[k]): tmp2.add(y[j,k])
        m.addConstr(tmp + tmp1 <=1) #constraint 13
        m.addConstr(tmp + tmp2 <=1) #constraint 14

    for a in range(len(bondlist_P)):
        k = bondlist_P[a][0]
        l = bondlist_P[a][1]
        tmp = LinExpr()
        for b in range(len(bondlist_R)):
            tmp.add(alpha[b,a])
        tmp1 = LinExpr()
        tmp2 = LinExpr()
        for i in range(N):
            if not (type_P[l] in neighbor_R[i]): tmp1.add(y[i,k])
            if not (type_P[k] in neighbor_R[i]): tmp2.add(y[i,l])
        m.addConstr(tmp + tmp1 <= 1) #15
        m.addConstr(tmp + tmp2 <= 1) #16

    #************constraint 17
    eq_R = GetEquivalentH(type_R, bondlist_R)
    eq_P = GetEquivalentH(type_P, bondlist_P)
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
                            m.addConstr(y[i,l]+y[j,k]<=1) #17'
    """
    #************constraint 19.. 
    # first, find equivalent pair (i, j)
    visit = np.zeros((N),dtype=int) 
    parentH = np.zeros((N),dtype=int)
    child = {}
    for bond in bondlist_R:
        if type_R[bond[0]]==1:
            visit[bond[0]]+=1
            parentH[bond[0]]=bond[1]
        if type_R[bond[1]]==1:
            visit[bond[1]]+=1
            parentH[bond[1]]=bond[0]
    for i in range(N):
        if visit[i]==1 and (type_R[parentH[i]]!=1 or parentH[i]>i):
            if not(parentH[i] in child): child[parentH[i]] = []
            child[parentH[i]].append(i)
    for a in child:
        if type_R[a]==1:
            child[a].append(a)
        for b in range( len(child[a])):
            i = child[a][b]
            for c in range(b+1, len(child[a])):
                j = child[a][c]
                print((i,j))
                tmp1 = LinExpr()
                tmp2 = LinExpr()
                for k in range(N):
                    tmp1.add(k*y[i,k])
                    tmp2.add(k*y[j,k])
                m.addConstr(tmp1 <= tmp2) # 19
    """

    #objective eta, expression 1, # of (breaked bond + formed bond)
    eta = LinExpr()
    for i in range(len(bondlist_R)):
        for j in range(len(bondlist_P)):
            eta.add(alpha[i,j],-2)
    eta.add(len(bondlist_R)+len(bondlist_P))
    m.setObjective(eta, GRB.MINIMIZE)
    m.setParam(GRB.Param.PoolSearchMode, 2) #solution pool
    #m.setParam(GRB.Param.PoolSolutions,maxSolutions+1) # maximum # of solutions
    m.setParam(GRB.Param.PoolSolutions,maxSolutions) # maximum # of solutions
    m.setParam(GRB.Param.PoolGap, 0)
    m.setParam(GRB.Param.TimeLimit, timelimit)
    m.optimize()
    if m.status == GRB.Status.TIME_LIMIT:
        print('Timeout')
        return -1
    if m.status != GRB.Status.OPTIMAL:
        print('no sol')
        return -1

    nSolutions = m.SolCount
    if nSolutions > maxSolutions:
        print ('more than max solutions')
        return -1

    m.setParam(GRB.Param.SolutionNumber,0)
    tmp = np.empty((N),dtype = int)
    for i in range(N):
        for j in range(N):
            #if np.isclose(m.getVarByName('y[{0},{1}]'.format(i,j)).Xn,1): tmp[j]=i
            if np.isclose(m.getVarByName('y[{0},{1}]'.format(i,j)).Xn,1): tmp[i]=j

    return tmp

    '''
    permutation = []
    stri = ''
    print(m.PoolObjVal)
    for k in range(nSolutions):
        m.setParam(GRB.Param.SolutionNumber,k)
        tmp = np.empty((N),dtype = int)
        for i in range(N):
            for j in range(N):
                if np.isclose(m.getVarByName('y[{0},{1}]'.format(i,j)).Xn,1): tmp[j]=i
                
        permutation.append(tmp)
    return permutation
    '''



