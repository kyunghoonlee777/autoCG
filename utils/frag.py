# -*- coding: utf-8 -*-
'''
--frag.py--
Analyzing/calculating various properties of molecules

It contains routines that
 - Read a xyz file and save the geometry
 - Find bonds in molecules, make adjacency matrices and bond order matrices
 - Distinguish molecules by defining Coulomb matrices
and so on..
'''

import os
import subprocess
import numpy as np

from scipy import spatial

from autoCG.utils import am
from autoCG import chem

visit=[] 
visited_V=[]


def makeAdjacency(atom_num,Blist): 
    """ It calculates an adjacency matrix.

    atom_num - the number of atoms in a molecule
    Blist - bondlist obtained from frag.findbond (list of bonds in a molecule. Each bond is given as a tuple of two atom indices.)

    returns: adjacency matrix
    """
    Adj=np.zeros((atom_num,atom_num))
    for i in range(len(Blist)):
        atom1=Blist[i][0]; atom2=Blist[i][1]
        Adj[atom1][atom2]=1; Adj[atom2][atom1]=1
    return Adj


def getDistanceMatrix(Adj,GetPath=None):
    """ It calculates a distance matrix from the adjacency matrix, using Floyd-Warshall algorithm.

    Adj - adjacency matrix
    GetPath - None -> do not return a path
              [u,v] -> return a path from u to v

    returns: Distance matrix, and optionally, a path from u to v 
    """
    D=1000*np.ones((len(Adj),len(Adj)))
    Next=1000*np.ones((len(Adj),len(Adj)))
    for i in range(len(D)): D[i][i]=0
    for i in range(len(D)):
        for j in range(i+1,len(D)):
            if Adj[i][j]!=0: D[i][j]=D[j][i]=1; Next[i][j]=j; Next[j][i]=i
    for k in range(len(D)):
        for i in range(len(D)):
            for j in range(len(D)):
                if D[i][j]>D[i][k]+D[k][j]: D[i][j]=D[i][k]+D[k][j]; Next[i][j]=Next[i][k]
            #for j in range(i+1,len(D)):
                #if D[i][j]>D[i][k]+D[k][j]: D[i][j]=D[j][i]=D[i][k]+D[k][j]; Next[i][j]=Next[i][k]

    if GetPath==None: return D
    else:
        u=GetPath[0]; v=GetPath[1]
        Path=[u]
        while u!=v:
            u=int(Next[u][v])
            Path.append(u)
        return D,Path

def get_trunc_atom_list(BO,unsatuatom_list,Kthneighbor):
    """ It returns atom indices of active atoms and non-active ones placed at the K-th nearest neighboring position of active ones.

    BO - bond order (or adjacency) matrix
    unsatuatom_list - list of atom indices of active atoms
    Kthneighbor - Given as an integer value (>=1). In terms of molecular graph,
                  non-active atoms placed at the K-th nearest position with respect to active atoms are only considered during solving MILP

    returns: truncated_atom_list - atom indices of active atoms and non-active ones placed at the K-th nearest neighboring position of active ones.
    """
    truncated_atom_list=unsatuatom_list[:]
    D=getDistanceMatrix(BO)
    #Gather the indices of atoms which are the Kth neighbors of reactive atoms (default K=1)
    for unsatuA in unsatuatom_list:
        for a in range(len(D)):
            if D[unsatuA][a]<=Kthneighbor and truncated_atom_list.count(a)<1: truncated_atom_list.append(a)
    truncated_atom_list.sort()
    return truncated_atom_list

def setDelocalizedBonds(atom_list,BO):
    """ It analyzes bond order matrix and finds bonds whose order can be between one and two (1.5) or two and three (2.5)

    atom_list - list of atom class instances
    BO - bond order matrix 

    returns:
        NewBO - a new bond order matrix with delocalized bond orders applied
        DelocalizedBondCount - the number of delocalized bonds
    """
    Delocalized_Atoms=[]
    for i in range(len(atom_list)):
        if str.upper(atom_list[i].element)=='C' and np.count_nonzero(BO[i])==3 and int(np.sum(BO[i]))==4: Delocalized_Atoms.append(i)

    Delocalized_Bonds=[]
    for i in range(len(Delocalized_Atoms)):
        for j in range(i+1,len(Delocalized_Atoms)):
            atom1=Delocalized_Atoms[i]; atom2=Delocalized_Atoms[j]
            if BO[atom1][atom2]!=0: Delocalized_Bonds.append(set([atom1,atom2]))

    NewBO=np.copy(BO)
    DelocalizedBondCount=0
    Real_Delocalized_Bonds=[]
    for i in range(len(Delocalized_Bonds)):
        for j in range(len(Delocalized_Bonds)):
            if i==j: continue
            Bond1=Delocalized_Bonds[i]; Bond2=Delocalized_Bonds[j]
            if len(Bond1&Bond2)==1: 
                atom1=list(Bond1)[0]; atom2=list(Bond1)[1]
                NewBO[atom1][atom2]=NewBO[atom2][atom1]=1.5
                DelocalizedBondCount+=1
                Real_Delocalized_Bonds.append([atom1,atom2])
                break
    ## Handling (3.5,3.5) or (4.5,4.5) or (3.5,4.5), (4.5,3.5) Bonds - 170501
    Real_Deloc_Atoms=[]
    for deloc_bond in Real_Delocalized_Bonds: 
        Real_Deloc_Atoms.append(deloc_bond[0]); Real_Deloc_Atoms.append(deloc_bond[1])
    Real_Deloc_Atoms=list(set(Real_Deloc_Atoms))
    BOSums=[np.sum(NewBO[x]) for x in Real_Deloc_Atoms]

    for deloc_bond in Real_Delocalized_Bonds: 
        atom1=deloc_bond[0]; atom2=deloc_bond[1]
        if '%.1f' % BOSums[Real_Deloc_Atoms.index(atom1)]!='4.0' and '%.1f' % BOSums[Real_Deloc_Atoms.index(atom2)]!='4.0':
            NewBO[atom1][atom2]=NewBO[atom2][atom1]=BO[atom1][atom2]

    return NewBO, DelocalizedBondCount


def BOmatrixtoList(adj):
    """ It converts bond order matrix (or adjacency matrix) into adjacency list
    adj - bond order (or adjacency) matrix

    return: Adjlist - adjacency list. list of bonds in a molecule. Each bond is given as a tuple of two atom indices.
    """
    Adjlist=[]
    for i in range(len(adj)):
        for j in range(i+1,len(adj)):
            if adj[i][j]!=0: Adjlist.append((i,j))
    return Adjlist


def FindResonances(atom_list,Adj,Unsatu_set,DU,Total_charge):
    """ It finds all possible resonances (bond order matrices) for the given molecule

    atom_list - list of atom class instances
    Adj - adjacency matrix
    Unsatu_set - list of indices of unsaturated atoms to which multiple bond orders can be assigned.
    DU - list of degrees of unsaturation of atoms
    Total_charge - total charge of a molecule

    returns: AllBOs - list of bond order matrices of possible resonances
    """
    AllBOs=[]
    for start_point in range(len(Unsatu_set)):
        BOs=[np.copy(Adj)]
        DUs=[DU[:]]; 
        unsatu_set=Unsatu_set[:]
        unsatu_set.insert(0,unsatu_set.pop(start_point)) #change start point
        Unsatu_sets=[unsatu_set[:]]

        while True:
            NewBOs=[]; NewDUs=[]; NewUns=[]
            for i in range(len(BOs)):
                Curr_BO=BOs[i]; Curr_DU=DUs[i]; Curr_Un=Unsatu_sets[i]
                while True:
                    if Curr_Un==[]: break

                    now=Curr_Un[0]
                    num_neighbors,neighbors=Find_Neighbors(Curr_BO,now,Curr_Un)

                    if num_neighbors==0: Curr_Un.pop(0) #no neighbors to make multiple bonds with
                    else: break

                    if Curr_Un==[]: break
                if Curr_Un==[]: 
                    AllBOs.append(Curr_BO)
                    continue

                if num_neighbors!=[]:
                    for target in neighbors:
                        NewBO=np.copy(Curr_BO)
                        NewDU=Curr_DU[:]
                        NewUn=Curr_Un[:]

                        Additional_bond=min(Curr_DU[now],Curr_DU[target])

                        NewBO[now][target]+=Additional_bond
                        NewBO[target][now]+=Additional_bond
                        NewDU[now]-=Additional_bond; NewDU[target]-=Additional_bond

                        if NewDU[now]<=0: NewUn.pop(NewUn.index(now))
                        if NewDU[target]<=0: NewUn.pop(NewUn.index(target))
                        NewBOs.append(NewBO); NewDUs.append(NewDU); NewUns.append(NewUn)
            BOs=NewBOs; DUs=NewDUs; Unsatu_sets=NewUns
            if BOs==[]: break
    return AllBOs

def Obtain_Candidate_BOs(atom_list,Adj,Unsatu_set,DUs,Total_charge):
    """ It finds the 'best' candidate bond order (BO) matrix that fits the criteria, 
    from the given set of degrees of unsaturation of atoms. 
    Various BO matrices are generated with varying starting atoms and directions of bond order assignment,
    and the most appropriate one is selected.

    atom_list - list of atom class instances
    Adj - adjacency matrix
    Unsatu_set - list of indices of unsaturated atoms to which multiple bond orders can be assigned.
    DUs - list of degrees of unsaturation of atoms
    Total_charge - total charge of a molecule

    returns:
        Current_BO - the most appropriate BO matrix
    """
    Current_BO=np.copy(Adj)
    Deloc_BO,Current_Deloc_BondCount=setDelocalizedBonds(atom_list,Current_BO)

    for start_point in range(len(Unsatu_set)):
        BO=np.copy(Adj)
        DU=DUs[:]; unsatu_set=Unsatu_set[:]
        sumDU=sum(DU)
        unsatu_set.insert(0,unsatu_set.pop(start_point)) #change start point

        now=unsatu_set[0]
        while unsatu_set!=[]:
            num_neighbors,neighbors=Find_Neighbors(Adj,now,unsatu_set)

            if num_neighbors==0: unsatu_set.pop(unsatu_set.index(now)) #no neighbors to make multiple bonds with
            else:
                Num_neighbors_2nd=[]
                for i in range(num_neighbors):
                    num_neighbors_2nd,dummy=Find_Neighbors(Adj,neighbors[i],unsatu_set)
                    Num_neighbors_2nd.append(num_neighbors_2nd-1) #omitting 'now' vertex

                target=neighbors[Num_neighbors_2nd.index(min(Num_neighbors_2nd))] #which has minimum number of 2nd neighbor is picked as target
                Additional_bond=min(DU[now],DU[target])
                BO[now][target]+=Additional_bond
                BO[target][now]+=Additional_bond
                DU[now]-=Additional_bond; DU[target]-=Additional_bond

                if DU[now]<=0: unsatu_set.pop(unsatu_set.index(now))
                if DU[target]<=0: unsatu_set.pop(unsatu_set.index(target))

            if unsatu_set==[]: break

            remaining_atoms_neighbors=[]
            for remaining_atom in unsatu_set:
                num_neighbors,neighbors=Find_Neighbors(Adj,remaining_atom,unsatu_set)
                remaining_atoms_neighbors.append(num_neighbors)

            minimum_neighbor=remaining_atoms_neighbors.index(min(remaining_atoms_neighbors))
            now=unsatu_set[minimum_neighbor]

        Deloc_BO,New_Deloc_BondCount=setDelocalizedBonds(atom_list,BO)

        if np.sum(BO)-np.sum(Adj)==sumDU: return BO
        if New_Deloc_BondCount > Current_Deloc_BondCount: Current_BO=np.copy(BO); Current_Deloc_BondCount=New_Deloc_BondCount
        if np.sum(BO) > np.sum(Current_BO): Current_BO=np.copy(BO);                Current_Deloc_BondCount=New_Deloc_BondCount 

    return Current_BO


def AdjtoBO(atom_list,Adj,Total_charge,TotalChargeMethod,ObtainAllResonances=False):
    """ It converts adjacency matrix into bond order matrix.

    atom_list - list of atom class instances
    Adj - adjacency matrix
    Total_charge - total charge of a molecule
    TotalChargeMethod - 'SumofFragments' or 'Ionic'.
                        'SumofFragments': total charge is calculated as sum of fragment charges which are given in the input file
                        'Ionic': total charge is calculated by information on bond orders and valences
    ObtainAllResonances - True or False

    returns:
        a bond order matrix - if ObtainAllResonances==False
        list of possible bond order matrices - otherwise
    """
    Elemlist=[str.upper(x.get_element()) for x in atom_list]
    #Zlist=[am.getZ(x) for x in Elemlist]
    Zlist = [x.get_atomic_number() for x in atom_list]

    DUs=[]; list_count=0
    Terminal_oxygens=[]
    for i in range(len(atom_list)): Terminal_oxygens.append([])

    for i in range(len(Elemlist)):
        a=Elemlist[i]
        if a=='CO' or a=='RH' or a=='NI' or a=='TI' or a=='PD' or  a=='FE' or a=='HF':
            if int(np.sum(Adj[i]))>=6: du=0
            else: du=6-int(np.sum(Adj[i]))  #For make M-L multiple bonds possible
        elif a=='I':
            if int(np.sum(Adj[i]))>=6: du=0
            else: du=3-int(np.sum(Adj[i]))  #For make M-L multiple bonds possible
        elif a=='P' or a=='N':
            if a=='P' and int(np.sum(Adj[i]))==4: du=[0,1]; list_count+=1
            elif a=='P' and int(np.sum(Adj[i]))==3: du=[0,2]; list_count+=1
            elif a=='N' and int(np.sum(Adj[i]))==1 and (Elemlist[list(Adj[i]).index(1)]=='N' or Elemlist[list(Adj[i]).index(1)]=='C'): du=[1,2]; list_count+=1
            elif int(np.sum(Adj[i]))==3: 
                if list(map(lambda x:Elemlist[x],[a for a,b in enumerate(list(Adj[i])) if b==1])).count('N')>=1: du=0
                else: du=[0,1]; list_count+=1
            elif int(np.sum(Adj[i]))==2: du=[1,2]; list_count+=1
            elif int(np.sum(Adj[i]))==1: du=[2,3]; list_count+=1
            else: du=am.MaxV(atom_list[i].get_element())-int(np.sum(Adj[i]))
        elif a=='B':
            du=3-int(np.sum(Adj[i]))
        elif a=='S':
            if int(np.sum(Adj[i]))==1: du=1 
            elif int(np.sum(Adj[i]))==2: du=[0,2]; list_count+=1 
            elif int(np.sum(Adj[i]))==3: du=[1,3]; list_count+=1 
            elif int(np.sum(Adj[i]))==4: du=[0,2]; list_count+=1 
            else: du=am.MaxV(atom_list[i].get_element())-int(np.sum(Adj[i]))
        elif a=='O' and int(np.sum(Adj[i]))==1:
            bonding_atom=list(Adj[i]).index(1)
            if str.upper(atom_list[bonding_atom].get_element())=='C' and int(np.sum(Adj[bonding_atom]))==1: du=2 #carbon monoxide
            elif str.upper(atom_list[bonding_atom].get_element())=='N' and int(np.sum(Adj[bonding_atom]))==3: Terminal_oxygens[bonding_atom].append(i); du=0 #for N-O single bond 
            else: du=1
        else:
            max_v = am.MaxV(atom_list[i].get_element())
            if max_v is None:
                print (atom_list[i].get_element())
            du = max_v-int(np.sum(Adj[i]))
            if du<0: 
                du=0
        DUs.append(du)

    New_Adj=np.copy(Adj)

    for i in range(len(Terminal_oxygens)):
        NObonds=Terminal_oxygens[i]
        if len(NObonds)==2: #NO2 terminal
            oxygen1=NObonds[0]
            New_Adj[i][oxygen1]=New_Adj[oxygen1][i]=2
            list_count-=1; DUs[i]=0

    for i in range(len(Elemlist)):
        a=Elemlist[i]
        if type(DUs[i])==int and  DUs[i]>0:
            neighbor_DUs=[]
            for j in range(len(New_Adj)):
                if New_Adj[i][j]!=0:
                    if type(DUs[j])==list: neighbor_DUs.append(1)
                    else: neighbor_DUs.append(DUs[j])
            if sum(neighbor_DUs)==0: 
                if type(DUs[i])==list: list_count-=1
                DUs[i]=0

    multiple_valence_count=sum(isinstance(i, list) for i in DUs)
    if multiple_valence_count!=list_count: list_count=multiple_valence_count

    if list_count==0 and sum(DUs)<=0:
        if ObtainAllResonances: return [New_Adj]
        return New_Adj #No unsatu atoms

    #Create possibilities for valences
    Candidate_DUs=[DUs]
    if list_count>0:
        while True:
            new_candidate=[]
            for DUs in Candidate_DUs:
                for i in range(len(DUs)):
                    if type(DUs[i])==list:
                        newDU1=DUs[:i]+[DUs[i][0]]+DUs[i+1:]
                        newDU2=DUs[:i]+[DUs[i][1]]+DUs[i+1:]
                new_candidate.append(newDU1); new_candidate.append(newDU2)
            Candidate_DUs=new_candidate
            if len(Candidate_DUs)==2**list_count: break

    ##sort candidate DUs
    for i in range(len(Candidate_DUs)):
        for j in range(i+1,len(Candidate_DUs)):
            sumDUi=sum(Candidate_DUs[i]); sumDUj=sum(Candidate_DUs[j])
            if sumDUi<sumDUj:
                temp=Candidate_DUs[i][:]
                Candidate_DUs[i]=Candidate_DUs[j][:]
                Candidate_DUs[j]=temp[:]

    Candidate_Unsatu_set=[]
    for DUs in Candidate_DUs:
        Unsatu_set=[]
        for i in range(len(DUs)):
            if DUs[i]>0: Unsatu_set.append(i)
        Candidate_Unsatu_set.append(Unsatu_set)

    AllResonances=[]

    Current_BO=np.copy(New_Adj)
    Current_FC=getFC(atom_list,Current_BO,Total_charge,TotalChargeMethod)
    #print ('Candidate_DU',Candidate_DUs)
    for i in range(len(Candidate_DUs)):
        if len(Candidate_Unsatu_set[i])>0:
            if ObtainAllResonances:
                AllBOs=FindResonances(atom_list,New_Adj,Candidate_Unsatu_set[i],Candidate_DUs[i],Total_charge)
                AllResonances+=AllBOs
            else:
                Candidate_BO=Obtain_Candidate_BOs(atom_list,New_Adj,Candidate_Unsatu_set[i],Candidate_DUs[i],Total_charge)
                Candidate_FC=getFC(atom_list,Candidate_BO,Total_charge,TotalChargeMethod)
                #print ('FC',Candidate_FC)
                charge_balance=(Total_charge==sum(Candidate_FC))

                if charge_balance: 
                    NumberOfNeutralAtoms=Candidate_FC.count(0)
                    if len(Candidate_FC)-abs(int(Total_charge))==NumberOfNeutralAtoms: 
                        return Candidate_BO

                    if NumberOfNeutralAtoms > Current_FC.count(0) or np.sum(Candidate_BO) > np.sum(Current_BO):
                        Current_BO=np.copy(Candidate_BO)
                        Current_FC=Candidate_FC[:]
                elif np.sum(Candidate_BO) > np.sum(Current_BO): 
                    Current_BO=np.copy(Candidate_BO)
                    Current_FC=Candidate_FC[:]
    if not ObtainAllResonances:
        MetalElements=['SC','TI','V' ,'CR','MN','FE','CO','NI','CU','ZN',\
                        'Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD',\
                        'LU','HF','TA','W' ,'RE','OS','IR','PT','AU','HG'] 

        if len(set(Elemlist) & set(MetalElements))!=0: Current_BO=Detect_MetalCarbonyl(Elemlist,Current_BO)

    if ObtainAllResonances: 
        Final_Resonances=[]
        Ceiglist=[]
        for oneBO in AllResonances:
            Resonance_FC=[int(x) for x in getFC(atom_list,oneBO,Total_charge,TotalChargeMethod)]
            if sum(Resonance_FC)!=int(Total_charge): 
                continue
            if 2 in Resonance_FC or -2 in Resonance_FC: 
                continue

            samecheck=False
            for prev_BO in Final_Resonances:
                if np.array_equal(prev_BO,oneBO): 
                    samecheck=True
                    break
            if not samecheck: 
                Final_Resonances.append(oneBO)
        return Final_Resonances
    return Current_BO

def getResonanceBO(IM,Fragments_charge,TotalChargeMethod):
    """ It returns a bond order (BO) matrix containing 1.5,2.5-order bonds based on consideration of all possible resonance BO matrices

    IM - an Intermediate instance
    Fragments_charge - charge of fragments given in the input file
    TotalChargeMethod - 'SumofFragments' or 'Ionic'.
                        'SumofFragments': total charge is calculated as sum of fragment charges which are given in the input file
                        'Ionic': total charge is calculated by information on bond orders and valences

    returns:  a bond order (BO) matrix containing 1.5,2.5-order bonds 
    """
    reducedlist=IM.atom_indices_each_molecule  #atom indices of each molecule in an intermediate
    IM_BO=IM.Adj
    bondlist=BOmatrixtoList(IM_BO)
    IM_Adj=makeAdjacency(len(IM_BO),bondlist)
    All_atom_list=IM.atom_list
    
    NewBO_withResonance=np.copy(IM_BO)
    molecule_index=0
    for indices_of_a_molecule in reducedlist:
        A_Adj=np.array([[[l[a] for a in indices_of_a_molecule] for l in IM_Adj][b] for b in indices_of_a_molecule])
        A_BO=np.array([[[l[a] for a in indices_of_a_molecule] for l in IM_BO][b] for b in indices_of_a_molecule])
        A_atom_list=[All_atom_list[x] for x in indices_of_a_molecule]
        A_bondlist=BOmatrixtoList(A_Adj)

        if TotalChargeMethod=='SumofFragments':
            Fragnumlist=list(set([x.molecule_index for x in IM.molecule_list[molecule_index]]))
            Total_Charge=sum([Fragments_charge[x] for x in Fragnumlist])
        if TotalChargeMethod=='Ionic':
            Total_Charge=GetTotalCharge_Ionic(A_atom_list,A_BO)

        All_Resonances=AdjtoBO(A_atom_list,A_Adj,Total_Charge,TotalChargeMethod,True)

        for A_bond in A_bondlist:
            atom1=A_bond[0]; atom2=A_bond[1]

            BondOrders_EachResonance=sorted(list(set([int(x[atom1][atom2]) for x in All_Resonances])))
            Original_index_atom1=indices_of_a_molecule[atom1]
            Original_index_atom2=indices_of_a_molecule[atom2]

            if BondOrders_EachResonance==[1,2]: NewBO_withResonance[Original_index_atom1][Original_index_atom2]=NewBO_withResonance[Original_index_atom2][Original_index_atom1]=1.5
            if BondOrders_EachResonance==[2,3]: NewBO_withResonance[Original_index_atom1][Original_index_atom2]=NewBO_withResonance[Original_index_atom2][Original_index_atom1]=2.5
        molecule_index+=1

    return NewBO_withResonance


def getFC(atom_list,BO,chg,TotalChargeMethod='SumofFragments'):
    """ It calculates formal charges of atoms in a molecule.

    atom_list - list of atom class instances 
    BO - a bond order matrix
    chg - total charge of a molecule 
    TotalChargeMethod - 'SumofFragments' or 'Ionic'.
                        'SumofFragments': total charge is calculated as sum of fragment charges which are given in the input file
                        'Ionic': total charge is calculated by information on bond orders and valences

    returns: a list of formal charges of atoms in a molecule
    """
    TotalElectrons=sum([x.get_atomic_number() for x in atom_list])-chg
    
    FC=[]
    bo_matrix_original = np.copy(BO)
    for i in range(len(atom_list)):
        BondCount=0
        for j in range(len(atom_list)):
            BondCount+=BO[i][j]
        if len(atom_list)==1: fc=chg 
        elif str.upper(atom_list[i].get_element())=='H':
            if BondCount==1: fc=0
            else:
                if len(atom_list)==1: fc=chg
                else: fc=0
        elif str.upper(atom_list[i].get_element())=='NA': fc=1-BondCount
        elif (str.upper(atom_list[i].get_element())=='F' or str.upper(atom_list[i].get_element())=='CL') and len(atom_list)==1: fc=-1
        elif str.upper(atom_list[i].get_element())=='S' and BondCount>=4: 
            if BondCount==4 or BondCount==6: fc=0 #valence expansion
            if BondCount==5: fc=1 #valence expansion
        elif str.upper(atom_list[i].get_element())=='P' and BondCount>=4: fc=5-BondCount;  #fc=0 #valence expansion
        elif str.upper(atom_list[i].get_element())=='RH' or str.upper(atom_list[i].get_element())=='CO' or str.upper(atom_list[i].get_element())=='NI' or str.upper(atom_list[i].get_element())=='TI' or str.upper(atom_list[i].get_element())=='PD' or str.upper(atom_list[i].get_element())=='FE': fc=0 # Transition metals - tricky..
        elif str.upper(atom_list[i].get_element())=='C' and list(BO[i]).count(1)==2 and list(BO[i]).count(2)==0: fc=0 #carbene
        elif str.upper(atom_list[i].get_element())=='C' and list(BO[i]).count(1)==0 and list(BO[i]).count(2)==1: fc=0 #=C terminus
        elif str.upper(atom_list[i].get_element())=='B':
            if np.sum(BO[i])==0: fc=0
            else: fc=3-np.sum(BO[i])
        elif BondCount==0: fc=0 #For a single atom, regarded as zero
        elif str.upper(atom_list[i].get_element())=='I' and BondCount>1: #valence shell expansion for I
            if BondCount==3 or BondCount==5: fc=0
            if BondCount==2 or BondCount==4: fc=1
        else:
            period,group = atom_list[i].get_period_group()
            if group >= 4:
                lonepair=(4-BondCount)*2
            else:
                lonepair=(group-BondCount) * 2

            fc=group-BondCount-lonepair
        FC.append(int(fc))
        '''
        try: FC.append(int(fc))
        except UnboundLocalError:
            SMILES,atom_labels,NumRings=smiles.GetSMILES(atom_list,BO,     BOmatrixtoList(BO),      np.zeros(( len(atom_list) )),'N')
            print(SMILES)
        '''
    for i in range(len(atom_list)):
        if str.upper(atom_list[i].get_element())=='C' and list(BO[i]).count(1)==3: #Three single bonds. carbocation or carbanion
            if chg>0: FC[i]=1
            elif chg<0: FC[i]=-1
            else: FC[i]=0
    
    if TotalChargeMethod=='SumofFragments' and TotalElectrons%2!=0 and sum(FC)!=chg:
        for i in range(len(FC)):
            if FC[i]!=0:  FC[i]=0 #RADICAL
            if sum(FC)==chg: break

    #Exceptional
    if len(atom_list)==2 and sorted([str.upper(atom_list[0].get_element()),str.upper(atom_list[1].get_element())])==['C','O'] and int(BO[0][1])==3: FC=[0,0]

    return FC


def GetTotalCharge_Ionic(atom_list,BO):
    """ It calculates total charge of a molecule based on bond order and valence information.

    atom_list - list of atom class instances
    BO - bond order matrix of a molecule 
    bondlist - list of bonds in a molecule. Each bond is given as a tuple of two atom indices.

    returns: Total charge of a molecule.
    """
    MetalElements=['SC','TI','V' ,'CR','MN','FE','CO','NI','CU','ZN',\
            'Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD',\
               'LU','HF','TA','W' ,'RE','OS','IR','PT','AU','HG']

    Element_list=[str.upper(x.get_element()) for x in atom_list]
    Total_Charge=0
    bondlist = BOmatrixtoList(BO)
    if Element_list==['H','H']: 
        return 0

    for i in range(len(Element_list)):
        if Element_list[i]=='H': Total_Charge+=1
        elif Element_list[i]=='C': 
            if list(BO[i]).count(1)==2 and list(BO[i]).count(2)==0: #carbene
                Total_Charge-=2 
            else: Total_Charge-=4
        elif Element_list[i]=='N': Total_Charge-=3
        elif Element_list[i]=='O': Total_Charge-=2
        elif Element_list[i]=='F': Total_Charge-=1
        elif Element_list[i]=='CL': Total_Charge-=1
        elif Element_list[i]=='BR': Total_Charge-=1
        elif Element_list[i]=='I': Total_Charge-=1
        elif Element_list[i]=='S': 
            if sum(BO[i])<=3: Total_Charge-=2
            elif sum(BO[i])<=5: Total_Charge-=4
            elif sum(BO[i])<=6: Total_Charge-=6
        elif Element_list[i] in MetalElements: Total_Charge-=sum(BO[i])

    for bond in bondlist:
        if Element_list[bond[0]]!='H' and Element_list[bond[1]]!='H':
            BondOrder=BO[bond[0]][bond[1]]
            Total_Charge+=2*BondOrder
    return Total_Charge



def getunsatuatom_list(FragAdj,Elemlist,reactive_atoms,UserDefinedMaxV):
    """ It returns information on active atoms.

    FragAdj - list of adjacency matrices of input fragments
    Elemlist - List of atom elements of all fragments (e.g. ['C','O','C','O',...])
    reactive_atoms - 
        'auto' -> find active atoms automatically by calculating degrees of unsaturation of each atom
        lists of atom indices -> active atoms are given from the input file
                                (e.g. [[0,1],[0,2],...] - atoms 0,1 in fragment 0 and atoms 0,2 in fragment 1 are assigned as active ones)
    UserDefinedMaxV - list of maximum valences of elements defined by user
                      (e.g. ['N',4,'O',3])

    returns:
         unsatuatom_list - list of indices of active atoms 
         unsatuDSD - list of degrees of unsaturation of active atoms
         Fragnum - list of indices of fragments to which active atoms belong.
         Zlist - list of atomic numbers of active atoms
    """
    UserMaxV=dict()
    for i in range(0,len(UserDefinedMaxV),2):
        try:
            UserMaxV[str.upper(UserDefinedMaxV[i])]=UserDefinedMaxV[i+1]
        except TypeError:  
            UserMaxV[UserDefinedMaxV[i]-1]=UserDefinedMaxV[i+1]

    Elem_MaxV=[]

    for i in range(len(Elemlist)):
        if i in list(UserMaxV.keys()): Elem_MaxV.append(UserMaxV[i]); continue 

        try:
            Elem_MaxV.append(UserMaxV[str.upper(Elemlist[i])])
        except KeyError:
            Elem_MaxV.append(am.MaxV(Elemlist[i]))


    FragsZlist=[]
    for i in range(len(Elemlist)): 
        FragsZlist.append(am.getZ(Elemlist[i]))
        #FragsZlist.append()
    AtomCount=0; unsatuatom_list=[]; unsatuDSD=[]; Fragnum=[]; Zlist=[]; GroupClist=[]
    Zlist_atomicnumber=[]
    FragZ=1.5

    for i in range(len(FragAdj)):
        for j in range(len(FragAdj[i])):
            if reactive_atoms!='auto' and reactive_atoms[i].count(j)==0: continue #skip if it is not user-defined atom
            SNCount=0
            for k in range(len(FragAdj[i])):
                if FragAdj[i][j][k]==1: SNCount+=1

            if Elem_MaxV[AtomCount+j]>SNCount: 
                unsatuatom_list.append(AtomCount+j); unsatuDSD.append(Elem_MaxV[AtomCount+j]-SNCount); Fragnum.append(i)
                if len(FragAdj[i])==1: Zlist.append(am.getZ(Elemlist[AtomCount+j]))
                else:
                    temp=np.copy(FragAdj[i]) 
                    Groups_adj=np.delete(np.delete(temp,(j),axis=0),(j),axis=1)
                    Groups_Zlist=FragsZlist[AtomCount:AtomCount+len(FragAdj[i])]
                    Real_atomicnumber=Groups_Zlist[j]
                    Groups_Zlist[j:j+1]=[]

                    Groups_eigenval=getCoulombic(Groups_Zlist,Groups_adj)

                    Newgroupcheck=1
                    for ii in range(len(GroupClist)):
                        if isSameMolecule(GroupClist[ii],Groups_eigenval)==True and Real_atomicnumber==Zlist_atomicnumber[ii]: Newgroupcheck=0; Zlist.append(FragZ+ii)
                    if Newgroupcheck==1:
                        Zlist.append(FragZ+len(GroupClist))
                        GroupClist.append(Groups_eigenval)
                        Zlist_atomicnumber.append(Real_atomicnumber)
        AtomCount+=len(FragAdj[i])


    return unsatuatom_list, unsatuDSD, Fragnum, Zlist

def Find_Neighbors(Adj,vertex,unsatu_set):
    """ It finds neighboring unsaturated atoms of the given atom. 
    It is used in generating bond order matrix from adjacency matrix.

    Adj - adjacency matrix
    vertex - index of the current atom of which neighboring atoms will be investigated
    unsatu_set - set of atom indices of unsaturated atoms
    
    returns:
        num_neighbors - number of neighboring atoms
        neighbors - list of indices of neighboring atoms
    """
    num_neighbors=0
    neighbors=[]
    for i in range(len(unsatu_set)):
        vertex2=unsatu_set[i]
        if Adj[vertex][vertex2]!=0: 
            num_neighbors+=1
            neighbors.append(vertex2)

    return num_neighbors,neighbors




def permuteBO(BO,atom_labels):
    """ It permutes atom indices of the bond order matrix, so that it matches with those of 3D geometry generated from SMILES.

    BO - The (original) bond order matrix
    atom_labels - It contains information on correspondence of atom indices of the bond order matrix with those of SMILES string
                  (e.g. [3,'H','H',0,'H'] indicates that the atom 3 in BO matrix corresponds to the atom 0 in SMILES,
                  and atoms 1, 2, 4 are indices of hydrogen atoms in SMILES)
                  This is used for matching atom numberings between BO matrix and the 3D geometry generated by pybel.
                  (Please refer to def.permuteBO and def.UFFrelax) 
    returns:
        newBO - permuted bond order matrix
    """
    newBO=np.zeros((len(BO),len(BO)))
    for i in range(len(BO)):
        try:
            a=atom_labels.index(i)
        except ValueError: 
            continue #not in the atom_labels list because it is hydrogen.
        for j in range(i+1,len(BO)):
            try:
                b=atom_labels.index(j)
            except ValueError: 
                continue #not in the atom_labels list because it is hydrogen.
            if BO[i][j]!=0: 
                newBO[a][b]=newBO[b][a]=BO[i][j]

    bonding_nonhydrogen_index=None
    for i in range(len(atom_labels)):
        if atom_labels[i]!='H': 
            bonding_nonhydrogen_index=i
        else: 
            newBO[bonding_nonhydrogen_index][i]=newBO[i][bonding_nonhydrogen_index]=1

    '''
    hydrogen_index=len(atom_labels)
    for i in range(len(numH)):
        if type(numH[i])==str:
            Num_Explicit_H=int(numH[i][1:])
            for H in range(Num_Explicit_H): newBO[i+H+1][i]=newBO[i][i+H+1]=1
            continue
        for j in range(numH[i]):
            newBO[hydrogen_index][i]=newBO[i][hydrogen_index]=1
            hydrogen_index+=1
    '''
    return newBO


