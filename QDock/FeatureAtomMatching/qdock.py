# -*- coding: utf-8 -*-
"""
@author: Jinyin Zha
"""
#Basic Packages
import os
import sys
import copy
import numpy as np
import prody
#QC simulator
import neal 
from pyqubo import Binary, Constraint
#Receprot and Ligand
sys.path.append(os.path.abspath(os.path.realpath(__file__))[:-8])
from famreceptor import Receptor
from famligand import Ligand

class FAMDock():
    def __init__(self):
        self.__step = "receptor"
        self.w_dict = {'H': 2.2,
         'C': 2.55,
         'N': 3.04,
         'O': 3.44,
         'F': 3.98,
         'Si': 1.9,
         'P': 2.19,
         'S': 2.58,
         'Cl': 3.16,
         'As': 2.18,
         'Se': 2.48,
         'Br': 2.96,
         'I': 2.66,
         'B': 2.04} 
        
    def get_step(self):
        return self.__step
    
    #Step 1, Prepare receptor file. 
    def make_receptor(self,receptor_path):
        if receptor_path[-4:] != ".pdb":
            print("Receptor shoud be PDB format!")
            return
        self.receptor = Receptor(receptor_path,"receptor.pdbqt")
        self.receptor_autodock_atom_types = np.unique(
            self.receptor.autodock_atom_types)
        self.__step = "ligand"
        
    #Step 2, Prepare Ligand files to be docked.     
    def make_ligand(self,ligand_paths):
        if self.__step != "ligand":
            print("Please run make_%s first!"%(self.__step))
            return
        self.ligands = []
        if not os.path.exists("Ligands"):
            os.mkdir("Ligands")
        for path in ligand_paths:
            self.ligands.append(Ligand(path,"Ligands"))
        self.ligands_autodock_atom_types = np.unique(np.hstack(
            [i.autodock_atom_types for i in self.ligands]))
        self.__step = "box"
    
    #Step 3, Create Docking Box with a ligand
    def make_box_ligand(self,path,
                           center_length=8,grid_length=1.0):
        if self.__step != "box":
            print("Please run make_%s first!"%(self.__step))
            return
        if not os.path.exists("Box_Rawligand"):
            os.mkdir("Box_Rawligand")
        self.raw_ligand = Ligand(path,"Box_Rawligand")
        self.box_center = np.mean(self.raw_ligand.coords,axis=0)
        self.box_lengths = center_length + 2 * np.max(np.abs(
            self.raw_ligand.coords - self.box_center),axis=0)
        self.grid_length = grid_length
        self.__autosite()
        
    #Step 3, Create Docking Box with manual input of center coordiante and lengths
    def make_box_input(self,x,y,z,dx,dy,dz,grid_length=1.0):
        self.box_center = np.array([x,y,z])
        self.box_lengths = np.array([dx,dy,dz])
        self.grid_length = grid_length
        self.__autosite()
    
    #Step 4, Docking (serial)
    def dock(self,edge_cutoff,K_dist,K_mono,n_pos=30,
             save_qubo=True,sim_dock=True,save_match=True,save_pose=True):
        if self.__step != "dock":
            print("Please run make_%s first!"%(self.__step))
            return
        return [self.indiv_dock(ligand,edge_cutoff,K_dist,K_mono,n_pos,
                                save_qubo,sim_dock,save_match,save_pose) \
                for ligand in self.ligands]
            
    #Step 4, Docking(single)
    def indiv_dock(self,ligand,edge_cutoff,K_dist,K_mono,n_pos=30,
                   save_qubo=True,sim_dock=True,save_match=True,save_pose=True):
        if self.__step != "dock":
            print("Please run make_%s first!"%(self.__step))
            return
        #4.1 Create QUBO
        H = 0
        Vs = []
        #4.1.1 Make Vertexes (quadratic term in QUBO)
        for i in range(ligand.n):
            t_lig = ligand.autodock_atom_types[i]
            if t_lig not in self.w_dict.keys():
                t_lig = t_lig[0]
            if t_lig == "A":
                t_lig = "C"
            for j,fs in enumerate(self.feature_atoms):
                t_fs = fs.getElement()
                weight = abs(self.w_dict[t_lig] - self.w_dict[t_fs]) - 0.5
                x = Binary("%d_%d_%5.5f"%(i,j,weight))
                Vs.append([x,i,j,weight])
                H += weight * x**2
        if len(Vs)*(len(Vs)-1)/2 > 100000000:
            print("Too Big!!!")
            return
        #4.1.2, Make Edges (cross terms in QUBO)
        C_dist = 0
        C_mono = 0
        for p in range(len(Vs)):
            x1,i1,j1,w1 = Vs[p]
            for q in range(p+1,len(Vs)):
                x2,i2,j2,w2 = Vs[q]
                dd = abs(ligand.d_matrix[i1,i2] - self.box_Dmatrix[j1,j2])
                if dd > edge_cutoff:
                    C_dist += x1 * x2
                if  i1==i2:
                    C_mono += x1 * x2 
        H += K_dist * Constraint(C_dist, label='link') 
        H += K_mono * Constraint(C_mono, label='mono')
        model = H.compile()       
        qubo, offset = model.to_qubo()
        #4.1.3, Save QUBO Models
        if not os.path.exists("QUBOs"):
            os.mkdir("QUBOs")
        if save_qubo:
            np.save("QUBOs/%s.npy"%(ligand.name),qubo)
        #4.2 Docking by PyQUBO simulator (simulated annealing)
        if not sim_dock:
            return 
        sampler = neal.SimulatedAnnealingSampler()
        raw_solution = sampler.sample_qubo(qubo,num_reads=n_pos,seed=42)
        samples = []
        for sample in raw_solution.samples():
            #if sample not in samples:
                samples.append(sample)
        #4.3, Convert matches back to docking pose
        news = []
        match = []
        for sample in samples:
            ls = []
            gs = []
            ws = []
            this_match = []
            for info in filter(lambda x:1== x[1],sample.items()):
                i,j,w = info[0].split("_")
                i = int(i)
                j = int(j)
                w = float(w)
                ls.append(ligand.coords[i])
                gs.append(self.feature_atoms.getCoords()[j])
                ws.append(w)
                this_match.append(info[0])
            try:
               ls = np.vstack(ls)
               gs = np.vstack(gs)
               match.append(copy.deepcopy(this_match))
            except:
               continue
            trans=prody.superpose(ls,gs)[1]
            news.append(prody.applyTransformation(trans,ligand.coords))
        if not os.path.exists("Matches"):
            os.mkdir("Matches")
        if save_match and match:
            np.save("Matches/%s_match.npy"%(ligand.name),match) 
        if not os.path.exists("Poses"):
                os.mkdir("Poses")
        if save_pose and news:
            tmp = ligand.ligand.copy()
            tmp.setCoords(news[0])
            if len(news) > 1:
                tmp.addCoordset(np.array(news[1:]))
            prody.writePDB("Poses/%s_poses.pdb"%(ligand.name),tmp)
        return np.array(news)
    
    
    def __autosite(self):
        self.dims = (self.box_lengths / self.grid_length).astype(np.int)
        for i in range(3):
            if self.dims[i] % 2 == 0:
                self.dims[i] = self.dims[i] - 1
        os.system("mkdir pocs")
        os.system("autosite -r %s --boxcenter [%1.3f,%1.3f,%1.3f]\
                  --boxdim [%d,%d,%d] -o pocs"%(
            self.receptor.pdbqt_path,
            self.box_center[0],self.box_center[1],self.box_center[2],
            self.dims[0]-1,self.dims[1]-1,self.dims[2]-1))        
        pocs = []
        poc_files = os.listdir("pocs")
        for poc_file in poc_files:
            if "_fp_" not in poc_file:
                continue
            else:
                pocs.append(prody.parsePDB("pocs/%s"%(poc_file)))
        self.feature_atoms = pocs[0]
        if len(pocs) > 1:
            for i in range(1,len(pocs)):
                self.feature_atoms += pocs[i] 
        self.box_Dmatrix = prody.buildDistMatrix(self.feature_atoms)
        self.__step = "dock"

