import os
import numpy as np
import prody
from scipy.spatial.distance import pdist,squareform

class Ligand():
    def __init__(self,path,pdbqt_dir_path):
        self.name = path.split("/")[-1].split(".")[0]
        self.pdbqt_path = "%s/%s.pdbqt"%(pdbqt_dir_path,self.name)
        self.pdb_path = "%s/%s.pdb"%(pdbqt_dir_path,self.name)
        self.make_pdbqt(path,self.pdbqt_path)
        self.read_pdbqt()
        self.make_pdb(self.pdbqt_path,self.pdb_path)
        self.n = len(self.coords)
        self.d_matrix =squareform(pdist(self.coords))
        self.ligand = prody.parsePDB(self.pdb_path)
        
    def make_pdbqt(self,in_path,out_path):
        os.system("prepare_ligand -l %s -o %s"%(in_path,out_path))
        return os.path.exists(out_path)
    
    def make_pdb(self,in_path,out_path):
        os.system("obabel %s -O %s"%(in_path,out_path))
        
    def read_pdbqt(self):
        coords = []
        qs = []
        autodock_atom_types = []
        f = open(self.pdbqt_path)
        for line in f.readlines():
            if line[0:6] in ["ATOM  ","HETATM"]:
                coords.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
                qs.append(float(line[66:76]))
                autodock_atom_types.append(line[77:79].strip())
        f.close()
        self.coords = np.array(coords)
        self.qs = np.array(qs)
        self.autodock_atom_types = np.array(autodock_atom_types) 
        
    
        
    
        
