import os
import numpy as np

class Receptor():
    def __init__(self,path,pdbqt_path):
        self.name = path.split("/")[-1].split(".")[0]
        self.pdbqt_path = pdbqt_path
        self.make_pdbqt(path,self.pdbqt_path)
        self.read_pdbqt()
        
    def make_pdbqt(self,in_path,out_path):
        os.system("prepare_receptor -r %s -o %s -A hydrogens"%(in_path,out_path))
        return os.path.exists(out_path)
    
    def read_pdbqt(self):
        autodock_atom_types = []
        f = open(self.pdbqt_path)
        for line in f.readlines():
            if line[0:6] in ["ATOM  ","HETATM"]:
                autodock_atom_types.append(line[77:79].strip())
        f.close()
        self.autodock_atom_types = np.array(autodock_atom_types) 
        
