import numpy as np

class Grid():
    def __init__(self,path,cutoff):
        f = open(path)
        es = np.array([float(i) for i in f.readlines()[6:]])
        if type(cutoff) == type(0.0):
            self.pos = np.argwhere(es<cutoff).flatten()
            self.es = es[self.pos]
        else:
            self.es = es
        self.n = len(self.es)
        
       
        
