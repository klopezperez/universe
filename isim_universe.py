from utils import *
from isim_comp import *
from FPSim2 import FPSim2Engine

# Script to calculate the isim of databases of molecules and get information about their similarity across time 

class ChemGalaxy:
    def __init__(self, name, fpe_data, n_ary='JT'):
        self.name = name
        self.fpe = FPSim2Engine(fpe_data)
        self.n_ary = n_ary
        self.fingerprints = unpack_fpe(self.fpe)
        self.size = len(self.fingerprints)
        self.indexes = self.fpe.fps[:, 0]
        self.isim = calculate_isim(self.fingerprints, n_objects=self.size, n_ary=n_ary)


    def __str__(self):
        return self.name
    
    def __repr__(self):
        return self.name

    def calc_comp_isim(self):
        self.comp_isim = calculate_comp_sim(self.fingerprints, n_ary=self.n_ary)   
        
    def get_outliers(self, percentage = 0.05):
        # Get the indexes of the outliers (mols with highest comp_isim values)
        if not self.comp_isim: self.calc_comp_isim()
        num_outliers = int(self.size * percentage)
        threshold = np.partition(self.comp_isim, -num_outliers)[-num_outliers]
        self.outliers = self.indexes[self.comp_isim > threshold]
        return self.outliers
    
    def get_medoids(self, percentage = 0.05):
        # Get the indexes of the medoids (mols with lowest comp_isim values)
        if not self.comp_isim: self.calc_comp_isim()
        num_medoids = int(self.size * percentage)
        threshold = np.partition(self.comp_isim, num_medoids)[num_medoids]
        self.medoids = self.indexes[self.comp_isim < threshold]
        return self.medoids
    
    def get_fingerprints(self, indexes_list = None):
        if indexes_list is None: indexes_list = self.indexes
        return self.fingerprints[indexes_list]
    
    def get_ctotal(self, indexes_list = None):
        if indexes_list is None: indexes_list = self.indexes
        return np.sum(self.fingerprints[indexes_list], axis = 0)
    
    def get_isim(self, indexes_list):
        if indexes_list is None: return self.isim
        else: return calculate_isim(self.fingerprints[indexes_list], n_objects=len(indexes_list), n_ary=self.n_ary)

import glob
for file in glob.glob('*.h5'):
    G = ChemGalaxy(file[:-3], file)
    print(G.name, G.size, G.isim)





