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
        self.indexes = np.array(self.fpe.fps[:, 0]) - 1
        self.isim = calculate_isim(self.fingerprints, n_objects=self.size, n_ary=n_ary)
        self.comp_isim = None

    def __str__(self):
        return self.name
    
    def __repr__(self):
        return self.name

    def calc_comp_isim(self):
        self.comp_isim = calculate_comp_sim(self.fingerprints, n_ary=self.n_ary)   
        
    def get_outliers(self, percentage = 0.05):
        # Get the indexes of the outliers (mols with highest comp_isim values)
        if self.comp_isim is None: self.calc_comp_isim()
        num_outliers = int(self.size * percentage)
        comp_sim = [(index, value) for index, value in zip(self.indexes, self.comp_isim)]
        comp_sim.sort(key = lambda x: x[1])
        self.outliers = [index for index, value in comp_sim[-num_outliers:]]
        return self.outliers
    
    def get_medoids(self, percentage = 0.05):
        # Get the indexes of the medoids (mols with lowest comp_isim values)
        if self.comp_isim is None: self.calc_comp_isim()
        num_medoids = int(self.size * percentage)
        comp_sim = [(index, value) for index, value in zip(self.indexes, self.comp_isim)]
        comp_sim.sort(key = lambda x: x[1])
        self.medoids = [index for index, value in comp_sim[:num_medoids]]
        return self.medoids
    
    def get_comp_isims(self, indexes_list = None):
        if self.comp_isim is None: self.calc_comp_isim()
        if indexes_list is None: indexes_list = self.indexes
        indexes_list = np.array(indexes_list)
        indexes_list = np.where(np.isin(self.indexes, indexes_list))[0]
        print(indexes_list)
        return self.comp_isim[indexes_list]
    
    def get_fingerprints(self, indexes_list = None):
        if indexes_list is None: indexes_list = self.indexes
        indexes_list = np.array(indexes_list)
        indexes_list = np.where(np.isin(self.indexes, indexes_list))[0]
        return self.fingerprints[indexes_list]
    
    def get_ctotal(self, indexes_list = None):
        if indexes_list is None: indexes_list = self.indexes
        indexes_list = np.array(indexes_list)
        indexes_list = np.where(np.isin(self.indexes, indexes_list))[0]
        return np.sum(self.fingerprints[indexes_list], axis = 0)
    
    def get_isim(self, indexes_list):
        if indexes_list is None: return self.isim
        if type(indexes_list) == int: raise ValueError('Indexes_list must be a list of indexes, isim is not calculated for a single molecule.')
        indexes_list = np.array(indexes_list)
        indexes_list = np.where(np.isin(self.indexes, indexes_list))[0]
        return calculate_isim(self.fingerprints[indexes_list], n_objects=len(indexes_list), n_ary=self.n_ary)





