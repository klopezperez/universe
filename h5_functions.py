from isim_comp import calculate_isim
from FPSim2 import FPSim2Engine
import numpy as np
import pickle

def h5_to_fp(h5file, fp='ECFP4'):
    """ Convert H5 fingerprints to numpy array """
    fpe = FPSim2Engine(h5file)
    fps = fpe.fps[:,1:-1].view('uint8')
    bits = np.unpackbits(fps[:, np.newaxis], axis=1).ravel()
    fps = bits.reshape(int(bits.size / (fps.shape[1]*8)), fps.shape[1]*8)
    if fp == 'MACCS':
        fps = fps[:, 1:167]
    return fps

