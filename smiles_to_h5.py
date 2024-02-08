### smiles to H5 fingerprint conversion ###
""" Use FPSim2 to convert smiles to H5 fingerprints.
    Input: csv file with smiles
    Output: H5 file with fingerprints"""

### Imports ###
import argparse 
from FPSim2.io import create_db_file
import numpy as np
import pandas as pd
import pickle

### Parsing arguments ###
args = argparse.ArgumentParser(description='Convert SMILES to H5 fingerprints')
args.add_argument('-i', '--input', help='Input file with smiles', required=True)
args.add_argument('-c', '--column', help='Column with smiles', required=True)
args.add_argument('-fp', '--fingerprint', help='Fingerprint type', required=True)
args.add_argument('-r', '--radius', help='Radius for Morgan fingerprints', required=False, default=2, type=int)
args.add_argument('-b', '--bits', help='Number of bits for Morgan fingerprints', required=False, default=2048, type=int)
args.add_argument('-o', '--output', help='Output H5 file', required=True)
args = args.parse_args()

### Read input file ###
if args.input.endswith('.csv'):
    data = pd.read_csv(args.input)
    smiles = data[args.column].tolist()
elif args.input.endswith('.pkl'):
    file = open(args.input, 'rb')
    data = pickle.load(file)
    file.close()
    smiles = data[args.column]

### Convert smiles to H5 fingerprints ###
if args.fingerprint != 'Morgan':
    create_db_file(smiles, args.output, args.fingerprint, gen_ids=True)
else:
    create_db_file(smiles, args.output, args.fingerprint, {'radius': args.radius, 'nBits': args.bits}, gen_ids=True)




