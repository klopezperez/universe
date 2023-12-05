### smiles to H5 fingerprint conversion ###
""" Use FPSim2 to convert smiles to H5 fingerprints.
    Input: csv file with smiles
    Output: H5 file with fingerprints"""

### Imports ###
import numpy as np
from FPSim2.io import create_db_file
import argparse 
import pandas as pd

### Parsing arguments ###
args = argparse.ArgumentParser(description='Convert SMILES to H5 fingerprints')
args.add_argument('-i', '--input', help='Input csv file with smiles', required=True)
args.add_argument('-c', '--column', help='Column with smiles', required=True)
args.add_argument('-fp', '--fingerprint', help='Fingerprint type', required=True)
args.add_argument('-r', '--radius', help='Radius for Morgan fingerprints', required=False, default=2, type=int)
args.add_argument('-b', '--bits', help='Number of bits for Morgan fingerprints', required=False, default=2048, type=int)
args.add_argument('-o', '--output', help='Output H5 file', required=True)
args = args.parse_args()

### Read csv input file ###
df = pd.read_csv(args.input, sep=',')
smiles = df[args.column].tolist()
#ids = df['id'].tolist()
#fpsim2list = [[ids[i], smiles[i]] for i in range(len(smiles))]

### Convert smiles to H5 fingerprints ###
if args.fingerprint != 'Morgan':
    create_db_file(smiles, args.output, args.fingerprint, gen_ids=True)
else:
    create_db_file(smiles, args.output, args.fingerprint, {'radius': args.radius, 'nBits': args.bits}, gen_ids=True)




