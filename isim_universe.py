from galaxy_class import *
import glob
import pandas as pd
import pickle

# Get the isim for each of the databases
isim_universe = []
for file in glob.glob('*.h5'):
    G = ChemGalaxy(file.split('/')[-1][:-3], file)
    isim_universe.append([G.name, G.size, G.isim, G.get_outliers(), G.get_medoids(), G.get_isim(G.get_outliers()), G.get_isim(G.get_medoids()), G.get_ctotal(G.get_outliers()), G.get_ctotal(G.get_medoids())])

# Create a dataframe with the isim of each database
isim_df = pd.DataFrame(isim_universe, columns = ['database', 'size', 'isim', 'outliers', 'medoids', 'outliers_isim', 'medoids_isim', 'c_outliers', 'c_medoids'])

# Sort alphabetically by database name
isim_df = isim_df.sort_values(by = 'database')
isim_df = isim_df.reset_index(drop = True)

# Safe the dataframe
isim_df.to_csv('isim_universe.csv', index = False)