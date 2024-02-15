from galaxy_class import *
import glob
import pandas as pd

# Get the isim for each of the databases
isim_universe = []
for file in glob.glob('*.h5'):
    G = ChemGalaxy(file[:-3], file)
    isim_universe.append([G.name, G.size, G.isim, G.get_outliers(), G.get_medoids()])

# Create a dataframe with the isim of each database
isim_df = pd.DataFrame(isim_universe, columns = ['database', 'size', 'isim', 'outliers', 'medoids'])

# Sort alphabetically by database name
isim_df = isim_df.sort_values(by = 'database')

# Save the dataframe to a csv file
#isim_df.to_csv('isim_universe.csv', index = False)

print(isim_df)