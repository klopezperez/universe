from galaxy_class import *
import glob
import pandas as pd
import pickle

# Get the isim for each of the databases
isim_universe = []
for file in glob.glob('*.h5'):
    G = ChemGalaxy(file[:-3], file)
    isim_universe.append([G.name, G.size, G.isim, G.get_outliers(), G.get_medoids(), G.get_ctotal(G.get_outliers()), G.get_ctotal(G.get_medoids())])

# Create a dataframe with the isim of each database
isim_df = pd.DataFrame(isim_universe, columns = ['database', 'size', 'isim', 'outliers', 'medoids', 'c_outliers', 'c_medoids'])

# Sort alphabetically by database name
isim_df = isim_df.sort_values(by = 'database')
isim_df = isim_df.reset_index(drop = True)

# Safe the dataframe
#isim_df.to_csv('isim_universe.csv', index = False)

# Create a matrix with the intersection of outliers
# Read files containing the smiles 

outliers, medoids = isim_df['outliers'].tolist(), isim_df['medoids'].tolist()
outliers, medoids = [set(x) for x in outliers], [set(x) for x in medoids]
outlier_matrix, medoids_matrix = [], []
for i in range(len(outliers)):
    row_outlier = []
    row_medoid = []
    for j in range(len(outliers)):
        row_outlier.append(len(outliers[i].intersection(outliers[j]))/len(outliers[i].union(outliers[j])))
        row_medoid.append(len(medoids[i].intersection(medoids[j]))/len(medoids[i].union(medoids[j])))
    outlier_matrix.append(row_outlier)
    medoids_matrix.append(row_medoid)

# Create an isim matrix of outliers unions
outlier_union = []
for i in range(len(outliers)):
    row = []
    for j in range(len(outliers)):
        c_total = np.array(isim_df['c_outliers'][i] + isim_df['c_outliers'][j])
        n = len(outliers[i].union(outliers[j]))
        row.append(calculate_isim(c_total, n_objects = n, n_ary = 'JT'))
    outlier_union.append(row)

# Create an isim matrix of medoids unions
medoid_union = []
for i in range(len(medoids)):
    row = []
    for j in range(len(medoids)):
        c_total = np.array(isim_df['c_medoids'][i] + isim_df['c_medoids'][j])
        n = len(medoids[i].union(medoids[j]))
        row.append(calculate_isim(c_total, n_objects = n, n_ary = 'JT'))
    medoid_union.append(row)

# Save matrices as npy files
#np.save('outlier_matrix.npy', np.array(outlier_matrix))
#np.save('medoids_matrix.npy', np.array(medoids_matrix))
#np.save('outlier_union.npy', np.array(outlier_union))
#np.save('medoid_union.npy', np.array(medoid_union))
