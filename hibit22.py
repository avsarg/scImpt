# -*- coding: utf-8 -*-
"""
scRNA-seq imputation for HIBIT22
"""
import os
os.chdir("E:\\HIBIT22")
os.getcwd()

from scripts import func_file as ff


################################################################################################################
######################################## PDAC-A dataset ########################################
################################################################################################################
seqA = ff.dset_pdacA()

seqA_magic = ff.load_impDset('seqA', 'magic', seqA)
seqA_saver = ff.load_impDset('seqA', 'saver', seqA)
seqA_alra = ff.load_impDset('seqA', 'alra', seqA)

pdacA = [seqA, seqA_alra, seqA_magic, seqA_saver]
del seqA, seqA_magic, seqA_saver, seqA_alra
labels = ['noimp', 'alra', 'magic', 'saver']

for i in range(len(pdacA)):
    # markers_plot(data, dName, suff_title)
    ff.markers_plot(pdacA[i], 'seqA', labels[i])


for i in range(len(pdacA)):
    # kmeans_clust(data, nclust, dName, suff_title)
    ff.kmeans_clust(pdacA[i], len(pdacA[0].obs['cell_labels'].cat.categories), 'seqA', labels[i])

# marker_in_deg(data, dName, labels)
ff.marker_in_deg(pdacA, 'seqA', labels)


################################################################################################################
######################################## PDAC-B dataset ########################################
################################################################################################################
seqB = ff.dset_pdacB()

seqB_magic = ff.load_impDset('seqB', 'magic', seqB)
seqB_saver = ff.load_impDset('seqB', 'saver', seqB)
seqB_alra = ff.load_impDset('seqB', 'alra', seqB)

pdacB = [seqB, seqB_alra, seqB_magic, seqB_saver]
del seqB, seqB_magic, seqB_saver, seqB_alra
labels = ['noimp', 'alra', 'magic', 'saver']

for i in range(len(pdacB)):
    # markers_plot(data, dName, suff_title)
    ff.markers_plot(pdacB[i], 'seqB', labels[i])


for i in range(len(pdacB)):
    # kmeans_clust(data, nclust, dName, suff_title)
    ff.kmeans_clust(pdacB[i], len(pdacB[0].obs['cell_labels'].cat.categories), 'seqB', labels[i])

# marker_in_deg(data, dName, labels)
ff.marker_in_deg(pdacB, 'seqB', labels)


################################################################################################################
######################################## pancreas dataset ########################################
################################################################################################################
panc = ff.dset_panc()

panc_magic = ff.load_impDset('panc', 'magic', panc)
panc_saver = ff.load_impDset('panc', 'saver', panc)
panc_alra = ff.load_impDset('panc', 'alra', panc)

pancs = [panc, panc_alra, panc_magic, panc_saver]
del panc, panc_magic, panc_saver, panc_alra
labels = ['noimp', 'alra', 'magic', 'saver']

for i in range(len(pancs)):
    # markers_plot(data, dName, suff_title)
    ff.markers_plot(pancs[i], 'panc', labels[i])


for i in range(len(pancs)):
    # kmeans_clust(data, nclust, dName, suff_title)
    ff.kmeans_clust(pancs[i], len(pancs[0].obs['cell_labels'].cat.categories), 'panc', labels[i])


# marker_in_deg(data, dName, labels)
ff.marker_in_deg(pancs, 'panc', labels)





