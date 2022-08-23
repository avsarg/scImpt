# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 12:37:08 2022

@author: Gulben AVSAR
"""
import magic
import scanpy as sc
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

def dset_pdacA():
    adata = sc.read_h5ad('./datasets/seqA_preprocessed.h5ad')
    del adata.uns['log1p']
    return(adata)

def dset_pdacB():
    adata = sc.read_h5ad('./datasets/seqB_preprocessed.h5ad')
    del adata.uns['log1p']
    return(adata)

def dset_panc():
    adata = sc.read_h5ad('./datasets/GSM2230757_human1_umifm_counts.h5ad')
    del adata.uns['log1p']
    return(adata)




def apply_magic(data, dName):
    X = data.copy()
    magic_operator = magic.MAGIC()
    all_magic = magic_operator.fit_transform(X, genes='all_genes')
    all_magic.to_df().to_csv('./datasets/{}_magic.txt'.format(dName), sep=' ', index=False)
    return()

def load_impDset(dName, suff_title, mainData):
    data = pd.read_csv('./datasets/'+dName+'_'+suff_title+'.txt', sep=' ')
    # data = pd.read_csv('./datasets/'+'seqA'+'_'+'alra'+'.txt', sep=' ')
    if suff_title == 'saver':
        data = data.T
    genes = data.columns.to_list()
    genes = pd.DataFrame(genes, index = genes)
    genes.rename(columns={genes.columns[0]:'Genes'}, inplace=True)
    cells = mainData.obs['cell_labels']
    cells = cells.rename('cell_labels')
    # cells = pd.DataFrame(cells, index =list(range(len(cells))))
    # cells['Cells'] = cells.index
    # cells.rename(columns={cells.columns[0]:'cell_labels'}, inplace=True)
    
    data.index = cells.index
    obs_idx = cells.index.to_frame()
    adata = sc.AnnData(data, var=genes, obs=obs_idx)
    adata.obs['cell_labels'] = cells
    del adata.obs[adata.obs.columns[0]]
    return(adata)


def markers_plot(data, dName, suff_title):
    if dName == 'seqA':
        markers ={'Acinar':['PRSS1','CTRB2','REG1A'],
                  'Alpha & Myeloid':['MRC1','MAFB','APOE','HLA-DRA','C1QA','CD14'],
                  'Beta':['NKX6-1','HEPACAM2','DKK3'],
                  'Cancer':['CEACAM5','MSLN','KRT17','LAMC2','KRT16'],
                  'Ductal':['SLC4A4','FXYD2','SPP1','TFF1', 'TFF2','TFF3'],
                  'Endothelial':['PLVAP','VWF'],
                  'Mast':['TPSAB1','CPA3'],
                  'RBCs':['HBB','HBA2'],
                  'T & NK cells':['CD3D','CD2','NKG7'],
                  'pDCs':['IRF7','GZMB'],
                  }
    elif dName == 'seqB':
            markers ={'Acinar':['PRSS1','CTRB2','REG1A'],
                      'Alpha & Myeloid':['MRC1','MAFB','APOE','HLA-DRA','C1QA','CD14'],
                      'Beta':['NKX6-1','HEPACAM2','DKK3'],
                      'Cancer':['CEACAM5','MSLN','KRT17','LAMC2','KRT16'],
                      'Ductal':['SLC4A4','FXYD2','SPP1','TFF1', 'TFF2','TFF3'],
                      'Endothelial':['PLVAP','VWF'],
                      'Fibroblast':['COL1A1','COL3A1','LUM','DCN'],
                      'Mast':['TPSAB1','CPA3'],
                      }
    else:
            markers = {'acinar':'PRSS1',
                       'as': ['PDGFRA', 'SPARC'],
                       'alpha':'GCG',
                       'beta':'INS',
                       'delta':'SST',
                       'duct':['SPP1', 'KRT19'],
                       'endo': ['FLT1', 'VWF'],
                       'epsilon':['GHRL', 'ANXA13'],
                       'gamma':'PPY',
                       'macrophage': 'SDS',
                       'mast': 'TPSAB1',
                       'qs': 'RGS5',
                       'schwann': 'SOX10',
                       't_cell':'TRAC',
                       }
    sc.set_figure_params(scanpy=True, fontsize=20)
    sc.pl.dotplot(data, markers,
                      groupby='cell_labels', cmap='BuPu',
                      # standard_scale='var',
                      save='{}_markers_{}.png'.format(dName, suff_title),
                      )
    return()




def kmeans_clust(data, nclust, dName, suff_title):
    sc.set_figure_params(scanpy=True, fontsize=20)
    if dName == 'panc' and suff_title == 'noimp':
        sc.tl.pca(data, svd_solver='arpack')
    if suff_title != 'noimp':
        sc.tl.pca(data, svd_solver='arpack')
    else:
        d = data.X.copy()
        s = silhouette_score(d, data.obs['cell_labels'])
        sc.pl.pca(data, color='cell_labels',
                  # frameon = False,
                  annotate_var_explained=True,
                  title='SI: {:0.2f}'.format(s),
                  save='_{}_Cells_{}.png'.format(dName, suff_title))
    # kmeans clustering
    d = data.X.copy()
    kmeans = KMeans(n_clusters=nclust, random_state=0).fit(d)
    data.obs['kmeans{}'.format(nclust)] = kmeans.labels_.astype(str)
    s = silhouette_score(d, data.obs['kmeans{}'.format(nclust)])
    # data
    sc.pl.pca(data, color='kmeans{}'.format(nclust),
              # frameon = False,
              annotate_var_explained=True,
              title='SI: {:0.2f}'.format(s),
              save='_{}_kClusters_{}.png'.format(dName, suff_title),
              )
    
    return()


def marker_in_deg(data, dName, labels):
    if dName == 'seqA':
        markers ={'Acinar':['PRSS1','CTRB2','REG1A'],
                  'Alpha & Myeloid':['MRC1','MAFB','APOE','HLA-DRA','C1QA','CD14'],
                  'Beta':['NKX6-1','HEPACAM2','DKK3'],
                  'Cancer':['CEACAM5','MSLN','KRT17','LAMC2','KRT16'],
                  'Ductal':['SLC4A4','FXYD2','SPP1','TFF1', 'TFF2','TFF3'],
                  'Endothelial':['PLVAP','VWF'],
                  'Mast':['TPSAB1','CPA3'],
                  'RBCs':['HBB','HBA2'],
                  'T & NK cells':['CD3D','CD2','NKG7'],
                  'pDCs':['IRF7','GZMB'],
                  }
    elif dName == 'seqB':
            markers ={'Acinar':['PRSS1','CTRB2','REG1A'],
                      'Alpha & Myeloid':['MRC1','MAFB','APOE','HLA-DRA','C1QA','CD14'],
                      'Beta':['NKX6-1','HEPACAM2','DKK3'],
                      'Cancer':['CEACAM5','MSLN','KRT17','LAMC2','KRT16'],
                      'Ductal':['SLC4A4','FXYD2','SPP1','TFF1', 'TFF2','TFF3'],
                      'Endothelial':['PLVAP','VWF'],
                      'Fibroblast':['COL1A1','COL3A1','LUM','DCN'],
                      'Mast':['TPSAB1','CPA3'],
                      }
    else:
            markers = {'acinar':'PRSS1',
                       'activated_stellate': ['PDGFRA', 'SPARC'],
                       'alpha':'GCG',
                       'beta':'INS',
                       'delta':'SST',
                       'ductal':['SPP1', 'KRT19'],
                       'endothelial': ['FLT1', 'VWF'],
                       'epsilon':['GHRL', 'ANXA13'],
                       'gamma':'PPY',
                       'macrophage': 'SDS',
                       'mast': 'TPSAB1',
                       'quiescent_stellate': 'RGS5',
                       'schwann': 'SOX10',
                       't_cell':'TRAC',
                       }
    
    for n in range(len(data)):
        print(labels[n])
        sc.tl.rank_genes_groups(data[n], 'cell_labels', method='wilcoxon', key_added = "wilcoxon")
        for i in range(len(data[0].obs['cell_labels'].cat.categories)):
            # print(i)
            ct = data[0].obs['cell_labels'].cat.categories[i]
            degs = sc.get.rank_genes_groups_df(data[n], group=str(ct),
                                        key='wilcoxon',
                                        pval_cutoff=0.05)['names'].squeeze().str.strip().tolist()
            k = list(markers.keys())
            k.remove(ct)
            for j in k:
                m = markers[j]
                c = [i for indx,i in enumerate(m) if list(map(lambda x: x in degs[:50], m))[indx] == True]
                if len(c) != 0:
                    print(ct+'  contains:'+ str(c))
    return()





