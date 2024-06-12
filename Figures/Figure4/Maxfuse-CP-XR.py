import numpy as np
import pandas as pd
from scipy.io import mmread

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (6, 4)

from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

import anndata as ad
import scanpy as sc
import maxfuse as mf

import seaborn as sns

cp = sc.read_h5ad('/mnt/scratch1/Touchstone-Analysis/CP_StJude_anndata.h5ad')
xr = sc.read_10x_mtx('/mnt/scratch1/Touchstone_data/new_data/TOUCHSTONE_STJ_XR_FFPE_PR_1_UOA_C_R1/cell_feature_matrix')

xr_meta = pd.read_csv('/mnt/scratch1/Touchstone-Analysis/MetaAnnotation_TOUCHSTONE_STJ_XR_FFPE_PR_1_UOA_C_R1.csv')
xr_meta = xr_meta.set_index('cellid')

xr = xr[xr.obs.index.isin(xr_meta.index)]


xr.obs['celltype'] = xr_meta['celltype']

correspondence = pd.read_csv('/mnt/scratch1/Touchstone-Analysis/protein_gene_conversion.csv')


rna_protein_correspondence = []

for i in range(correspondence.shape[0]):
    curr_protein_name, curr_rna_names = correspondence.iloc[i]
    if curr_protein_name not in cp.var_names:
        continue
    if curr_rna_names.find('Ignore') != -1:  # some correspondence ignored eg. protein isoform to one gene
        continue
    curr_rna_names = curr_rna_names.split('/')  # eg. one protein to multiple genes
    for r in curr_rna_names:
        if r in xr.var_names:
            rna_protein_correspondence.append([r, curr_protein_name])

rna_protein_correspondence = np.array(rna_protein_correspondence)

# Columns rna_shared and protein_shared are matched.
# One may encounter "Variable names are not unique" warning,
# this is fine and is because one RNA may encode multiple proteins and vice versa.
rna_shared = xr[:, rna_protein_correspondence[:, 0]].copy()
protein_shared = cp[:, rna_protein_correspondence[:, 1]].copy()

# Make sure no column is static
mask = (
    (rna_shared.X.toarray().std(axis=0) > 0.1) 
    & (protein_shared.X.std(axis=0) > 0.1)
)
rna_shared = rna_shared[:, mask].copy()
protein_shared = protein_shared[:, mask].copy()
print([rna_shared.shape,protein_shared.shape])

# process rna_shared
sc.pp.normalize_total(rna_shared)
sc.pp.log1p(rna_shared)
sc.pp.scale(rna_shared)

# plot UMAP of rna cells based only on rna markers with protein correspondence

sc.pp.neighbors(rna_shared, n_neighbors=15)
sc.tl.umap(rna_shared)
#sc.pl.umap(rna_shared, color='celltype')

rna_shared = rna_shared.X.copy()
protein_shared = protein_shared.X.copy()

# process all RNA features
sc.pp.normalize_total(xr)
sc.pp.log1p(xr)
sc.pp.highly_variable_genes(xr, n_top_genes=5000)
# only retain highly variable genes
xr = xr[:, xr.var.highly_variable].copy()
sc.pp.scale(xr)


# plot UMAPs of rna cells based on all active rna markers

sc.pp.neighbors(xr, n_neighbors=15)
sc.tl.umap(xr)
#sc.pl.umap(xr, color='celltype')

# make sure no feature is static
rna_active = xr.X
protein_active = cp.X
rna_active = rna_active[:, rna_active.std(axis=0) > 1e-5] # these are fine since already using variable features
protein_active = protein_active[:, protein_active.std(axis=0) > 1e-5] # protein are generally variable


# inspect shape of the four matrices
print(rna_active.shape)
print(protein_active.shape)
print(rna_shared.shape)
print(protein_shared.shape)


# call constructor for Fusor object
# which is the main object for running MaxFuse pipeline
fusor = mf.model.Fusor(
    shared_arr1=rna_shared,
    shared_arr2=protein_shared,
    active_arr1=rna_active,
    active_arr2=protein_active,
    labels1=np.array(xr.obs['celltype']),
    labels2=np.array(cp.obs['leiden'])
)

fusor.split_into_batches(
    max_outward_size=5000,
    matching_ratio=4,
    metacell_size=2,
    verbose=True
)

fusor.construct_graphs(
    n_neighbors1=15,
    n_neighbors2=15,
    svd_components1=40,
    svd_components2=15,
    resolution1=2,
    resolution2=2,
    # if two resolutions differ less than resolution_tol
    # then we do not distinguish between then
    resolution_tol=0.1,
    verbose=True
)


fusor.find_initial_pivots(
    wt1=0.3, wt2=0.3,
    svd_components1=10, svd_components2=12
)

fusor.refine_pivots(
    wt1=0.3, wt2=0.3,
    svd_components1=40, svd_components2=None,
    cca_components=25,
    n_iters=1,
    randomized_svd=False, 
    svd_runs=1,
    verbose=True
)


fusor.filter_bad_matches(target='pivot', filter_prop=0.5)

pivot_matching = fusor.get_matching(order=(2, 1),target='pivot')

lv1_acc = mf.metrics.get_matching_acc(matching=pivot_matching, 
    labels1=np.array(xr.obs['celltype']), 
    labels2=np.array(cp.obs['leiden']),
    order = (2,1)
)
lv1_acc

fusor.propagate(
    svd_components1=40, 
    svd_components2=None, 
    wt1=0.7,
    wt2=0.7,
)

fusor.filter_bad_matches(
    target='propagated',
    filter_prop=0.3
)


full_matching = fusor.get_matching(order=(2, 1), target='full_data')



pd.DataFrame(list(zip(full_matching[0], full_matching[1], full_matching[2])), 
             columns = ['mod1_indx', 'mod2_indx', 'score'])


rna_cca, protein_cca_sub = fusor.get_embedding(
    active_arr1=fusor.active_arr1,
    active_arr2=fusor.active_arr2[full_matching[1],:] # cells in codex remained after filtering
)


dim_use = 15 # dimensions of the CCA embedding to be used for UMAP etc

cca_adata = ad.AnnData(
    np.concatenate((rna_cca[:, :dim_use], protein_cca_sub[:, :dim_use]), axis=0),
    dtype=np.float32
)

cca_adata.obs['data_type'] = ['rna'] * rna_cca.shape[0] + ['protein'] * protein_cca_sub.shape[0]
cca_adata.obs['cell_type'] = list(np.concatenate((np.array(xr.obs['celltype']), np.array(cp.obs['leiden'])[full_matching[1]]), axis = 0))




sc.pp.neighbors(cca_adata, n_neighbors=15)
sc.tl.umap(cca_adata)
sc.pl.umap(cca_adata, color='data_type', save = "~/MaxFuse-UMAP-XR-CP.png")


