import os
import glob
import pickle
import pandas as pd
import numpy as np
from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from arboreto.algo import genie3
from numpy.core.umath_tests import inner1d
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune, prune2df, df2regulons
from pyscenic.aucell import aucell
import seaborn as sns

DATA_FOLDER="./"
RESOURCES_FOLDER="../resources"
DATABASE_FOLDER = "../databases/"
SCHEDULER="123.122.8.24:8786"
DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "hg19*.feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'hh_total_tfs.txt')
SC_EXP_FNAME = os.path.join(RESOURCES_FOLDER, "Human.pse20.txt")
REGULONS_FNAME = os.path.join(DATA_FOLDER, "regulons.p")
MOTIFS_FNAME = os.path.join(DATA_FOLDER, "motifs.csv")
AUC_FNAME=os.path.join(DATA_FOLDER, "aucell.csv")
Co_FNAME=os.path.join(DATA_FOLDER, "adjacencies.csv")

ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0).T
ex_matrix.shape
tf_names = load_tf_names(MM_TFS_FNAME)
db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
	return os.path.basename(fname).split(".")[0]

dbs=[RankingDatabase(fname=fname,name=name(fname)) for fname in db_fnames]

#Phase1:co-expression module
adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)
adjacencies.to_csv(Co_FNAME)
adjacencies = pd.read_csv("/home/jingjingw/Jingjingw/Project/2018-MH-new/2019-1-18-TotalFig2-new/2_SCENIC/Human_total/adjacencies.csv")
modules = list(modules_from_adjacencies(adjacencies, ex_matrix))

#Phase2: RcisTarget [Prune modules for targets with cis regulatory footprints]
# Calculate a list of enriched motifs and the corresponding target genes for all modules.
with ProgressBar():
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

# Create regulons from this table of enriched motifs.
regulons = df2regulons(df)
# Save the enriched motifs and the discovered regulons to disk.
df.to_csv(MOTIFS_FNAME)
with open(REGULONS_FNAME, "wb") as f:
    pickle.dump(regulons, f)

regulons = prune(dbs, modules, MOTIF_ANNOTATIONS_FNAME)
    
#Phase3: AUCell    
auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
auc_mtx.to_csv(AUC_FNAME)
sns_plot=sns.clustermap(auc_mtx, figsize=(12,12))   
sns_plot.savefig("sns.clustermap.png") 
    
