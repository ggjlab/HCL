load("/home/ggj/HCA/RData/pbmc/old/AdultBoneMarrowCD34P1_Seurat.RData")
setwd("/home/ggj/NEW/DifferentiationForce/Data/20200111_embryo/wt/paga/")
dim(pbmc@data)
#17364 11781
setwd("/home/ggj/NEW/DifferentiationForce/Data/20200111_embryo/wt/paga")
aa<-as.data.frame(as.matrix(pbmc@assays$RNA@counts))
write.csv(aa[,rownames(anno)],file = "./wt.dge.csv",quote = F)
anno<-FetchData(pbmc,vars = "ident")
write.csv(anno,file = "./wt_embryo.anno.csv",quote = F)
table(pbmc@ident)

vargene<-pbmc@var.genes
gene<-rownames(pbmc@raw.data)
genefilter<-gene%in%vargene
write.csv(genefilter,file = "./genefilter.csv",quote = F,row.names = F)


anno<-FetchData(pbmc,vars="ident")
anno$shuchu<-ifelse(anno$ident%in%c(1,2,6,12,14,18,21),"FALSE","TRUE")
write.csv(anno$shuchu,file="nonimmunecell.csv",quote=F,row.names=F)

anno$shuchu<-ifelse(anno$ident%in%c(1,2,6,12,14,18,21),"TRUE","FALSE")
write.csv(anno$shuchu,file="immunecell.csv",quote=F,row.names=F)

#####################################python
import numpy as np
import scanpy.api as sc
import pandas as pd
import os
import pandas as pd

os.chdir("/home/ggj/NEW/DifferentiationForce/Data/20200111_embryo/wt/paga/")
adata=sc.read_csv("wt.dge.csv",delimiter=',').transpose()
adata.var_names = pd.read_csv('gene.csv', header=None)[0]

datause= pd.read_table("wt.dge.csv",sep=",",index_col=0)
adata=sc.AnnData(datause.T)


mito_genes = [name for name in adata.var_names if name.startswith('mt-')]
#adata[:, mito_genes]=0



adata.obs['cluster']= pd.read_csv('wt_embryo.anno.csv',sep=",",header=None)[0].values
adata.obs['type']= pd.read_csv('wt_embryo.anno.csv',sep=",",header=None)[1].values
adata.obs['batch']= pd.read_csv('wt_embryo.anno.csv',sep=",",header=None)[2].values
adata.obs['cluster1']= pd.read_csv('wt_embryo.anno.csv',sep=",",header=None)[3].values



sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=0)
adata.obs['n_counts'] = adata.X.sum(axis=1)


sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
adata.raw = adata

adata.write('./wt1.h5ad', compression='gzip')


gene_filter1=pd.read_csv("genefilter.csv")
gene_filter=gene_filter1['x'].values
sc.pl.filter_genes_dispersion(filter_result)
adata = adata[:, gene_filter]

cc=adata.var_names
cc=cc.values
bb=list(set(cc).intersection(set(gene_filter)))
len(bb)


filter_result = sc.pp.filter_genes_dispersion(adata.X, min_mean=0.01, max_mean=15, min_disp=0.4)
import collections
collections.Counter(filter_result.gene_subset)
#Counter({False: 14819, True: 1750})
sc.pl.filter_genes_dispersion(filter_result)
adata = adata[:, filter_result.gene_subset]


sc.pp.regress_out(adata, ['n_counts'])

## scale the data
sc.pp.scale(adata, max_value=10)


### PCA
sc.tl.pca(adata, n_comps=50)
sc.pl.pca_loadings(adata)
# visualize
adata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
sc.pl.pca_scatter(adata, color='COL1A1')
# PC
sc.pl.pca_variance_ratio(adata, log=True,  show=50,n_pcs=50)
## 25
adata





sc.pp.neighbors(adata, n_neighbors=10,n_pcs=25)
sc.tl.louvain(adata, resolution=1)
sc.tl.tsne(adata,use_fast_tsne=True,n_jobs=20,perplexity=200,n_pcs=25)
sc.pl.tsne(adata, color='louvain',size=8,legend_loc="on data")


adata.obs['type']=adata.obs['type'].astype('category')


sc.pl.tsne(adata, color='louvain',size=8)
sc.pl.tsne(adata, color='tissue',size=8,legend_loc="on data")
sc.pl.tsne(adata, color='tissue',size=8)
adata.write('./wt_cluster.h5ad', compression='gzip')



sc.pp.neighbors(adata,n_pcs=25)
sc.tl.umap(adata, min_dist=0.1)
sc.pl.umap(adata, color='type', title='UMAP', legend_loc='on data', legend_fontsize=5)
sc.pl.umap(adata, color='louvain', title='UMAP', legend_loc='on data', legend_fontsize=5)



#sc.tl.paga(adata, groups='type', model='v1.0')
sc.tl.paga(adata, groups='type')

sc.pl.paga(
    adata,
    layout='fr',
    threshold=0.01,
    fontsize=8,
    node_size_scale=1,
    node_size_power=0.7,
    max_edge_width=0.7)

import matplotlib.pyplot as pl

sc.tl.paga(adata, groups='type')

sc.tl.draw_graph(adata,init_pos="paga",layout="fa",maxiter=500)
sc.pl.draw_graph(adata, color='type',title='Force Atlas 2', legend_loc='on data', legend_fontsize=5,palette=sc.pl.palettes.default_20)
#,save="iter1000.pdf")
sc.pl.tsne(adata, color='louvain',size=8,palette=sc.pl.palettes.godsnot_64,legend_loc="on data")

sc.pl.draw_graph(adata, color='cluster',title='Force Atlas 2', legend_loc='on data', legend_fontsize=5,palette=sc.pl.palettes.default_20,size=8)

adata.write('./wt_paga.h5ad', compression='gzip')



adata1=adata.copy()
sc.tl.paga(adata1, groups='louvain')



################################ remove immune cells
adata.obs['donor_tf']= pd.read_csv("/home/ggj/NEW/DifferentiationForce/Data/20200111_embryo/wt/paga/nonimmunecell.csv",sep=",",header=None)[0].values
#adata.obs['cluster_num']= pd.read_csv('/home/ggj/NEW/HCA/SPRING/PBCD34/cluster.csv',sep=",",header=None)[0].values

cluster = adata[adata.obs["donor_tf"]]
#6118 × 1750 

sc.pp.neighbors(cluster,n_pcs=25)
sc.tl.umap(cluster, min_dist=0.1)

sc.tl.paga(cluster, groups='type')
sc.pl.paga(
    cluster,
    layout='fr',
    threshold=0.01,
    fontsize=8,
    node_size_scale=1,
    node_size_power=0.7,
    max_edge_width=0.7)

import matplotlib.pyplot as pl
sc.tl.draw_graph(cluster,init_pos="paga",layout="fa",maxiter=500)
sc.pl.draw_graph(cluster, color='type',  legend_loc='on data',legend_fontsize=5,palette=sc.pl.palettes.default_26,size=8)

cluster.write('./wt_nonimmune_paga.h5ad', compression='gzip')


################################  immune cells
adata.obs['donor_tf']= pd.read_csv("/home/ggj/NEW/DifferentiationForce/Data/20200111_embryo/wt/paga/immunecell.csv",sep=",",header=None)[0].values
immune = adata[adata.obs["donor_tf"]]
#5663 × 1750 

sc.pp.neighbors(immune,n_pcs=25)
sc.tl.umap(immune, min_dist=0.1)

sc.tl.paga(immune, groups='type')
sc.pl.paga(
    immune,
    layout='fr',
    threshold=0.01,
    fontsize=8,
    node_size_scale=1,
    node_size_power=0.7,
    max_edge_width=0.7)

import matplotlib.pyplot as pl
sc.tl.draw_graph(immune,init_pos="paga",layout="fa",maxiter=500)
sc.pl.draw_graph(immune, color='type',  legend_loc='on data',legend_fontsize=5,palette=sc.pl.palettes.default_26,size=8)
immune.write('./wt_immune_paga.h5ad', compression='gzip')
