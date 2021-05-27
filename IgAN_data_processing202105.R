# load raw data
rawdata=read.csv("IgAN_3143cell_raw_counts.txt.gz", sep="\t", row.names=1)
cell.ids=colnames(rawdata)


# Pagoda2 clustering

library(velocyto.R)
library(pagoda2)

ldat = read.loom.matrices("raw_data.loom")
emat<- ldat$spliced
colnames(emat)=gsub("IgA_10p:|_S.*", "", colnames(emat))

emat.f=emat[, cell.ids]

rownames(emat.f)=make.unique(rownames(emat.f))

r.f <- Pagoda2$new(emat.f,modelType='plain',log.scale=T)

r.f$adjustVariance(plot=T,gam.k=10)

r.f$calculatePcaReduction(nPcs=50,n.odgenes=3e3)

set.seed(123)
r.f$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=30,verbose=T)

r.f$makeKnnGraph(k=20,type='PCA',center=T,distance='cosine')
r.f$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')

# Pagoda2 cluster result
pagoda.cluster= r.f$clusters$PCA[[1]] 
save(pagoda.cluster, file="pagoda.cluster.RData")

# Pagoda2 tSNE layout
r.f$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,min.group.size=10,shuffle.colors=F,mark.cluster.cex=2,alpha=0.5,main='clusters')


# Seurat analysis (V2)

library(Seurat)
library(dplyr)

dat.seurat <- CreateSeuratObject(raw.data = rawdata, min.cells = 3, min.genes = 200, 
                                 project = "glo")

mito.genes <- grep(pattern = "^mt-", x = rownames(x = dat.seurat@data), value = TRUE)
percent.mito <- Matrix::colSums(dat.seurat@raw.data[mito.genes, ])/Matrix::colSums(dat.seurat@raw.data)

dat.seurat <- AddMetaData(object = dat.seurat, metadata = percent.mito, col.name = "percent.mito")

pdf("plot_plates_stats_vln.pdf", width = 20, height = 5 )
VlnPlot(object = dat.seurat, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

#Normalizing the data
dat.seurat <- NormalizeData(object = dat.seurat, normalization.method = "LogNormalize", 
                            scale.factor = 500000)

#load the Pagoda2 cluster result into Seurat object

dat.seurat@ident=pagoda.cluster[names(dat.seurat@ident)]
table(dat.seurat@ident)

#identify marker genes for each cluster

dat.pagoda.markers <- FindAllMarkers(object = dat.seurat, only.pos = TRUE, min.pct = 0.2)


#identify differentially expressed genes for each cluster

dat.seurat.genotype=rep("", nrow(dat.seurat@meta.data))
dat.seurat.genotype[grep("318|316|313|315|325|322", rownames(dat.seurat@meta.data) )]="IgAN"
dat.seurat.genotype[grep("317|314|326|321", rownames(dat.seurat@meta.data) )]="Control"

dat.seurat@meta.data$geno <- dat.seurat.genotype

dat.seurat@meta.data$celltype.geno <- paste0(dat.seurat@ident, "_", 
                                             dat.seurat@meta.data$geno)

dat.seurat <- StashIdent(dat.seurat, save.name = "celltype")
dat.seurat <- SetAllIdent(dat.seurat, id = "celltype.geno")

nk=length(unique(pagoda.cluster))
pagoda.k.de.gene.list=list()
for(i in 1:nk){
  print(i)
  pagoda.k.de.gene.list[[i]]=FindMarkers(dat.seurat, ident.1 = paste0(i, "_IgAN"), ident.2 = paste0(i, "_Control"),print.bar = FALSE)
}


