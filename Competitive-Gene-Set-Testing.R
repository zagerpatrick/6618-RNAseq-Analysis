## 01.0-Package Import ---------------------------------------------
setwd("P:/Experiments/021_rnaseq_(Christin)/6618_Analysis_Patrick")

library(magrittr)

# Load Bioconductor Packages
library(edgeR)
library(limma)
library(org.Mm.eg.db)

# Load Functions File
source("6618_Func.R")

## 02.0-Data Import ---------------------------------------------

setwd("P:/Experiments/021_rnaseq_(Christin)/6618_Analysis_Patrick")
load(file = "6618_Analysis_Pre_Process.RData")

# Load in Mouse Gene Sets 

setwd("./Gene_Set_Database")
load("mouse_c2_v5p2.rdata")
load("mouse_c5_v5p2.rdata")


## 03.0-Model Setup ---------------------------------------------


anno.entrez <- select(org.Mm.eg.db,
                  keys = rownames(voom_transformed), # rownames
                  keytype= "ENSEMBL", # rownames are ENSEMBL identifiers
                  columns= c("ENTREZID"), uniqueRows=TRUE)

anno.trim.entrez <- aggregate(anno.entrez['ENTREZID'], anno.DE['ENSEMBL'], 
                          FUN = function(X) paste(unique(X), collapse=", "))

voom_gene_set <- voom_transformed

rownames(voom_gene_set) <- anno.trim.entrez$ENTREZID

idx.c2 <- ids2indices(Mm.c2, id=as.character(unlist(anno.trim.entrez['ENTREZID'])))

idx.c5 <- ids2indices(Mm.c5, id=as.character(unlist(anno.trim.entrez['ENTREZID'])))


## 04.0-Competitive Gene Set Testing ---------------------------------------------

cam.curated.noCNVvsCNV_GalRC <- camera(voom_gene_set, idx.c2, design, contrast=contr.matrix[,5])
head(cam.curated.noCNVvsCNV_GalRC,5)

cam.GO.noCNVvsCNV_GalRC <- camera(voom_gene_set, idx.c5, design, contrast=contr.matrix[,5])
head(cam.GO.noCNVvsCNV_GalRC,5)


library(xlsx)
write.xlsx(cam.curated.noCNVvsCNV_GalRC, "cam.curated-GAL.noCNV.RC-GAL.CNV.RC.xlsx")
write.xlsx(cam.GO.noCNVvsCNV_GalRC, "cam.GO-GAL.noCNV.RC-GAL.CNV.RC.xlsx")



barcodeplot(voomed.fitted$t[,3], index=idx.c5$GO_NEGATIVE_REGULATION_OF_AMINE_TRANSPORT, 
            index2=idx.c5$GO_NEGATIVE_REGULATION_OF_AMINE_TRANSPORT, main="WTvsIHH_CNVRC")

