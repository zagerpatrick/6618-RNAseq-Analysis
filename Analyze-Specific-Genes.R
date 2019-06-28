## 01.0-Package Import ---------------------------------------------

library(magrittr)

# Load Bioconductor Packages
library(edgeR)
library(limma)
library(org.Mm.eg.db)

# Load Functions File
source("6618_Func.R")

## 02.0-Data Import ---------------------------------------------

load(file = "6618_Analysis_Pre_Process.RData")

## 03.0-Model Setup ---------------------------------------------

temp <- rlog.norm.counts
better.names <- treatfactor <- factor(paste(factor.matrix[,2], 
                                factor.matrix[,3], factor.matrix[,4], sep="."))

colnames(temp) <- better.names

anno.human.counts <- select(org.Mm.eg.db,
                  keys = rownames(temp), # rownames
                  keytype= "ENSEMBL", # rownames are ENSEMBL identifiers
                  columns= c("SYMBOL")) # what to return

human.counts <- merge(as.data.frame(temp), 
                      anno.human.counts, by.x = "row.names", by.y = "ENSEMBL")

ggplot(human.counts, aes(x = factor(year), y = case, fill = code))+
  geom_dotplot(binaxis = 'y', stackdir = 'center',
               position = position_dodge())

ggplot(human.counts, aes(x = 1, y = mpg)) +
  geom_dotplot(binaxis = "y", stackdir = "center")



mtcars