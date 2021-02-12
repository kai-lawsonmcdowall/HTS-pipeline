if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")
BiocManager::install("GenomicFeatures")
BiocManager::install('EnhancedVolcano')
install.packages("readr")

library(tximport)
library(GenomicFeatures)
library(readr)
library(EnhancedVolcano)


#importing the salmon quantification: 


dir <- "/localhome/bs20kjlm/Precision Medicine/HTS project/Data"
list.files(dir)
samples <- read.table(file.path(dir, "samplesL20.csv"),header=TRUE, sep = ",")
samples <- samples[-c(7, 8, 9, 10),] # removes additional rows created
files <- file.path(dir, samples$Sample)
files
names(files) <- samples$Sample
tx2gene <- read.csv("/localhome/bs20kjlm/Precision Medicine/HTS project/Data/gencode.v35.tx2gene.csv") #import their file
head(tx2gene)
a<- gsub("\\..*","",tx2gene[,1]) # these two lines remove the version. 
b<- gsub("\\..*","",tx2gene[,2])
c<-cbind(a,b)
colnames(c)=colnames(tx2gene)
tx2gene <- as.data.frame(c)
head(tx2gene)
txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreTxVersion = T,ignoreAfterBar = T)

# Using DESeq2: 

#downloaded the tar.gz file and then imported it into the working directory and installed via
#Tools -> Install Packages -> Package Archive File (.tar.gz) and then in install 
#link to install the tar.gz file: https://bioconductor.org/packages/release/bioc/html/DESeq2.html


library(DESeq2)
?DESeqDataSetFromTximport 
dds <- DESeqDataSetFromTximport(txi, colData =samples, design = ~Condition)
dds <- DESeq(dds)
dds$condition <- factor(dds$Condition, levels = c("PAR","I10"))
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
res <- results(dds)
res

#most differntially expressed genes
head(res[order(res$padj),], 10)

head(res) #see how the results currently stand
nms <- row.names(res)
nms15 <- substr(nms,1,15)
rownames(res) <- nms15
res <- na.omit(res)
head(res) 
# the results following the removal of Ensembl gene identifiers at end of "ENSG" 
# omission of res results that had NA values when expression values are very low. 
summary(res)

plotDispEsts(dds, main="Dispersion plot") #A dispersion plot
#DEseq2 assumes a negative binomial distribution
# dispersion is how muhc variance deviates from mean

plotMA(dds)

#volcano plot

alpha <- 0.1 # Threshold on the adjusted p-value

cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
plot(res$log2FoldChange, 
     -log10(res$padj), 
     col=cols, panel.first=grid(),
     main="Volcano plot of differentially expressed genes", 
     xlab="Effect size: log 2(fold-change)", 
     ylab="-log 10(adjusted p-value)", 
     ylim=c(0,15),
     pch=16, cex=0.6)

abline(v=0)
abline(v=c(-1,1), col="red")
abline(h=-log10(alpha), col="red")
gn.selected <- abs(res$log2FoldChange) > 2.5 & res$padj < alpha 


up <- rownames(res)[which(res$log2FoldChange>1)]
down <- rownames(res)[which(res$log2FoldChange<(1*-1))]

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", ylim = c(0,15), col = "grey"))

with(subset(res , padj<.1), points(log2FoldChange, -log10(pvalue), pch=20,col="grey"))
with(res[which(rownames(res) %in% up),], points(log2FoldChange, -log10(pvalue), pch=20, col="#0F7F7F"))
with(res[which(rownames(res) %in% down),], points(log2FoldChange, -log10(pvalue), pch=20, col="lightcoral"))
abline(h=-log10(alpha), col="black")

res$log2FoldChange

#clustering and heatmaps: 

rld <- rlogTransformation(dds) #log transform the data
head(assay(rld))



# this code connects to the ensembl server

library(biomaRt)
BiocManager::install("topGO")

bm <- useMart("ensembl") # # biomaRt connects to the Ensmebl resource, biomart, which gives info on Ensembl related info.
bm <- useDataset("hsapiens_gene_ensembl", mart=bm) # use the human database 
EG2GO <- getBM(attributes=c("ensembl_gene_id","external_gene_name", "go_id"), mart=bm, useCache = FALSE)
# fetches the different gene id's


head(EG2GO) # EG2GO is a data frame where ENSG corresponds to HGNC names
gid2GO <- by (EG2GO$go_id, EG2GO$ensembl_gene_id, function(x) as.character(x)) # takes some time
# "by" aggregates all the HGNC terms by gene into a list of lists  format (list of Go terms categorized by ENSG identifier) 
head(gid2GO, 20) #prints the first 20 values of this list of lists
gL <- res$padj
gL # this is the gene list which contains a list of adjusted differental expression p values. 
names(gL) <- row.names(res) #

names(gL)  #A list of the ENSG terms (gene list)


library(topGO)

topDiffGenes <- function(allScore) {return (allScore<0.01)} # returns p>0.01

GOdata <- new ("topGOdata", ontology="BP", allGenes=gL, 
               annot = annFUN.gene2GO, gene2GO=gid2GO, 
               description = "All diff genes", geneSel=topDiffGenes)


genelist <- as.data.frame(gL)
allGO <- genesInTerm(GOData)
#Getting all annotated genes for a GO ID
allGO["GO:0019221"]
RetrivedGenes <- lapply(allGO,function(x) x[x %in% genelist$V1]) # where INT.GENES$V1 is my list of interesting genes
# Your significant genes for GO:0019221
RetrivedGenes[["GO:0019221"]]




# The above code creates a topGOdata object, with ontology = Biological Process,
# The (ENSG) genes that need processing and their respective P-values
# AnnFUN.gene2Go are used when we provide our own annaotations as genes-to-GO (i.e. gid2GO)
# genesel is a gene selection process where only genes with a p value of 0.01  or less are selected

#Finally, we can perform the GO enrichment tests: 


test.stat <- new ("classicCount", testStatistic=GOFisherTest, name="Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
allRes <- GenTable(GOdata,classic=resultFisher, topNodes=5)




#Annotated : number of genes in database which are annotated with the GO-term.
#Significant : number of genes belonging to your input which are annotated with the GO-term.
#Expected : show an estimate of the number of genes a node of size Annotated would have if the significant genes were to be randomly selected from the gene universe.



#creating a readable excel file
install.packages("xlsx")
library(xlsx)

write.table(allRes, "/localhome/bs20kjlm/allRes.csv", row.names = FALSE)


# Write the first data

getwd
print(res)


dir2 <- "/localhome/bs16saj/precision_med_2020/HTS_project_data"
list.files(dir2)

