library(SNPRelate)
library(rrBLUP)
library(vcfR)
library(data.table)

setwd('/home/mathieu/postdoc_heasley/wine_subclade_spores/gwas/')

vcf.path <- 'wine_subclade.biallelic.vcf.gz'
gds.path <- 'wine_subclade.gds'
vcf <- read.vcfR(vcf.path)
vcf.fix <- as.data.frame(getFIX(vcf))
F_segr_parents <- sapply(read.table('F_segr_parents.txt'), as.logical)


## SNPRelate PCA

#snpgdsVCF2GDS(vcf.path, gds.path, method='biallelic.only', ignore.chr.prefix='chromosome')

gds <- snpgdsOpen(gds.path)

S <- read.gdsn(index.gdsn(gds, "sample.id"))

prunedLD <- snpgdsLDpruning(gds, slide.max.bp=50000, ld.threshold = 0.2)

pca <- snpgdsPCA(gds, snp.id=unlist(prunedLD), num.thread=2)

df.pca <- data.frame(sample.id = pca$sample.id, 
                     EV1 = pca$eigenvect[,1],
                     EV2 = pca$eigenvect[,2],
                     EV3 = pca$eigenvect[,3],
                     EV4 = pca$eigenvect[,4],
                     EV5 = pca$eigenvect[,5],
                     EV6 = pca$eigenvect[,6],
                     stringsAsFactors = FALSE)

plot(df.pca$EV2, df.pca$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

snpgdsClose(gds)

## Format for rrBLUP

geno.raw.path <- 'rrblup_geno.raw'
geno.raw <- read.table(geno.raw.path, header=TRUE, stringsAsFactors=FALSE)
S <- geno.raw[, 1]
geno.matrix <- t(geno.raw[, 7:ncol(geno.raw)])
colnames(geno.matrix) <- S

geno.matrix[geno.matrix == 0] <- -1   # Homozygote 1
geno.matrix[geno.matrix == 1] <- 0    # Heterozygote
geno.matrix[geno.matrix == 2] <- 1    # Homozygote 2
geno.matrix[is.na(geno.matrix)] <- NA

markers <- rownames(geno.matrix)
chromosomes <- as.integer(vcf.fix$CHROM)

position <- as.integer(vcf.fix$POS)
map.df <- data.frame(marker=markers, chrom=chromosomes, pos=position, stringsAsFactors=FALSE)

final.geno <- cbind(map.df, geno.matrix)
final.geno <- final.geno[F_segr_parents,]

## Import phenotypic data
GC_LOGISTIC <- read.csv('/home/mathieu/postdoc_heasley/experiments/growth_curves/HU_pheno/tables/GC_LOGISTIC.csv')
GC_LOGISTIC <- GC_LOGISTIC[GC_LOGISTIC$dataset=='wt_ypd_sc' & GC_LOGISTIC$medium=='YPD' & GC_LOGISTIC$HU==0 & GC_LOGISTIC$OD=='log_OD600', ]
rownames(GC_LOGISTIC) <- GC_LOGISTIC$strain

pheno <- GC_LOGISTIC[S, c('strain', 'mu')]
pheno$tri8 <- pheno$strain=='YJM948'

## rrBLUP

gwas <- GWAS(pheno, final.geno, fixed='tri8', K=NULL, n.PC=0,
             min.MAF=0.05, n.core=1, P3D=TRUE, plot=TRUE)
gwas.path <- 'rrblup_gwas.csv'
write.csv(gwas, file = gwas.path)

## loop for score threshold

random.gwas.permut <- list()

for(i in c(1:500)){
  random.pheno <- pheno
  random.pheno$mu <- sample(random.pheno$mu)
  random.gwas.iter <- GWAS(random.pheno, final.geno, fixed='tri8', K=NULL, n.PC=0,
                           min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
  random.gwas.iter$iter = i
  random.gwas.permut[[i]] = random.gwas.iter

  print(i)
  
}

random.gwas.permut <- rbindlist(random.gwas.permut)

random.gwas.path <- 'random.rrblup_gwas.csv'
write.csv(random.gwas.permut, file = random.gwas.path)