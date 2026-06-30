library(DESeq2)
library(Rsubread)

base.dir <- '/home/mathieu/postdoc_heasley/rna_seq/'
setwd(base.dir)

strains.info.path = paste(base.dir, 'script/strains_info.csv', sep='')
strains.info <- read.csv(strains.info.path)
strains.info <- strains.info[strains.info$project=='MH',]
row.names(strains.info) <- strains.info$sample

## FeatureCounts

bam.list <- as.character(lapply(strains.info$sample, (function(x) paste(base.dir, 'data/bam/', x,'.dedup.bam', sep=''))))
gtf.path <- paste(base.dir, '/data/ref/featurecounts.gtf', sep='')

fc <- featureCounts(bam.list, annot.ext=gtf.path,
              isGTFAnnotationFile = TRUE, GTF.featureType = 'gene', GTF.attrType = 'Name')

fc.path <- paste(base.dir, 'tables/counts.csv', sep='')
#write.csv(fc$counts, file = fc.path)
# re-import counts
counts <- read.csv(fc.path, row.names = 1)

fc.stats.path <- paste(base.dir, 'tables/fc_stats.csv', sep='')
write.csv(fc$stat, file = fc.stats.path)


## DESeq
col.data <- strains.info[,c('rep','IC','y_prime_kb')]
# Normalize independent variables to 1
col.data$y_prime_kb = col.data$y_prime_kb/max(col.data$y_prime_kb)
col.data$IC = col.data$IC/max(col.data$IC)

S <- row.names(col.data)
#row.names(col.data) <- S
#counts <- fc$counts
colnames(counts) <- lapply(colnames(counts), function(x) strsplit(x, '.', fixed=TRUE)[[1]][[1]])
counts <- counts[, S]

# Interaction between CN and IC

dds.complete <- DESeqDataSetFromMatrix(countData = counts, 
                                       colData = col.data, 
                                       design = ~ IC * y_prime_kb)
dds.complete <- DESeq(dds.complete)

for(x in resultsNames(dds.complete)){
  res <- results(dds.complete, name = x, independentFiltering = FALSE)
  
  deseq.results.path <- paste(base.dir, 'tables/deseq.complete.', x, '.csv', sep='')
  write.csv(res, deseq.results.path)
  
}

# CN alone

S.IC0 <- row.names(col.data[col.data$IC==0, ])

dds.single <- DESeqDataSetFromMatrix(countData = counts[, S.IC0], 
                                     colData = col.data[S.IC0, ], 
                                     design = ~ y_prime_kb)
dds.single <- DESeq(dds.single)

for(x in resultsNames(dds.single)){
  res <- results(dds.single, name = x, independentFiltering = FALSE)
  
  deseq.results.path <- paste(base.dir, 'tables/deseq.single_yprime.', x, '.csv', sep='')
  write.csv(res, deseq.results.path)
}

# IC for low yprime alone
  
S.Y20 <- row.names(col.data[col.data$y_prime_kb<0.05, ])

dds.single2 <- DESeqDataSetFromMatrix(countData = counts[, S.Y20], 
                                      colData = col.data[S.Y20, ], 
                                      design = ~ IC)
dds.single2 <- DESeq(dds.single2)

for(x in resultsNames(dds.single2)){
  res <- results(dds.single2, name = x, independentFiltering = FALSE)
  
  deseq.results.path <- paste(base.dir, 'tables/deseq.single_IC.', x, '.csv', sep='')
  write.csv(res, deseq.results.path)
}

### Re-do the analysis on unmasked reference to look at TPE.

## FeatureCounts

bam.list.um <- as.character(lapply(strains.info$sample, (function(x) paste(base.dir, 'data/bam/', x,'_unmasked.dedup.bam', sep=''))))
gtf.path <- paste(base.dir, '/data/ref/featurecounts.gtf', sep='')

fc.um <- featureCounts(bam.list.um, annot.ext=gtf.path,
                       isGTFAnnotationFile = TRUE, GTF.featureType = 'gene', GTF.attrType = 'Name')

fc.path.um <- paste(base.dir, 'tables/counts.um.csv', sep='')
#write.csv(fc.um$counts, file = fc.path.um)
# re-import counts
counts.um <- read.csv(fc.path.um, row.names = 1)

fc.stats.path.um <- paste(base.dir, 'tables/fc_stats.um.csv', sep='')
#write.csv(fc.um$stat, file = fc.stats.path.um)


## DESeq
col.data <- strains.info[,c('rep','IC','y_prime_kb')]
# Normalize independent variables to 1
col.data$y_prime_kb = col.data$y_prime_kb/max(col.data$y_prime_kb)
col.data$IC = col.data$IC/max(col.data$IC)

S <- row.names(col.data)
#row.names(col.data) <- S
#counts <- fc$counts
colnames(counts.um) <- lapply(colnames(counts.um), function(x) gsub('_unmasked.dedup.bam', '', x))
counts.um <- counts.um[, S]

# CN alone

S.IC0 <- row.names(col.data[col.data$IC==0, ])

dds.single.um <- DESeqDataSetFromMatrix(countData = counts.um[, S.IC0], 
                                     colData = col.data[S.IC0, ], 
                                     design = ~ y_prime_kb)
dds.single.um <- DESeq(dds.single.um)

for(x in resultsNames(dds.single.um)){
  res <- results(dds.single.um, name = x, independentFiltering = FALSE)
  
  deseq.results.path.um <- paste(base.dir, 'tables/deseq.single_yprime.um.', x, '.csv', sep='')
  write.csv(res, deseq.results.path.um)
}