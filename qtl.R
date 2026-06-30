library(qtl2)

setwd('/home/mathieu/postdoc_heasley/wine_subclade_spores/rqtl2/')
geno_codes <- c(1,2)
names(geno_codes) <- c('AA','BB')
write_control_file('config.yaml', 
                   geno_file = 'geno.csv', 
                   pmap_file = 'pmap.csv', 
                   gmap_file = 'gmap.csv', 
                   pheno_file = 'pheno.csv',
                   crosstype = 'haploid', 
                   geno_codes=geno_codes, overwrite=TRUE)

cross <- read_cross2('config.yaml', quiet = FALSE)
gprob <- calc_genoprob(cross)
out <- scan1(gprob, cross$pheno)
peaks <- find_peaks(out, cross$pmap)

write.csv(out, file = 'lod.csv')
write.csv(peaks, file = 'peaks.csv')
