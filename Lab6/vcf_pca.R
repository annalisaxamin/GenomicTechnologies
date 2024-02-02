# load the package 
library("SNPRelate") 

# open the VCF file and convert it to GDS format 
vcf_file <- "AllSamplesMerged_chrM.vcf.gz" 
snpgdsVCF2GDS(vcf_file, "AllSamplesMerged_chrM.gds", method="biallelic.only") 

# load the GDS file 
genofile <- snpgdsOpen("AllSamplesMerged_chrM.gds") 

# compute the PCA
variantsPca <- snpgdsPCA(genofile, autosome.only=F)

# choose 10 colors, one per sample 
jColors <- c('chartreuse3', 'cornflowerblue', 'darkgoldenrod1',
             'peachpuff3', 'mediumorchid2', 'turquoise3', 'wheat4',
             'slategray2','red','blue') 

set.seed(33)
# compute a small random position offset for each sample (jitter) 
x_jitter = jitter(variantsPca$eigenvect[,1], .1) 
y_jitter = jitter(variantsPca$eigenvect[,2], .1) 

# plot the PCA 
plot(x_jitter, y_jitter, col=jColors, xlab = "PC1", ylab = "PC2", pch=16) 

# plot the labels of each sample 
text(x_jitter, y_jitter, labels=variantsPca$sample.id, 
     pos=runif(10, min=1, max=4)) # random uniform

