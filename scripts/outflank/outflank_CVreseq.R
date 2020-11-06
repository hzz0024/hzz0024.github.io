devtools::install_github("whitlock/OutFLANK")
library(OutFLANK)  # outflank package
library(vcfR)
library(bigsnpr)
library(ggplot2)

bedfile = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.bed"
# this will create a .rds file
snp_readBed(bedfile)
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.rds")
# See how it looks like
str(obj.bigSNP, max.level = 2, strict.width = "cut")
# Get aliases for useful slots
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# Note that most of the algorithms of this package don’t handle missing values. 
# I used snp_fastImputeSimple() to impute missing values of genotyped variants.
#G2 <- snp_fastImputeSimple(G, method = c("random"))
#G3 <- snp_fastImputeSimple(G, method = c("mode"))
SIZE <- c(20, 50, 100, 200, 500)
#SIZE <- c(50)
ind_keeps = list()
# perform snp_clumping on masked dataset at different window size (in KB)
for(i in seq(length(SIZE))){
  size_ = SIZE[i]
  ind.keep_ <- snp_clumping(
    G,
    infos.chr = CHR,
    infos.pos = POS,
    ind.row = rows_along(G),
    S = NULL,
    thr.r2 = 0.2,
    size = size_,
    exclude = NULL,
    ncores = 1
  )
  ind_keeps[[i]] = ind.keep_
}
# 1-5 the index list
index <- ind_keeps[[2]]
# output the index list
write.table(obj.bigSNP$map$marker.ID[index], file = "SNP_thinned.txt", sep = "\t",
            row.names = FALSE, quote = F, col.names=FALSE)
# ind.keepSVD <- snp_autoSVD(
#   G,
#   infos.chr = CHR,
#   infos.pos = POS,
#   thr.r2 = 0.2,
#   size = 50,
# )

################# start outflank #################

#data("sim1a")
#str(sim1a)
str(obj.bigSNP)
# create a new matrix to store the genotypes in bigsnpr
newm = matrix(nrow=length(obj.bigSNP$genotypes[,1]), ncol=length(obj.bigSNP$genotypes[1,]))
for(i in seq(length(obj.bigSNP$genotypes[,1]))){
  newm[i,]=obj.bigSNP$genotypes[i,]
}
# add pop information
#obj.bigSNP$pop <- rep(c(seq(1,5)), each = 6)

obj.bigSNP$pop <- rep(c(1,2,2,2,1), each = 6)
# calculate FST on all the loci in our dataset.
my_fst <- MakeDiploidFSTMat(newm, locusNames = obj.bigSNP$map$physical.pos, popNames = obj.bigSNP$pop)
# Using OutFLANK() function to estimate the parameters on the neutral FST distribution
out_trim <- OutFLANK(my_fst[ind_keeps[[2]],], NumberOfSamples=5, qthreshold = 0.05, Hmin = 0.1)
str(out_trim)

OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)

OutFLANKResultsPlotter(out_trim , withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         TRUE, RightZoomFraction = 0.15, titletext = NULL)

hist(out_trim$results$pvaluesRightTail)
# Using estimated neutral mean FST and df to calculate P-values for all loci
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, 
                                dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)
length(P1$OutlierFlag[P1$OutlierFlag==TRUE])
my_out <- P1$OutlierFlag==TRUE
plot(P1$He, P1$FST, pch=19, col=rgb(0,0,0,0.1))
points(P1$He[my_out], P1$FST[my_out], col="blue")
#my_fst <- MakeDiploidFSTMat(t(sim1a$G), locusNames = sim1a$position, popNames = sim1a$pop)
head(my_fst)
plot(my_fst$He, my_fst$FST)
plot(my_fst$FST, my_fst$FSTNoCorr)
abline(0,1)


dat <- read.vcfR("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.recode.vcf.gz")# , verbose = FALSE 
geno <- extract.gt(dat) # Character matrix containing the genotypes
position <- getPOS(dat) # Positions in bp
chromosome <- getCHROM(dat) # Chromosome information
pop <- rep(c(seq(1,5)), each = 6)
G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
G[geno %in% c("0/0", "0|0")] <- 0
G[geno %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
table(as.vector(G))
dat_fst <- MakeDiploidFSTMat(t(G), locusNames = position, popNames = pop)
head(dat_fst)
#Data checks: Heterozygosity vs. FST
plot(dat_fst$He, dat_fst$FST)

plot(dat_fst$FST, dat_fst$FSTNoCorr)
abline(0, 1, col = "red")

#### Evaluating OutFLANK with trimmed SNPs ####
out_trim <- OutFLANK(dat_fst, NumberOfSamples=5, qthreshold = 0.05, Hmin = 0.1)
str(out_trim)
head(out_trim$results)

OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
OutFLANKResultsPlotter(out_trim , withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         TRUE, RightZoomFraction = 0.15, titletext = NULL)
hist(out_trim$results$pvalues)
hist(out_trim$results$pvaluesRightTail)
#### Run OutFLANK on all loci ####
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar,
                                dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)

################# start pcadapt #################
install.packages("pcadapt")
library("pcadapt")
path_to_file <- "./pcadapt.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
# Scree plot
x <- pcadapt(input = filename, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
plot(x, option = "screeplot")
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
poplist.names <- c(rep("Louisiana_SL", 6),rep("Louisiana_Sel", 6), rep("DelBay_Sel_NEH", 6),rep("Ches_Sel_DEBY", 6), rep("DelBay_CS", 6))
print(poplist.names)
plot(x, option = "scores", pop = poplist.names)
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)
# Computing the test statistic based on PCA
x <- pcadapt(filename, K = 3)
# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
# Manhattan Plot
plot(x , option = "manhattan")
