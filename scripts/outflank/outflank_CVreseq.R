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
ind_keeps = list()
# perform snp_clumping on masked dataset
for(i in seq(length(SIZE))){
  size_ = SIZE[i]
  ind.keep_ <- snp_clumping(
    G2,
    infos.chr = CHR,
    infos.pos = POS,
    ind.row = rows_along(G2),
    S = NULL,
    thr.r2 = 0.2,
    size = size_,
    exclude = NULL,
    ncores = 1
  )
  ind_keeps[[i]] = ind.keep_
}

index <- ind_keeps[[4]]
write.table(obj.bigSNP$map$marker.ID[index], file = "SNP_thinned.txt", sep = "\t",
            row.names = FALSE, quote = F, col.names=FALSE)
#ind.keep <- snp_autoSVD(
#  G2,
#  infos.chr = CHR,
#  infos.pos = POS,
#  thr.r2 = 0.2,
#  size = 200,
#)

################# start outflank #################

# data("sim1a")
# str(sim1a)
# sim1a$pop
# sim1a$envi
# # calculate Fst
# my_fst <- MakeDiploidFSTMat(t(sim1a$G), locusNames = sim1a$position, popNames = sim1a$pop)
# head(my_fst)
# plot(my_fst$He, my_fst$FST)
# plot(my_fst$FST, my_fst$FSTNoCorr)
# abline(0,1)
# data("which_pruned")
# head(which_pruned)

dat <- read.vcfR("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.thinned.dom_wild.recode.vcf" , verbose = FALSE )
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