devtools::install_github("whitlock/OutFLANK")
library(OutFLANK)  # outflank package
library(vcfR)
library(bigsnpr)
library(ggplot2)

bedfile = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.noLM.bed"
# this will create a .rds file
obj.bed <- bed(bedfile)
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.noLM.rds")
# See how it looks like
str(obj.bigSNP, max.level = 2, strict.width = "cut")
# Get aliases for useful slots
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
# Note that most of the algorithms of this package don’t handle missing values. 
# I used snp_fastImputeSimple() to impute missing values of genotyped variants.
G2 <- snp_fastImputeSimple(G, method = c("random"))
#G4 <- snp_fastImputeSimple(G)
#G3 <- snp_fastImputeSimple(G, method = c("mode"))
SIZE <- c(20, 50, 100, 200, 500)
ind_keeps = list()

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

index <- ind_keeps[[5]]
write.table(obj.bigSNP$map$marker.ID[index], file = "SNP_thinned.txt", sep = "\t",
            row.names = FALSE, quote = F, col.names=FALSE)
#ind.keep <- snp_autoSVD(
#  G2,
#  infos.chr = CHR,
#  infos.pos = POS,
#  thr.r2 = 0.2,
#  size = 200,
#)










data("sim1a")
str(sim1a)
sim1a$pop
sim1a$envi
# calculate Fst
my_fst <- MakeDiploidFSTMat(t(sim1a$G), locusNames = sim1a$position, popNames = sim1a$pop)
head(my_fst)
plot(my_fst$He, my_fst$FST)
plot(my_fst$FST, my_fst$FSTNoCorr)
abline(0,1)
data("which_pruned")
head(which_pruned)



DB_1_vcf <- read.vcfR("DB_1_prune.recode.vcf" , verbose = FALSE )
geno <- extract.gt(DB_1_vcf) # Character matrix containing the genotypes
position <- getPOS(DB_1_vcf) # Positions in bp
chromosome <- getCHROM(DB_1_vcf) # Chromosome information
pop <- c(rep(1, 5),rep(2, 7))
G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
G[geno %in% c("0/0", "0|0")] <- 0
G[geno %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
table(as.vector(G))
my_fst <- MakeDiploidFSTMat(t(G), locusNames = position, popNames = pop)
head(my_fst)
#Data checks: Heterozygosity vs. FST
plot(my_fst$He, my_fst$FST)

plot(my_fst$FST, my_fst$FSTNoCorr)
abline(0,1)

#### Evaluating OutFLANK with trimmed SNPs ####
out_trim <- OutFLANK(my_fst, NumberOfSamples=2, qthreshold = 0.1, Hmin = 0.1)
str(out_trim)
