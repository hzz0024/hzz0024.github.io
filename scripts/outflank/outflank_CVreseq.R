devtools::install_github("whitlock/OutFLANK")
library(OutFLANK)  # outflank package
library(vcfR)
library(bigsnpr)

bedfile = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.noLM.bed"
#snp_readBed(bedfile)

obj.bed = bed(bedfile)

ind.keep <- bed_clumping(
  obj.bed,
  ind.row = rows_along(obj.bed),
  S = NULL,
  thr.r2 = 0.2,
  size = 5,
  exclude = NULL,
  ncores = 1
)




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
