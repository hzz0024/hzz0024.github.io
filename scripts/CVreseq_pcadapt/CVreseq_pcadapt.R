library(OutFLANK)
library(vcfR)
library(bigsnpr)
library(ggplot2)
library(plyr)
#library(pcadapt)
setwd("~/Documents/Ryan_workplace/CVreseq_clumping")
bedfile = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.bed"
# this will create a .rds file
snp_readBed(bedfile)

# run pcadapt
bedfile = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.bed"
bed.obj <- read.pcadapt(bedfile, type = "bed")
x <- pcadapt(bed.obj, K = 3)

bed_pcadapt(bed.obj, newpc$u)
  

pcadapt()


# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.rds")
# See how it looks like
str(obj.bigSNP, max.level = 2, strict.width = "cut")
# Get aliases for useful slots
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos


Gm = matrix(nrow=length(obj.bigSNP$genotypes[,1]), ncol=length(obj.bigSNP$genotypes[1,]))
for(i in seq(length(obj.bigSNP$genotypes[,1]))){
  Gm[i,]=obj.bigSNP$genotypes[i,]
}

a_freq2 <- rowSums(Gm)/(2*ncol(Gm))
keep_loci = which(a_freq2 > 0.01 & a_freq2 < 0.99)
G_sub2 <- G_sub[keep_loci,]
training <- list(G = G_sub2, position = final_df$position, 
                 chromosome = final_df$chrom)


newpc<-snp_autoSVD(G=G,infos.chr = CHR,infos.pos = POS)

newpc <- big_SVD(G, k=3) 
snp_pcadapt(G, U.row = newpc$u[,1])

toMatrix <- function(G){
  Gm = matrix(nrow=length(G[,1]), ncol=length(G[1,]))
  for(i in seq(length(G[,1]))){
    Gm[i,]=G[i,]
  }
  return(Gm)
}

remove.packages("bigsnpr")
install.packages("/Users/ryan/Downloads/bigsnpr", repos = NULL, type="source")
library(bigsnpr)
setwd("~/Documents/Ryan_workplace/CVreseq_clumping")
obj.bigSNP <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.rds")
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos

newpc <- big_SVD(G, k=3) 
test <- snp_gc(snp_pcadapt(G, U.row = newpc$u[,1:2]))
result <- -predict(test,log10=T)

na_idx = c(102120,432510,452617,470761,492280,492281,506875,696154,696155,947788,1045640,1281332,1285977,1292155,1387378,1387379,1489399,1502236,1502237,1566254,1571089,1704001,1704049,1705835,1705836,1706681,1751737,1933633,2015470,2040086,2040108,2086156,2503937,2503981,2556820,2556822,2569353,2650550,2730314,2730315,2730316,2771093,2785790,2842281,2862515,2873343,2978492,3023826,3023827,3076380,3076381,3076382,3126955,3134065,3261668,3428280,3636806,3977518,3977519,4240259,4272699,4374683,4374745,4635250,4635251,4735050,4836671,4836672,4991820,4991990,4991991,4991992,5000999,5095498,5104855,5128374,5189110,5189111,5410894)
#na_idx = c(102120,432510,452617,470761,492280,492281,506875,696154,696155,947788,1045640,1281332,1285977,1292155,1387378,1387379,1489399,1502236,1502237,1566254,1571089,1704001,1704049,1705835,1705836,1706681,1751737,1933633,2015470,2040086,2040108,2086156,2503937,2503981,2556820,2556822,2569353,2650550,2730314,2730315,2730316,2771093,2785790,2842281,2862515,2873343,2978492,3023826,3023827,3076380,3076381,3076382,3126955,3134065,3261668,3428280,3636806,3977518,3977519,4240259,4272699,4374683,4374745,4635250,4635251,4735050,4836671,4836672,4991820,4991990,4991991,4991992,5000999,5095498,5104855,5128374,5189110,5189111,5410894,5635537,5965927,5986034,6004178,6025697,6025698,6040292,6229571,6229572,6481205,6579057,6814749,6819394,6825572,6920795,6920796,7022816,7035653,7035654,7099671,7104506,7237418,7237466,7239252,7239253,7240098,7285154,7467050,7548887,7573503,7573525,7619573,8037354,8037398,8090237,8090239,8102770,8183967,8263731,8263732,8263733,8304510,8319207,8375698,8395932,8406760,8511909,8557243,8557244,8609797,8609798,8609799,8660372,8667482,8795085,8961697,9170223,9510935,9510936,9773676,9806116,9908100,9908162,10168667,10168668,10268467,10370088,10370089,10525237,10525407,10525408,10525409,10534416,10628915,10638272,10661791,10722527,10722528,10944311)
Gm = toMatrix(G)
Gm = Gm[,-na_idx]
dim(Gm)
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
newpc1 <- big_SVD(G_coded, k=3) 
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u))
result <- -predict(test,log10=T)

tmp = as.integer(CHR[-na_idx])*104144668+POS[-na_idx]
(maxp <- max(result, na.rm=TRUE))
plot(tmp, result,
     ylim= c(0, maxp*1.19), pch=19, col=rgb(0,0,0,0.5), xaxs="i", yaxs="i",
     cex=0.5, ylab= "-log10(P) PCAdapt B.P.", xaxt="n", type="h", bty="l", main="B) PCAdapt Best Practice", adj=0)
plot_layers(y_head=maxp*100, y_arrows=c(maxp*1.18, maxp*1.05), thisSim=thisim)

library('pcadapt')
filename = 'SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.bed'
G_all <- read.pcadapt(filename, type = "bed")
pca_all <- pcadapt(G_all,K=3)#, data.type="genotype")
all_data <- -log10(pca_all$pvalues)
tmp = as.integer(CHR)*104144668+POS
(maxp <- max(all_data, na.rm=TRUE))
plot(tmp, all_data,
     ylim= c(0, maxp*1.19), pch=19, col=rgb(0,0,0,0.5), xaxs="i", yaxs="i",
     cex=0.5, ylab= "-log10(P) PCAdapt Naive", xaxt="n", type="h", bty="l", main="B) PCAdapt Best Practice", adj=0)
plot_layers(y_head=maxp*100, y_arrows=c(maxp*1.18, maxp*1.05), thisSim=thisim)





Gm1 = Gm[,which_pruned][,10000:11000]
Gm2 = Gm[,which_pruned][,10000:11000]
sum(is.na(Gm1))
sum(is.na(Gm2))
G_coded <- add_code256(big_copy(Gm1,
                                type="raw"),
                       code=bigsnpr:::CODE_012)

G_coded1 <- add_code256(big_copy(Gm2,
                                type="raw"),
                       code=bigsnpr:::CODE_012)


obj.svd<-big_SVD(G_coded, fun.scaling = snp_scaleBinom(), k = 10)
test <- snp_gc(snp_pcadapt(G_coded1, U.row = obj.svd$u))
result <- -predict(test,log10=T)



#test case
test <- snp_attachExtdata()
GG <- test$genotypes
matrix = toMatrix(GG)

newpc<-snp_autoSVD(G=GG, infos.chr = test$map$chromosome, infos.pos = test$map$physical.pos)
which_pruned = attr(newpc, 'subset')
Gm1 = matrix[,which_pruned]
sum(is.na(Gm1))
G_coded <- add_code256(big_copy(Gm1,
                                type="raw"),
                       code=bigsnpr:::CODE_012)


obj.svd<-big_SVD(G_coded, fun.scaling = snp_scaleBinom(), k = 10)
result1 <- -predict(snp_gc(snp_pcadapt(GG, U.row = obj.svd$u)),log10=T)

obj.svd<-big_SVD(GG, fun.scaling = snp_scaleBinom(), k = 10)
result2 <- -predict(snp_gc(snp_pcadapt(GG, U.row = obj.svd$u)),log10=T)

plot(test$map$physical.pos, result1)
plot(test$map$physical.pos, result2)


filename = '/Users/ryan/Downloads/TTT_RecombinationGenomeScans-ed399e4d22f4b9907775c21e1f88cdb3ee2cb7c2/11109_Invers_VCFallFILT.vcf'
vcf <- read.vcfR(filename)
ends=c(0,  seq(50000, 450000, by=50000))
dim(vcf@gt)
vcf@fix[,"CHROM"] <- NA
POS <- as.numeric(vcf@fix[,"POS"])
for (i in 1:(length(ends)-1)){
  cond <- POS >= ends[i] &  POS < ends[i+1]
  print(c(ends[i], ends[i+1], sum(cond)))
  vcf@fix[cond,"CHROM"] = i
}
my_ord <- order(as.numeric(vcf@fix[,"POS"]))
vcf2 <- vcf
vcf2 <- vcf[my_ord,]
rm(my_ord, cond, ends)
invloc <- which(vcf2@fix[,"POS"]==320000)
check <- levels(factor(vcf2@gt[invloc ,]))
check
table(vcf2@gt[invloc ,],vcf2@gt[invloc ,])
table(vcf2@gt[invloc+1 ,],vcf2@gt[invloc+1 ,])
if((("0|0" %in% check)&("2|2" %in% check))|
   (("0|0" %in% check)&("0|2" %in% check))|
   (("0|0" %in% check)&("2|0" %in% check))|
   (("2|2" %in% check)&("0|2" %in% check))|
   (("2|2" %in% check)&("2|0" %in% check))
){
  vcf2@gt[invloc ,] <- sub("2|2", "1|1", vcf2@gt[invloc ,], fixed=TRUE)
  vcf2@gt[invloc ,] <- sub("2|0", "1|0", vcf2@gt[invloc ,], fixed=TRUE)
  vcf2@gt[invloc ,] <- sub("0|2", "0|1", vcf2@gt[invloc ,], fixed=TRUE)
}
if((("2|2" %in% check)&("3|3" %in% check))|
   (("2|2" %in% check)&("2|3" %in% check))|
   (("2|2" %in% check)&("3|2" %in% check))|
   (("3|3" %in% check)&("2|3" %in% check))|
   (("3|3" %in% check)&("3|2" %in% check))
){
  vcf2@gt[invloc ,] <- sub("2|2", "0|0", vcf2@gt[invloc ,], fixed=TRUE)
  vcf2@gt[invloc ,] <- sub("3|3", "1|1", vcf2@gt[invloc ,], fixed=TRUE)
  vcf2@gt[invloc ,] <- sub("2|3", "0|1", vcf2@gt[invloc ,], fixed=TRUE)
  vcf2@gt[invloc ,] <- sub("3|2", "1|0", vcf2@gt[invloc ,], fixed=TRUE)
}

filename = '/Users/ryan/Downloads/TTT_RecombinationGenomeScans-ed399e4d22f4b9907775c21e1f88cdb3ee2cb7c2/11109_Invers_VCFallFILT.vcf'


vcf <- vcf2
vcf@gt[invloc ,]
rm(vcf2,  check)
geno <- vcf@gt[,-1] 
position <- getPOS(vcf) # Positions in bp
chromosome <- getCHROM(vcf) # Chromosome information
G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
a_freq <- rowSums(G)/(2*ncol(G))
G_sub = G
final_df0 <- data.frame(position=position, chrom=chromosome, 
                        a_freq_old = a_freq, 
                        vcf_ord = 1:length(position),
                        unique = paste(position, format(a_freq, scientific=FALSE), sep="_")
)
a_freq2 <- rowSums(G_sub)/(2*ncol(G_sub))
final_df0$a_freq_final <- a_freq2
keep_loci = which(final_df0$a_freq_final > 0.01 & final_df0$a_freq_final < 0.99)
final_df0$keep_loci <- FALSE
final_df0$keep_loci[keep_loci] <- TRUE
sum(is.na(final_df0$keep_loci))

G_sub2 <- G_sub[keep_loci,]

final_df <- final_df0[which(final_df0$keep_loci),]
training <- list(G = G_sub2, position = final_df$position, 
                 chromosome = final_df$chrom)

options(bigstatsr.typecast.warning = FALSE)
G_coded <- add_code256(big_copy(t(training$G),
                                type="raw"),
                       code=bigsnpr:::CODE_012)
# puts it in the raw format and stores likelihood genotype probability
newpc<-snp_autoSVD(G=G_coded,
                   infos.chr = as.integer(training$chromosome),
                   infos.pos = training$position)
which_pruned =  attr(newpc, which="subset")
training$G_pruned <- training$G[which_pruned,]#

G_p_coded <- add_code256(big_copy(t(training$G_pruned),
                                type="raw"),
                       code=bigsnpr:::CODE_012)

newpc<-snp_autoSVD(G=G_p_coded,
                   infos.chr = as.integer(training$chromosome)[which_pruned],
                   infos.pos = training$position[which_pruned])

training$G_coded <- G_coded
test <- snp_gc(snp_pcadapt(training$G_coded, U.row = newpc$u[,1]))
final_df$pcadapt_3.0.4_PRUNED_log10p <- -predict(test,log10=T)
plot(final_df$position,final_df$pcadapt_3.0.4_PRUNED_log10p)

thisim = 10945
tmp = as.integer(final_df$chrom)*449947+final_df$position
(maxp <- max(final_df$pcadapt_3.0.4_PRUNED_log10p, na.rm=TRUE))
plot(tmp, final_df$pcadapt_3.0.4_PRUNED_log10p,
     ylim= c(0, maxp*1.19), pch=19, col=rgb(0,0,0,0.5), xaxs="i", yaxs="i",
     cex=0.5, ylab= "-log10(P) PCAdapt B.P.", xaxt="n", type="h", bty="l", main="B) PCAdapt Best Practice", adj=0)
plot_layers(y_head=maxp*100, y_arrows=c(maxp*1.18, maxp*1.05), thisSim=thisim)


filename = '/Users/ryan/Downloads/TTT_RecombinationGenomeScans-ed399e4d22f4b9907775c21e1f88cdb3ee2cb7c2/G_all.bed'
G_all <- read.pcadapt(filename, type = "bed")
pca_all <- pcadapt(G_all,K=3)#, data.type="genotype")
final_df$pcadapt_3.0.4_ALL_log10p <- -log10(pca_all$pvalues)

plot(tmp, final_df$pcadapt_3.0.4_ALL_log10p,
     ylim= c(0, maxp*1.19), pch=19, col=rgb(0,0,0,0.5), xaxs="i", yaxs="i",
     cex=0.5, ylab= "-log10(P) PCAdapt B.P.", xaxt="n", type="h", bty="l", main="B) PCAdapt Best Practice", adj=0)
plot_layers(y_head=maxp*100, y_arrows=c(maxp*1.18, maxp*1.05), thisSim=thisim)

