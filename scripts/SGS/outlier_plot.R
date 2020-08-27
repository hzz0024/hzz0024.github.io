#######################
#  Delta_p plot       #
#######################

library(ggplot2)
# load the neutral datasets
files <- list.files('.', pattern = "*.txt")
# load the observation dataset
observation_file = 'obs_deltap.output'
# obtain the deltap from obs dataset
observations = read.delim(observation_file, header = TRUE, sep = "\t", dec = ".")
observations = observations$deltaP
# create a matrix for neu deltap 
deltaP_matix = c()
for(file in files){
  dat <- read.delim(file, header = TRUE, sep = "\t", dec = ".")
  deltaP_matix = cbind(deltaP_matix,dat$deltaP)
}
dim(deltaP_matix)

# create variables for quantile values
mins = c()
maxs = c()
mids = c()
# create tags for deltap comparsion
tags = c()
# loop over snps
num_snp = 3664
for(snp in seq(1,num_snp)){
  delta_Ps = deltaP_matix[snp,]
  q0 = unname(quantile(delta_Ps, probs=0.25))
  mins = c(mins, q0)
  q1 = unname(quantile(delta_Ps, probs=0.99))
  maxs = c(maxs, q1)
  q2 = unname(quantile(delta_Ps, probs=0.5))
  mids = c(mids, q2)
  
  observation = observations[snp]
  if(observation>q1){
    tags = c(tags, 1)
  }
  else{
    tags = c(tags, 0)
  }
}

DATA = data.frame(X=seq(1,num_snp), MIN=mins, MAX=maxs, MID=mids,OBS=observations)

sp <- ggplot(DATA, aes(x=X, y=MID)) +
  geom_point(size=.5)+
  geom_point(aes(x=X, y=OBS),size=.5,color='red')+
  # draws the range bars
  geom_errorbar(data=DATA, aes(ymin=MIN, ymax=MAX), width=.001,color='yellow',alpha=.2)
# add x and y-axis titles
sp + scale_x_continuous(name="SNP", limits=c(0, 400)) +
  scale_y_continuous(name="absolute deltap", limits=c(0, 0.5))

#######################
#  Delta_p vs start p #
#######################

##################### reveal the relationship between deltap and start p, first plot #####################
#setwd("/Volumes/cornell/DelBay19_Hopper/permutation/4_deltap_plot/deltap_vs_p_plot1")
# Delta p plot
library(ggplot2)
# load the neutral datasets
files <- list.files('.', pattern = "*.txt")
# load the dataset
ch_file = 'CH_ref_98_ch_doMAF_filter.mafs.extracted'
ref_file = 'REF_ref_98_ref_doMAF_filter.mafs.extracted'
# obtain the deltap from obs dataset
ch = read.delim(ch_file, header = TRUE, sep = "\t", dec = ".")
ref = read.delim(ref_file, header = TRUE, sep = "\t", dec = ".")
p0 = ref$knownEM
p1 = ch$knownEM
DATA = data.frame(MIN=p0, MAX=p1)
DATA = DATA[order(DATA$MIN),]
num_snp = 3664
DATA$X = seq(1, num_snp)

sp <- ggplot(DATA, aes(x=X, y=MIN)) +
  geom_point(size=.5)+
  geom_point(aes(x=X, y=MAX),size=.5,color='red')+
  # draws the range bars
  geom_errorbar(data=DATA, aes(ymin=MIN, ymax=MAX), width=.001,color='yellow',alpha=.8)
# add x and y-axis titles
sp + scale_x_continuous(name="SNP", limits=c(0, 400)) +
  scale_y_continuous(name="Allele frequency", limits=c(0, 1)) +
  labs(title = "Minor allele changes for 386 SNP outliers",
       subtitle = "ref allele = black, ch allele = red, actual deltap = yellow")

##################### reveal the relationship between deltap and start p, second plot #####################
#setwd("/Volumes/cornell/DelBay19_Hopper/permutation/4_deltap_plot/deltap_vs_p_plot2")
library(ggplot2)
# load the dataset
ch_file = 'CH_maf0.05_pctind0.7_cv30.mafs.extracted'
ref_file = 'REF_maf0.05_pctind0.7_cv30.mafs.extracted'
#deltap_file = 'obs_deltap.output'
# obtain the deltap from obs dataset
ch = read.delim(ch_file, header = FALSE, sep = "\t", dec = ".")
ref = read.delim(ref_file, header = FALSE, sep = "\t", dec = ".")
deltap = ch$V6 - ref$V6
p0 = ref$V6
p1 = ch$V6
DATA = data.frame(p=p0, delta_p=deltap)
DATA = DATA[order(DATA$p),]
num_snp = 96
sp <- ggplot(DATA, aes(x=p, y=delta_p)) +
  geom_point(size=.5)
# add x and y-axis titles
sp + scale_x_continuous(name="p", limits=c(min(DATA$p), max(DATA$p))) +
  scale_y_continuous(name="Deltap (absolute values)", limits=c(min(DATA$delta_p), max(DATA$delta_p))) +
  labs(title = "Deltap against reference allele p for the observation data",
       subtitle = "note deltap ranges from 0-0.5")



# load the neutral datasets
file1 <- list.files('.', pattern = "*.txt")
file2 <- list.files('.', pattern = "*.extracted")
# create a matrix for neu deltap and p
deltaP_matix = c()
for(f1 in file1){
  dat1 <- read.delim(f1, header = TRUE, sep = "\t", dec = ".")
  deltaP_matix = cbind(deltaP_matix,dat1$deltaP)
}
p_matix = c()
for(f2 in file2){
  dat2 <- read.delim(f2, header = TRUE, sep = "\t", dec = ".")
  p_matix = cbind(p_matix,dat2$knownEM)
}
dim(deltaP_matix)
dim(p_matix)
# create variables for quantile values
mid1 = c()
mid2 = c()
# create tags for deltap comparsion
tags = c()
# loop over snps
num_snp = 386
for(snp in seq(1,num_snp)){
  delta_Ps = deltaP_matix[snp,]
  p_s = p_matix[snp,]
  q2 = unname(quantile(delta_Ps, probs=0.5))
  ref_p = unname(quantile(p_s, probs=0.5))
  mid1 = c(mid1, q2)
  mid2 = c(mid2, ref_p)
}

DATA = data.frame(X=seq(1,num_snp), neu_deltap = mid1, neu_p = mid2)

sp <- ggplot(DATA, aes(x=neu_p, y=neu_deltap)) +
  geom_point(size=.5)
sp + scale_x_continuous(name="p", limits=c(0, 0.5)) +
  scale_y_continuous(name="Deltap (absolute values)", limits=c(0, 0.1)) +
  labs(title = "Deltap against reference allele p for the neutral data",
       subtitle = "Deltap and p are medium (50% quantile) values for each SNP, note deltap ranges from 0-0.1")
