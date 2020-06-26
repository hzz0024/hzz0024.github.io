### step 1 reorder the major and minor allele ###
ref_file = 'ch_ref_98_ref_doMAF_filter.mafs'
template = read.delim(ref_file, header = TRUE, sep = "\t", dec = ".")
snp = data.frame(chr=template$chromo, chr=template$position)
write.table(snp, file = "snp_list.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
# then use python extract.py for data formatting

### step 2 calculate the theta

n=251*2
a1 = 0
for(i in seq(1,n-1)){
  a1 = a1 + 1/i
}

S = 2.185/85
M_hat = S/a1

### step 2 create null model

library(gtools)
n = 100
k = 25
M = 0.0037

draw_distribution <- function(n,k,M){
  tajima_points = c()
  for(i in seq(1:(n-1))){
    upper_right =  M * (1/i + 1/(n-i))
    tajima_points = c(tajima_points, upper_right)
  }
  tajima_points = tajima_points/sum(tajima_points) 
  points = c()
  for(i in seq(1:(n-1))){
    p = i/n
    upper_left = factorial(n)/(factorial(n - k)*factorial(k)) * p^k * (1-p)^(n-k)
    upper_right =  M * (1/i + 1/(n-i))/sum(tajima_points)
    res = upper_left * upper_right / (1/(n-1))
    points = c(points, res)
  }
  points=points/sum(points)
  return(list("tajima" = tajima_points, "points" = points))
}

### actual deltap ###
null_distribution <- function(n, k, obs_delta){
  N1 = draw_distribution(n, k, 0.0037)
  num_sample = 10000
  sample_p1 = sample(seq(1:(n-1))/n, num_sample, prob=N1$points, replace=TRUE)
  sample_p2 = sample(seq(1:(n-1))/n, num_sample, prob=N1$points, replace=TRUE)
  delta_ps = c()
  for(j in seq(1,num_sample)){
    p1 = sample_p1[j]
    p2 = sample_p2[j]
    delta_p = (p1-p2)
    delta_ps = c(delta_ps, delta_p)
  }
  hist(delta_ps, xlab="delta_p", main = "Null distribution of deltap from 10000 iterations ")
  abline(v=obs_delta, col='red')
}

filename = 'ch_ref_98_ref_doMAF_filter.mafs.extracted'
obs_file = 'obs_deltap.output'

dat <- read.delim(filename, header = TRUE, sep='\t')
obs_dat <- read.delim(obs_file, header = TRUE, sep='\t')

for(i in seq(1,dim(dat)[1])){
  null_distribution(n=dat$nInd[i]*2, k=floor(dat$nInd[i]*dat$knownEM[i]), obs_delta=obs_dat$deltaP[i])
}


### abs deltap ###
draw_distribution <- function(n,k,M){
  tajima_points = c()
  for(i in seq(1:(n-1))){
    upper_right =  M * (1/i + 1/(n-i))
    tajima_points = c(tajima_points, upper_right)
  }
  tajima_points = tajima_points/sum(tajima_points) 
  points = c()
  for(i in seq(1:(n-1))){
    p = i/n
    upper_left = factorial(n)/(factorial(n - k)*factorial(k)) * p^k * (1-p)^(n-k)
    upper_right =  M * (1/i + 1/(n-i))/sum(tajima_points)
    res = upper_left * upper_right / (1/(n-1))
    points = c(points, res)
  }
  points=points/sum(points)
  return(list("tajima" = tajima_points, "points" = points))
}

null_distribution <- function(n, k, obs_delta){
  N1 = draw_distribution(n, k, 0.01351005) #chr5 global theta
  num_sample = 10000
  sample_p1 = sample(seq(1:(n-1))/n, num_sample, prob=N1$points, replace=TRUE)
  sample_p2 = sample(seq(1:(n-1))/n, num_sample, prob=N1$points, replace=TRUE)
  delta_ps = c()
  for(j in seq(1,num_sample)){
    p1 = sample_p1[j]
    p2 = sample_p2[j]
    delta_p = abs(p1-p2)
    delta_ps = c(delta_ps, delta_p)
  }
  hist(delta_ps, xlab="delta_p", main = "Null distribution of deltap from 10000 iterations ")
  abline(v=obs_delta, col='red')
  return(length(delta_ps[delta_ps>obs_delta])/num_sample)
}

filename = 'ch_ref_98_ref_doMAF_filter.mafs.extracted'
obs_file = 'obs_deltap.output'

dat <- read.delim(filename, header = TRUE, sep='\t')
obs_dat <- read.delim(obs_file, header = TRUE, sep='\t')

for(i in seq(1,dim(dat)[1])){
  null_distribution(n=dat$nInd[i]*2, k=floor(dat$nInd[i]*dat$knownEM[i]), obs_delta=obs_dat$deltaP[i])
}

p_values = c()
for(i in seq(1,dim(dat)[1])){
  p_value = null_distribution(n=dat$nInd[i]*2, k=floor(dat$nInd[i]*dat$knownEM[i]), obs_delta=obs_dat$deltaP[i])
  p_values = c(p_values, p_value)
}

out = data.frame(chromo=dat$chromo, position=dat$position, p_value=p_values)
write.table(out, file = "p_value_list.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
