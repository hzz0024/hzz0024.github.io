### create null model
library(gtools)
library(hash)

filename = 'CHR_maf0.05_pctind0.7_cv30.mafs'
challenge = 'CH_maf0.05_pctind0.7_cv30.mafs'
obs_file = 'obs_deltap_cv30.output'
outputfile = 'p_values.txt'

dat2 <- read.delim(challenge, header = TRUE, sep='\t')
n2s <- dat2$nInd*2
dat <- read.delim(filename, header = TRUE, sep='\t')
obs_dat <- read.delim(obs_file, header = TRUE, sep='\t') # from ch - ref
dat$n2 = n2s

#dat = dat[dat$chromo==5,]
#obs_dat= obs_dat[obs_dat$chromo==5,]

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
    #upper_left = factorial(n)/(factorial(n - k)*factorial(k)) * p^k * (1-p)^(n-k)
    upper_left = exp(lfactorial(n)-(lfactorial(n - k)+lfactorial(k))) * p^k * (1-p)^(n-k)
    upper_right =  M * (1/i + 1/(n-i))/sum(tajima_points)
    res = upper_left * upper_right / (1/(n-1))
    points = c(points, res)
  }
  points=points/sum(points)
  return(list("tajima" = tajima_points, "points" = points))
}

null_distribution <- function(n, k, n2){
  theta = 0.0091 #global theta
  N1 = draw_distribution(n, k, theta) 
  n2_k = sample(seq(1:(n-1))/n, 1, prob=N1$points, replace=TRUE) * n2
  N2 = draw_distribution(n2, n2_k, theta)
  
  num_sample = 10000
  #set.seed(1)
  sample_p1 = sample(seq(1:(n-1))/n, num_sample, prob=N1$points, replace=TRUE)
  #set.seed(0)
  #sample_p2 = sample(seq(1:(n-1))/n, num_sample, prob=N1$points, replace=TRUE)
  sample_p2 = sample(seq(1:(n2-1))/n2, num_sample, prob=N2$points, replace=TRUE)
  
  delta_ps = c()
  for(j in seq(1,num_sample)){
    p1 = sample_p1[j]
    p2 = sample_p2[j]
    delta_p = (p2-p1) # note here I focus on the actual deltap, which requires two sided test, the order is same as ch - ref
    delta_ps = c(delta_ps, delta_p)
  }
  #return(length(delta_ps[delta_ps>obs_delta])/num_sample)
  return(delta_ps)
}

# main 
cnt = 0
dic <- hash()
p_values = c()
sink(outputfile)
for(i in seq(1,dim(dat)[1])){
  s = paste0(i,'/',dim(dat)[1])
  message(s,"\r",appendLF=FALSE)
  n=dat$nInd[i]*2
  k=floor(dat$nInd[i]*2*dat$knownEM[i])
  n2=dat$n2[i]
  
  obs_delta=obs_dat$deltaP[i]
  
  key = paste0(n, ' ', k, ' ', n2)
  if(has.key(key, dic)){
    delta_ps <- dic[[key]]
  }
  else{
    delta_ps = null_distribution(n=n, k=k, n2=n2)
    dic[[key]] <- delta_ps
  }
  
  if(FALSE){
     DT = data.frame(delta_ps)
     ggplot(data=DT, aes(delta_ps)) + 
      geom_histogram(bins = 10)+
      xlab("p-value")
     write.csv(DT, 'DT.txt')
     write.csv(obs_delta, 'obs_delta')
     DT = read.csv('DT.txt')
     break
  }
  
  #if(obs_delta>0){
  #  res = t.test(delta_ps, mu = obs_delta, alternative='less', conf.level = 0.99)
  #}else if(obs_delta<0){
  #  res = t.test(delta_ps, mu = obs_delta, alternative='greater', conf.level = 0.99)
  #}else{
  #  message('obs_delta is 0')
  #}

  #if(res$p.value==0){
  #  stopifnot(res$conf.int[2] < obs_delta)
  #}
  
  right_tail = quantile(delta_ps, .999)
  left_tail = quantile(delta_ps, .001)
  
  if (obs_delta < left_tail | obs_delta > right_tail){
    cnt = cnt +  1
    cat(dat$chromo[i])
    cat('\t')
    cat(dat$position[i])
    cat('\n')
  }
  #p_value <- length(delta_ps[delta_ps>obs_delta])/length(delta_ps)
  #p_value = res$p.value

  #if(p_value==0){
  #  message('p value is 0')
  #}else if(p_value==1){
  #  message('p value is 1')
  #}
   
  flush.console()
  #cat(p_value)
  #cat('\n')
  #hist(delta_ps, xlab="delta_p", main = paste0('p value: ',p_value))
  
  #abline(v=obs_delta, col='red')
  #abline(v=right_tail, col='pink')
  #abline(v=left_tail, col='blue')
  
}
sink()
message(cnt)

#p_values <- read.delim(outputfile,header=FALSE) #read p values from p_values.txt
#out = data.frame(chromo=dat$chromo, position=dat$position, p_value=p_values)
#write.table(out, file = "p_value_list_all.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
