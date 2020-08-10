args = commandArgs(trailingOnly=TRUE)

library(gtools)
library(hash)

options(scipen=999)

IDX = as.numeric(args[1])
message('Thread:')
message(IDX)

IDX = 0
total_n = 1732037
interval = total_n/10
idxs = seq(interval*IDX+1, interval*(IDX+1))
message(idxs[1])
message('-------')
message(idxs[length(idxs)])

filename = 'REF_maf0.05_pctind0.7_cv30.mafs'
challenge = 'CH_maf0.05_pctind0.7_cv30.mafs'
#obs_file = 'obs_deltap_cv30.output'
outputfile = paste0('p_values_theta', IDX, '.txt')
pi_file = 'pi_correct_all.txt'

message(outputfile)

dat2 <- read.delim(challenge, header = TRUE, sep='\t')
n2s <- dat2$nInd*2
n2ks <- round(dat2$knownEM*dat2$nInd*2)
dat <- read.delim(filename, header = TRUE, sep='\t')
#obs_dat <- read.delim(obs_file, header = TRUE, sep='\t') # from ch - ref
dat$n2 = n2s
dat$n2_k = n2ks
pi_dat <- read.delim(pi_file, header = FALSE, sep='\t') # from ch - ref

message(dim(dat)[1])
message(total_n)
#dat = dat[dat$chromo==5,]
#obs_dat= obs_dat[obs_dat$chromo==5,]

draw_distribution <- function(n,k,M){

 
  M = as.numeric(M)

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
  
  # for(i in seq(1:(n/2-1))){
  #   p = i/n
  #   #upper_left = factorial(n)/(factorial(n - k)*factorial(k)) * p^k * (1-p)^(n-k)
  #   upper_left = exp(lfactorial(n)-(lfactorial(n - k)+lfactorial(k))) * p^k * (1-p)^(n-k)
  #   upper_right =  M * (1/i + 1/(n-i))/sum(tajima_points) + M * (1/(n-i) + 1/i)/sum(tajima_points)
  #   res = upper_left * upper_right / (1/(n-1))
  #   points = c(points, res)
  # }
  points=points/sum(points)
  return(list("tajima" = tajima_points, "points" = points))
}


theta_map = hash()

save_theta <- function(pi_dat){
  for(row in 1:nrow(pi_dat)){
    
    if(is.na(pi_dat[row, "V6"])){
      next
    }
    
    cur_chorm = pi_dat[row, "V1"]
    cur_pos = pi_dat[row, "V2"]
    cur_theta = formatC(pi_dat[row, "V6"], digits = 5, format = "f")
    
    theta_map[[paste0(cur_chorm,' ',cur_pos)]] = cur_theta
    
  }
}

save_theta(pi_dat)

null_distribution <- function(n, k, n2, n2_k, theta){
  #theta = 0.0091 #global theta
  N1 = draw_distribution(n, k, theta) 
  
  #plot( seq(1:(n/2-1))/n, N1$points, xlab="p1", ylab="probability")
  #plot( seq(1:(n-1))/n, N11$points, xlab="p11", ylab="probability")

  # uncomment below if the expected allele count is drew from the N1 (i.e. best ref pop)
  n2_k = sample(seq(1:(n-1))/n, 1, prob=N1$points, replace=FALSE) * n2
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
  
  hist(delta_ps)
  #return(length(delta_ps[delta_ps>obs_delta])/num_sample)
  return(list('p1'=sample_p1, 'p2'=sample_p2, 'delta_p'=delta_ps))
}

# main 
cnt = 0
del_cnt = 0
dic <- hash()
p_values = c()
sink(outputfile)
idxs = seq(interval*IDX+1, interval*(IDX+1))
cat(idxs[1])
cat('-------')
cat(idxs[length(idxs)])
cat('\n')

for(i in idxs){
  s = paste0(i,'/',dim(dat)[1])
  message(s,"\r",appendLF=FALSE)
  n=dat$nInd[i]*2
  k=floor(dat$nInd[i]*2*dat$knownEM[i])
  n2=dat$n2[i]
  n2_k=dat$n2_k[i]

  cur_pos = dat$position[i] - dat$position[i]%%200
  theta_key = paste0(dat$chromo[i],' ',cur_pos)
  
  
  if(!has.key(theta_key, theta_map)){
    del_cnt = del_cnt + 1
    next
    #message(theta_key)
    #message(dat$position[i])
  }
  
  theta = theta_map[[theta_key]]
  if(theta<0.000001){
    next
  }
  #obs_delta=obs_dat$deltaP[i]
  obs_delta=dat2$knownEM[i] - dat$knownEM[i]

  key = paste0(n, ' ', k, ' ', n2, ' ', n2_k, ' ', theta)
  #if(has.key(key, dic)){
  #  delta_ps <- dic[[key]]
  #}
  #else{
  #  delta_ps = null_distribution(n=n, k=k, n2=n2, n2_k=n2_k, theta=theta)
  #  dic[[key]] <- delta_ps
  #}
  res = null_distribution(n=n, k=k, n2=n2, n2_k=n2_k, theta=theta)
  delta_ps = res$delta_p
  
  hist(delta_ps, breaks=30, xlab="delta_p")
  exit()
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
 
  #df = data.frame(p1=res$p1, p2=res$p2, delta_p=res$delta_p)
  #write.csv(df,paste0('storeP/',dat$chromo[i],'_',dat$position[i],'.csv'))
  if (obs_delta < left_tail | obs_delta > right_tail){
    # store p
    df = data.frame(p1=res$p1, p2=res$p2, delta_p=res$delta_p)
    write.csv(df,paste0('storeP/',dat$chromo[i],'_',dat$position[i],'.csv'))  
    
    cnt = cnt +  1
    cat(dat$chromo[i])
    cat('\t')
    cat(dat$position[i])
    cat('\n')
  }
  
  
  #p_value <- length(delta_ps[delta_ps>obs_delta])/length(delta_ps)
  
  flush.console()
  #cat(p_value)
  #cat('\n')
  
  
  #abline(v=obs_delta, col='red')
  #abline(v=right_tail, col='pink')
  #abline(v=left_tail, col='blue')
  
}
sink()
message(cnt)
message(del_cnt)


#p_values <- read.delim(outputfile,header=FALSE) #read p values from p_values.txt
#out = data.frame(chromo=dat$chromo, position=dat$position, p_value=p_values)
#write.table(out, file = "p_value_list_all.txt", sep = "\t", row.names = FALSE, col.names = FALSE)


