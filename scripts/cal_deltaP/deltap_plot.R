delta_ps = c()
for(j in seq(1,num_sample)){
  p = sample_p[j]
  delta_p = (0.25-p)
  delta_ps = c(delta_ps, delta_p)
}

hist(delta_ps, xlab="delta_p", main = "Null distribution of deltap from 10000 iterations ")
