Delta p plot

files <- list.files('.', pattern = "*.txt")
observation_file = 'obs_deltap.output'

observations = read.delim(observation_file, header = TRUE, sep = "\t", dec = ".")
observations = observations$deltaP

deltaP_matix = c()
for(file in files){
  dat <- read.delim(file, header = TRUE, sep = "\t", dec = ".")
  deltaP_matix = cbind(deltaP_matix,dat$deltaP)
}


dim(deltaP_matix)
mins = c()
maxs = c()
mids = c()
tags = c()

num_snap = 386
for(snap in seq(1,num_snap)){
  delta_Ps = deltaP_matix[snap,]
  q0 = unname(quantile(delta_Ps, probs=0.25))
  mins = c(mins, q0)
  q1 = unname(quantile(delta_Ps, probs=0.99))
  maxs = c(maxs, q1)
  q2 = unname(quantile(delta_Ps, probs=0.5))
  mids = c(mids, q2)
  
  observation = observations[snap]
  if(observation>q1){
    tags = c(tags, 1)
    print(length(tags))
  }
  else{
    tags = c(tags, 0)
  }
  
}




library(ggplot2)
DATA = data.frame(X=seq(1,num_snap), MIN=mins, MAX=maxs, MID=mids,OBS=observations)


sp <- ggplot(DATA, aes(x=X, y=MID)) +
  geom_point(size=.5)+
  geom_point(aes(x=X, y=OBS),size=.5,color='red')+
  #draws the CI error bars
  geom_errorbar(data=DATA, aes(ymin=MIN, ymax=MAX), width=.001,color='yellow',alpha=.2)

sp + scale_x_continuous(name="SNP", limits=c(0, 400)) +
  scale_y_continuous(name="absolute deltap", limits=c(0, 0.5))
