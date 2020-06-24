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


res1 = draw_distribution(100, 25, 0.0037)
res2 = draw_distribution(40, 10, 0.0037)
plot( seq(1:(n-1))/n, abs(res1$points-res2$points),  xlab="p", ylab="probability")

X = seq(1:(n-1))/n
Y = res1$points
plot( X, Y, xlab="p", ylab="probability")
X

n=100
plot( seq(1:(n-1))/n, res1$points, xlab="p", ylab="probability")
n = 40
plot( seq(1:(n-1))/n, res2$points, xlab="p", ylab="probability")


n = 100
N1 = draw_distribution(n, 25, 0.0037)

num_sample = 10000
sample_p = sample(seq(1:(n-1))/n, num_sample, prob=N1$points, replace=TRUE)

delta_ps = c()
for(j in seq(1,num_sample)){
  p = sample_p[j]
  delta_p = (0.25-p)
  delta_ps = c(delta_ps, delta_p)
}

hist(delta_ps, xlab="delta_p", main = "Null distribution of deltap from 10000 iterations ")

