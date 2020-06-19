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


res1 = draw_distribution(n, k, 0.0037)

plot( seq(1:(n-1))/n, res1$points-res2$points,  xlab="p", ylab="probability")
plot( seq(1:(n-1))/n, res1$points, xlab="p", ylab="probability")
plot( seq(1:(n-1))/n, res2$points, xlab="p", ylab="probability")

plot( seq(1:(n-1))/n, res1$tajima, xlab="p", ylab="probability")

