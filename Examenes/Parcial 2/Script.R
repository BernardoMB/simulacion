# Pregunta 1 
lambda <-  1/2 # Parametro de la doble exponencial
n <- 10000 # Tamaño de muestra
breaks <- 50 # Prticion del histograma
simulaDobleExponencial <- function(n, lambda) {
  f <- function(x) {(lambda/2)*exp(-lambda*abs(x))}
  x <- NULL
  x0 <- 3
  for(i in 0:n) {
    w <- ifelse(i==0,x0,x[i])
    y <- rnorm(1, mean = w, sd = 2.5)
    alfa <- min(1,(f(y)*dnorm(w,mean=y,sd=1))/(f(w)*dnorm(y,mean=w,sd=1)))
    x <- append(x,ifelse(runif(1)<alfa,y,w))
  }
  return(list(x=x,f=f(sort(x))))
}
data <- simulaDobleExponencial(n, lambda)

par(mfrow = c(3,1))
plot(data$x, type="l", main="Trayectoria del proceso", xlab="t", ylab="X(t)")
hist(data$x, probability=T, breaks = breaks, main="Histograma", xlab="x", ylab="Frec")
lines(sort(data$x), data$f, col="green")

li <- 1
ls <- n

# CDF
plot(ecdf(data$x[li:ls]))
x <- seq(-15,15,by=0.01)
library(rmutil)
lines(x, plaplace(x, m=0, s=1/lambda), col="green")

# One-sample Kolmogorov-Smirnov test
ks.test(data$x[li:ls], "plaplace", 0, 1/lambda)

# Pregunta 2
n <- 10000 # tamaño de la muestra
z <- 3

# Metodo 1
vars <- rexp(n, rate = 1)
vals <- 1-exp((-(z/(vars^2))))
mean(vals)

# Metodo 2
total <- 0
for (i in 1:n) {
  v <- rexp(1, rate = 1)
  w <- rexp(1, rate = 1/v)
  vw <- v*w
  if (vw <= z) {
    total <- total + 1
  }
}
total/n


# Pregunta 3
# Si la distribucion objetivo es la distribucion limite 
# de una caminata aleatoria, entonces se requiere que 

# Pregunta 4

# Pregunta 6
simulaPoisson <- function(n, lambda = 3) {
  x <- NULL
  x0 <- 3
  for(i in 0:n) {
    w <- ifelse(i==0,x0,x[i])
    y <- rgeom(1, prob=1/3)
    #alfa <- min(1,(dpois(y, lambda)*dgeom(w, prob=1/3))/(dpois(w, lambda)*dgeom(y, prob=1/3)))
    alfa <- (dpois(y, lambda)*dgeom(w, prob=1/3))/(dpois(w, lambda)*dgeom(y, prob=1/3))
    x <- append(x,ifelse(runif(1)<alfa,y,w))
  }
  return(list(x=x,f=dpois(sort(x), lambda)))
}
pois <- simulaPoisson(500)
par(mfrow = c(2,1))
plot(pois$x,type="l", main="Trayectoria del proceso", xlab="t", ylab="X(t)")
p1_pois <- hist(pois$x,probability=T, main="Histograma", xlab="x", ylab="Frec")
lines(sort(pois$x),pois$f,col="red")

# Prueba de bondad de ajuste Kolmogorov-Smirnov
y <- rpois(500, 3)
#install.packages("stats")
library(stats)
ks_pois <- ks.test(pois$x, y)





