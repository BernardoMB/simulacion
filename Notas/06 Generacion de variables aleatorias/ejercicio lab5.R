#Muestrear la mitad de una distribucion normal
fn<-function(x){2/sqrt(2*pi)*exp(-(x^2)/2)}
g<-function(x,lambda)lambda*exp(-x*lambda)
c<-sqrt(2/(pi*lambda^2))*exp((lambda^2)/2)

#Algoritmo
lambda<-1
x<-numeric()
for (i in 1:1000) {
  u1<-runif(1)
y<-(-lambda*log(u1))
u2<-runif(1)
x<-append(x,if(u2*c*g(y,lambda=lambda) <= fn(y)){y})
}
x
#Recordar que esta es una media normal


#Si quisieramos simular una normal
lambda<-1
x<-numeric()
for (i in 1:1000) {
  u1<-runif(1)
  y<-(-lambda*log(u1))
  u2<-runif(1)
  signo<-ifelse(runif(1)<0.5,-1,1)
  x<-append(x,if(u2*c*g(y,lambda=lambda) <= fn(y)){signo*y})
}
x


qqnorm(x)
hist(x)
ks.test(x,"pnorm")