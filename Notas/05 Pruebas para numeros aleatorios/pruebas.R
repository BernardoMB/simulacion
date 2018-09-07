# 02/09/2018
# Algunas pruebas para números pseudoaleatorios
# Autor: Bernardo Mondragón Brozon

# Prueba de autocorrelación para uniformes
cor.j <- function(x, j) {
  n <- length(x)
  xbar <- mean(x)
  if (j >= n) stop("j debe ser menor que length(x)")
  rhoj<-sum((x[1:(n-j)]-xbar)*(x[(1+j):n]-xbar))/(n-j)/var(x)
  rhoj
}

cor.unif.j<-function(x,j){
  n<-length(x)
  h<-floor((n-1)/j)-1
  v<-NULL
  rhoj<-12*sum(x[1+(0:h)*j]*x[1+((0:h)+1)*j])/(h+1)-3
  rhoj
}

prueba.correl<-function(x,sig=0.01){
  A<-NULL
  n<-length(x)
  for (i in 1:floor(n/4)){
    h<-floor((n-1)/i)-1
    A<-c(A,cor.unif.j(x,i)*(h+1)/sqrt(13*h+7))
  }
  pvals<-pnorm(abs(A),lower.tail=F)
  plot(pvals,xlab="j",ylab="p-vals",ylim=c(0,0.5),main="Prueba de independencia",cex=0.5)
  dep<-match(pvals[pvals<sig/2],pvals)
  text(dep*1.02,pvals[dep]*1.1,labels=dep,cex=0.7)
  abline(h=sig/2,col="red",lwd=3)
  return(A=A,pvals=pvals,j=dep)
}

prueba.rachas <- function(x){
  a <- matrix(c(4529.4,9044.9,13568,18091,22615,27892,
                9044.9,18097,27139,36187,45234,55789,
                13568,27139,40721,54281,67852,83685,
                18091,36187,54281,72414,90470,111580,
                22615,45234,67852,90470,113262,139476,
                27892,55789,83685,111580,139476,172860),nrow=6)
  b <- c(1/6,5/24,11/120,19/720,29/5040,1/840)
  n <- length(x)
  x1 <- c(1)
  x1[2:length(x)] <- sign(diff(x))
  cambios <- c((1:length(x1))[x1==-1],length(x)+1)
  tabla <- table(c(cambios[1]-1,diff(cambios)),exclude=NULL)
  aa <- tabla[match(1:length(x),as.numeric(names(tabla)))]
  aa <- ifelse(is.na(aa),0,aa)
  aa[6] <- sum(aa[6:n])
  r <- aa[1:6]
  names(r) <- c(1:5,">=6")
  R <- as.numeric((r-n*b)%*%a%*%t(t((r-n*b)))/n)
  return(list(x=head(x),R=R,r=r,pval=pchisq(R,6,lower.tail=F)))
}

prueba.chisq.uniforme<-function(x,k=ceiling(length(x)/5)){
  n<-length(x)
  part<-seq(0,1,length=k+1)
  z<-hist(x,breaks=part,plot=F)$counts
  ch<-(k/n)*sum((z-n/k)^2)
  pval<-pchisq(ch,k-1,lower.tail=F)
  return(part=part,freqs=z,estadistica=ch,pval=pval)
}



#Prueba de poker para numeros aleatorios.

separa <- function(x){
  w <- substring(as.character(x),3)
  if(nchar(w) < 4) for(i in 1:(4-nchar(w))) w <- paste(w, "0", sep = "")
  return(as.numeric(unlist(strsplit(w,""))))
}

tabla <- function(x,k){
  tabla <- table(x)
  aa <- tabla[match(0:(k-1),as.numeric(names(tabla)))]
  aa <- ifelse(is.na(aa),0,aa)
  r <- aa[1:k]
  names(r) <- 0:(k-1)
  return(r)  
}

prueba.poker<-function(vec){
  z <- round(vec,4) #redondeo a 4 decimales
  N <- length(z)
  #apaga los warnings por un momento
  ow <- options("warn")
  options(warn = -1)
  dim(z) <- c(N,1) #convierte el vector de datos a una matriz columna
  z1 <- apply(z, 1, separa) #separa cada numero aleatorio en los componentes
  z2 <- t(apply(z1, 2, tabla, k=10)) #crea una tabla de frecuencias
  #    z3<-matrix(unlist(z2),nrow=N,ncol=10,byrow=T) #crea matriz de frecuencias para digitos
  # frecuencias de ceros unos, dos y tres de la tabla de frecuencias anterior,
  # para caracterizar los posibles juegos en cada "mano":
  # pachuca: hay 6 ceros y 4 unos siempre.
  # un par : hay 7 ceros y 1 dos y 2 unos.
  # dos pares: hay 2 dos, 8 ceros
  # una tercia: hay 1 tres y 1 uno, 8 ceros
  # un pokar: hay 1 cuatro y 9 ceros        
  z1 <- apply(z2, 1, table) 
  pachuca <- sum(unlist(lapply(z1,function(x){ifelse(x[1]==6,1,0)})))
  unpar <- sum(unlist(lapply(z1,function(x){ifelse(x[1]==7,1,0)})))
  dospar <- sum(unlist(lapply(z1,function(x){ifelse((x[1]==8)&(x[2]==2),1,0)})))
  tercia <- sum(unlist(lapply(z1,function(x){ifelse((x[1]==8)&(length(x)==3),1,0)})))
  pokar <- sum(unlist(lapply(z1,function(x){ifelse(x[1]==9,1,0)})))
  #obs<-apply(z3,2,sum) #suma por columna
  esp <- N*c(0.504, 0.432, 0.027, 0.036, 0.001)
  obs <- c(pachuca, unpar, dospar, tercia, pokar)
  prueba <- sum((obs-esp)^2/esp)
  pval <- 1 - pchisq(prueba,4)
  #names(obs) <- 0:9
  options(ow)
  return(list(cbind(Esperado=esp,Observado=obs),Estadistica=prueba,pval=round(pval,5)))
}


#paradoja de San Petesburgo
#x<-NULL
#sde<-NULL
#for (i in 1:1000){
#    u<-1+rgeom(10000,0.55)
#    x[i]<-mean(2^u)
#    sde[i]<-sd(2^u)
#}
