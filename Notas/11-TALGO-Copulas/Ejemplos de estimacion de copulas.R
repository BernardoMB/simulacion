
library(copula)

# Ejemplo 1
copula.normal4 <- ellipCopula(family = "normal", dim = 4, dispstr = "un",
                              param = c(0.4,0.5,0.2,0,0.3,0.8))
copula.normal4 #objeto de clase normalCopula

u <- rCopula(200,copula.normal4) #genera observaciones de la cópula construida
cor(u)
pairs(u,pch=16, cex=0.5)

# Ejemplo 2
micopula.t3 <- ellipCopula(family = "t", dim = 3, dispstr = "toep",
                           param = c(0.8,0.5), df = 8)
micopula.t3 #objeto de clase tCopula
rCopula(5,micopula.t3) #genera cinco observaciones de la cópula construida

# Ejemplo 3 (Modelar eventos que no tienen tanto riesgo de manera simultanea. Ej siniestros de auto)
# Aqui los eventos se egeneran de manera simultanea en el inicio mas que en la cola de la distribucion
clayton2 <- archmCopula(family = "clayton", dim = 2, param = 2)
clayton2 #el programa llama alpha al parámetro

# Generemos una muestra de ésta cópula:
y <- rCopula(1000,clayton2)

# Graficando
par(mfrow=c(1,2))
contour(clayton2,dCopula) #gráfica de curvas de nivel
plot(y,cex=0.3)

# Ejemplo 4
frank2 <- archmCopula(family = "frank", dim = 2, param = 2)
frank2 #el programa llama alpha al parámetro

# Generemos una muestra de ésta cópula:
y.4 <- rCopula(1000,frank2)

# Graficando
par(mfrow=c(1,2))
contour(frank2,dCopula) #gráfica de curvas de nivel
plot(y.4,cex=0.3)

# Ejemplo 5
gumbel2 <- archmCopula(family = "gumbel", dim = 2, param = 2)
gumbel2 #el programa llama alpha al parámetro

# Generemos una muestra de ésta cópula:
y.5 <- rCopula(1000,gumbel2)

# Graficando
par(mfrow=c(1,2))
contour(gumbel2,dCopula) #gráfica de curvas de nivel
plot(y.5,cex=0.3)

# Combinando
par(mfrow=c(1,3))
contour(clayton2,dCopula)
contour(frank2,dCopula)
contour(gumbel2,dCopula)

# Otro ejemplo
copula.Frank5 <- archmCopula(family = "frank", dim = 3, param = 5)
micopula <- mvdc(copula = copula.Frank5, margins = c("norm","pois","gamma"),
                 paramMargins = list(list(mean=10,sd=2), list(lambda=3), list(shape=2,scale=4)))
u <- rMvdc(300,micopula) #muestra aleatoria
par(mar=c(1,1,1,1));pairs(u,pch=16,cex=0.5)

# Comparando
library(scatterplot3d)
par(mfrow=c(1,2),mar=c(1,2,1,1),oma=c(0,0,1,1),mgp=c(2,1,0))
u <- rMvdc(200,micopula)
scatterplot3d(u,cex.symbols=0.5,pch=16)
v <- rCopula(200,copula.Frank5)
scatterplot3d(v,cex.symbols=0.5,pch=16)

# Estimacion de coulas
suppressMessages(library(Ecdat)) # fuente de datos
library(copula)
suppressMessages(library(fGarch)) # función de densidad t estandarizada
suppressMessages(library(MASS)) # usa las funciones fitdistr y kde2d
suppressMessages(library(fCopulae)) # funciones adicionales de copula (pempiricalCopula y ellipticalCopulaFit)
data(CRSPday,package="Ecdat")
head(as.data.frame(CRSPday)) # muestra la estructura de los datos

ibm <- CRSPday[,5]; crsp <- CRSPday[,7]
n <- length(ibm); n #número de observaciones

par(pty = "s"); plot(ibm, crsp, cex = 0.4, pch = 16)
abline(h = 0, v = 0, col="red")

# Estimando las marginales
est.ibm <- as.numeric(fitdistr(ibm,"t")$estimate) #parámetros t: media, escala, gl
est.crsp <- as.numeric(fitdistr(crsp,"t")$estimate)
#Convierte los parámetros de escala a desviaciones estándar en el caso de la t
est.ibm[2] <- est.ibm[2]*sqrt(est.ibm[3]/(est.ibm[3]-2))
est.crsp[2] <- est.crsp[2]*sqrt(est.crsp[3]/(est.crsp[3]-2))
#Grados de libertad para cada caso
est.ibm[3]
est.crsp[3]

# Ahora estimo con que copula las voy a pegar
tau <- cor(ibm,crsp,method = "kendall")
omega <- 2/pi*asin(tau)
c(tau,omega)

copula2 <- tCopula(omega,dim=2,dispstr = "un",df = 2)

# Ajuste de la copula especifica
# La función pstd es la distribución estándar t
# por el método de máxima verosimilitud
d1 <- cbind(pstd(ibm, mean = est.ibm[1], sd = est.ibm[2], nu = est.ibm[3]),
            pstd(crsp, mean = est.crsp[1], sd = est.crsp[2], nu = est.crsp[3]))
fit1 <- fitCopula(copula2, method = "ml", optim.method = "L-BFGS-B", data = d1,
                  start = c(omega,5), lower = c(0,2.5), upper = c(0.5,15))
fit1

#Ajusta copula normal
fnorm <- fitCopula(data=d1,copula=normalCopula(-0.3,dim=2),
                   method="ml",optim.method="BFGS",start=0.5)
#Ajusta Gumbel
fgumbel <- fitCopula(data=d1,copula=gumbelCopula(3,dim=2),
                     method="ml",optim.method="BFGS",start=2)
#Ajusta Frank
ffrank <- fitCopula(data=d1,copula=frankCopula(3,dim=2),
                    method="ml",optim.method="BFGS",start=2)
#Ajusta Clayton
fclayton <- fitCopula(data=d1,copula=claytonCopula(3,dim=2),
                      method="ml",optim.method="BFGS",start=2)

# Comparacion grafica
u <- d1
dem <- pempiricalCopula(u[,1],u[,2])
par(mfrow=c(3,2),mar=c(2,2,2,2))
contour(dem$x,dem$y,dem$z,main="Cópula Empírica")
contour(tCopula(fit1@estimate[1],df=round(fit1@estimate[2],0)),pCopula,main="Cópula t")
contour(normalCopula(fnorm@estimate),pCopula,main="Cópula Normal")
contour(gumbelCopula(fgumbel@estimate),pCopula,main="Cópula Gumbel")
contour(frankCopula(ffrank@estimate),pCopula,main="Cópula Frank")
contour(claytonCopula(fclayton@estimate),pCopula,main="Cópula Clayton")
# Todas se parecen a la empirica

# Ahora las densidades bivariadas. Ajustamos la empirica a los datos
par(mfrow=c(3,2),mar=c(2,2,2,2))
contour(kde2d(u[,1],u[,2]),main="KDE")
contour(tCopula(fit1@estimate[1],df=fit1@estimate[2]),dCopula,
        main="Cópula t",nlevels=25)
contour(normalCopula(fnorm@estimate),dCopula,main="Cópula Normal",nlevels=25)
contour(gumbelCopula(fgumbel@estimate),dCopula,main="Cópula Gumbel",nlevels=25)
contour(frankCopula(ffrank@estimate),dCopula,main="Cópula Frank",nlevels=25)
contour(claytonCopula(fclayton@estimate),dCopula,main="Cópula Clayton",nlevels=25)

# Discernimos mediando el criterio de Akaike
#AIC Normal
2*length(fnorm@estimate)-2*fnorm@loglik
#AIC Gumbel
2*length(fgumbel@estimate)-2*fgumbel@loglik
#AIC frank
2*length(ffrank@estimate)-2*ffrank@loglik
#AIC Clayton
2*length(fclayton@estimate)-2*fclayton@loglik
#AIC t
2*length(fit1@estimate)-2*fit1@loglik







