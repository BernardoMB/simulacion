library(MASS)
data(cats)
summary(cats)
library(boot)
q95 <- function(x,ind) return(c(quantile(x[ind],0.95),var(x[ind])/length(ind)))

peso <- cats$Bwt
bootnp <- boot::boot(data=peso, statistic=q95,R=1000)
intconf <- boot.ci(bootnp, type=c("norm","basic","perc","stud"))

valores <- data.frame(normal=intconf$normal[c(2,3)],
                      basico=intconf$basic[c(4,5)],
                      t=intconf$student[c(4,5)],
                      perc=intconf$percent[c(4,5)])

dotchart(as.matrix(valores))

# Modelo de regresion lineal
lm(Hwt~Bwt, data=cats)
alfa1 <- function(x,ind) return(c(lm(x[ind,3]~x[ind,2])$coef[2],var(x[ind,2:3])[2,2]/length(ind)))
bootnp <- boot::boot(data=cats, statistic=alfa1,R=1000)
intconf <- boot.ci(bootnp, type=c("norm","basic","perc","stud"))
valores <- data.frame(normal=intconf$normal[c(2,3)],
                      basico=intconf$basic[c(4,5)],
                      t=intconf$student[c(4,5)],
                      perc=intconf$percent[c(4,5)])
dotchart(as.matrix(valores))
abline(v=4.0341, col="red")
# No se recomienda usar el intervalo basado en elos percentiles

# Ahora con Bootstrap parametrico
