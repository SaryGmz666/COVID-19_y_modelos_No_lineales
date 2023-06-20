library(ggplot2)
library(minpack.lm)
library(ggpubr)
library(MASS)
suppressPackageStartupMessages(library(tidyverse))
library(tidyverse)
library(ISLR)
library(MLmetrics)
library(survival)
library(fitdistrplus)


####################    Datos    ####################
# Lectura de incidencia acumulada de casos
Diasiniciopand <- 197
ff   <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master"
ff   <- paste(ff,"/csse_covid_19_data/csse_covid_19_time_series/",sep="")
ff   <- paste(ff,"time_series_covid19_confirmed_global.csv",sep="")
dat  <- read.csv2(ff,header=TRUE,sep=",")
dim(dat) # 266 x 151
nc   <- dim(dat)[2] 

# Paises en la base de datos
pais <- unique(dat[,2])
np   <- length(pais)

# Colapsar provincias y tener un sólo registro por país
# (primeras 4 columnas de dat contienen: Provincia, País, Lat, Lon)
datu <- matrix(0,np,nc-4) 
for(i in 1:np){ datu[i,] <- colSums(dat[dat[,2]==pais[i],-(1:4)]) }

# Período de tiempo en la base de datos
pd   <- names(dat)[5]   # primer día
nch  <- nchar(pd)   
pdia <- as.Date( substr(pd,2,nch), format= "%m.%d.%y")
ud   <- names(dat)[nc]  # último día
nch  <- nchar(ud)                                         
udia <- as.Date( substr(ud,2,nch), format= "%m.%d.%y")    
dd   <- seq.Date(from=pdia,to=udia,by=1) # secuencia de días (no hay días sin reporte)
nd   <- length(dd)                                        

id      <- which(pais=="Mexico")
pob2018 <- 128932753             # Población actual
hab     <- pob2018
unit     <- pais[id]
zz       <- datu[id,]
zzc      <- diff(zz)
nonn     <- min( which(zzc > 0) )
tiempoA  <- dd[-(1:nonn)]
casosA   <- zzc[-(1:(nonn-1))]    
casosaA  <- datu[id,-(1:nonn)]
N        <- length(tiempoA)

grafejesfijos <- list()
i             <- 0
for(MM in Diasiniciopand:N){
  i       <- i+1
  tiempo  <- tiempoA[1:MM]
  casos   <- casosA[1:MM]     
  casosa  <- casosaA[1:MM] 
  dat2    <- data.frame(Unidad=unit, dateRep=tiempo, casos=casos, 
                        casosa=casosa, pob=pob2018)
  #bb      <- Gomp1(dat2)   # Gomp1 es la misma
  h       <- 300 - i
  aux     <- Gomp2(h, coe=bb$Coeficientes, sigma=bb$sigma, dfdato=bb$DatosId)
  pred    <- aux$preds
  VVar    <- aux$Mvarianza
  GHH     <- 360
  ggra    <- Gomp3X(h, pred, coe=bb$Coeficientes, dfdato=bb$DatosId, VV=VVar, GHH=GHH)
  grafejesfijos[[i]] <- ggarrange(ggra$graf1,ggra$graf3, ncol = 1, nrow = 2)
}

#Columna time para simplificacion del modelo
nn       <- length(tiempo)
tt      <- 1:nn
Datos   <- cbind(dat2, Time = tt)
fechau <- tail(Datos$dateRep,1)
paisn       <- (Datos$Unidad)[1]
casosu <- tail(Datos$casosa,1)
casosf <- format( casosu, big.mark = ",")

####################    Modelo Gompetz Acumulado    ####################

#Valores iniciales para las estimaciones

coe1  <-  10 
coe2  <-  5
coe3  <-  0.1

esti.sp  <- nlsLM(casosa ~ coe1*exp(-coe2*exp(-coe3*Time)), start=list(coe1=coe1,coe2=coe2,coe3=coe3),data=Datos)
summary(esti.sp)

alpha   <- coefficients(esti.sp)[1]
beta    <- coefficients(esti.sp)[2]
kappa   <- coefficients(esti.sp)[3]

Gmpz    <-function(x)
{
 alpha*exp(-beta*exp(-kappa*x))
}

##### Datos interesantes del modelo Gompertz #####

# Cuantas personas infectadas habra al final de la pandemia?
CMax  <- round(alpha)

# Cuando va a ser el pico?
PMaxsr  <- log(beta)/kappa
PMax    <- round(PMaxsr)

#Que fecha?
FMax  <- tiempo[1] + PMax

#cuantos casos se esperaban en tal fecha?
CPMax <- alpha*exp(-beta*exp(-kappa*PMax))
CaMaxi<- round(CPMax)

#En que fecha llegara al 95% del maximo total?
cota        <- 0.95
tDuracion   <- round(-log(-log(cota)/beta)/kappa)
fDuracion   <- tiempo[1] + tDuracion

#AcuGmpz   <- ggplot(Datos, aes(tt, casosa)) + geom_point() + xlab("Tiempo") + ylab("Numero de casos acumulados") + 
#             ggtitle(paste("Modelo Gompertz,",paisn,":",fechau,"Casos acumulados = ",casosf," \n ( Total Predichos ",CMax,")")) 
#AcuGmpz   <- AcuGmpz + stat_function(fun = Gmpz, colour = "red")

AcuGmpz2  <- ggplot(Datos) + geom_step(mapping=aes(x=tt,y=casosa), linetype=1,size=.7) + xlab("Tiempo") + 
             ylab("Casos") + ggtitle(paste("Modelo Gompertz,",paisn,":",fechau,"Casos acumulados = ",casosf," 
             ( Total Predichos ",CMax,")")) + stat_function(fun = Gmpz, colour = "red", linetype=1,size=.9)
AcuGmpz2

##### Calclulo de los errores del modelo y prueba de bondad de ajuste #####

#Errores del modelo acumulado
ErrorGmpzA   <- replicate(Diasiniciopand, 0)
ValorModA    <- replicate(Diasiniciopand, 0)
ValoresA   <- Datos$casosa

for (i in 0:Diasiniciopand) 
{
  ValorModA[i] <- alpha*exp(-beta*exp(-kappa*i))
}
ErrorGmpzA   <- ValorModA - ValoresA
Datos        <- cbind(Datos, ErrorG = ErrorGmpzA)

#Prueba de bondad de ajuste
hist(Datos$ErrorG, main = "Histograma", xlab = "Error Modelo Gompertz")

ajuste <- fitdist(Datos$ErrorG, "norm")
ajuste$estimate

plot(ajuste)
prueba <- gofstat(ajuste)
prueba$kstest
prueba$chisqpvalue

# MSE
MSEGA = 0
  for (i in 1:length(ErrorGmpzA)) 
  {
    MSEGA = MSEGA + (ErrorGmpzA[i])^2
  }
MSEGA     <- MSEGA/Diasiniciopand
MSEGA

# RMSE
RMSEGA      <- sqrt(MSEGA)
RMSEGA

# MAE
MAEGA = 0
for (i in 1:length(ErrorGmpzA)) 
{
  MAEGA = MAEGA + abs(ErrorGmpzA[i])
}
MAEGA     <- MAEGA/Diasiniciopand
MAEGA

####################    Modelo Gompetz    ####################

GmpzF    <-function(x)
{
  alpha*exp(-beta*exp(-kappa*x))*beta*kappa*exp(-kappa*x)
}

yendp      <- alpha*exp(-beta*exp(-kappa*PMax))*beta*kappa*exp(-kappa*PMax)
yend9      <- alpha*exp(-beta*exp(-kappa*tDuracion))*beta*kappa*exp(-kappa*tDuracion)

GmpzFull  <- ggplot(Datos) + geom_step(mapping=aes(x=tt,y=casos), linetype=1,size=.5) + xlab("Tiempo") + 
             ylab("Casos") + ggtitle(paste("Modelo Gompertz: Casos diarios COVID en" ,paisn, "\nCasos Acum Pico" ,CaMaxi))
GmpzFull  <- GmpzFull + stat_function(fun = GmpzF, colour = "red", linetype=1,size=1) + geom_text(x=PMax,y=0,label=FMax) + 
             geom_text(x=tDuracion,y=0,label=fDuracion) + geom_segment(aes(x= PMax, y=0, xend= PMax, yend= yendp), colour= "orange",
             linetype=1,size=1) + geom_segment(aes(x= tDuracion, y=0, xend= tDuracion, yend= yend9), colour= "orange", linetype=1,size=1)
GmpzFull

#GmpzFull2  <- ggplot(Datos) + geom_linerange(mapping=aes(x=tt,y=casos), ymin=0, ymax=casos, linetype=1,size=.5) + xlab("Tiempo") +
#              ylab("Casos diarios") + ggtitle(paste("Modelo Gompertz: Casos diarios COVID en" ,paisn, "\nCasos Acum Pico" ,CaMaxi))
#GmpzFull2  <- GmpzFull2 + stat_function(fun = GmpzF, colour = "red", linetype=1,size=1) + geom_text(x=PMax,y=0,label=FMax) + 
#              geom_text(x=tDuracion,y=0,label=fDuracion) + geom_segment(aes(x= PMax, y=0, xend= PMax, yend= yendp), colour= "orange",
#              linetype=1,size=1) + geom_segment(aes(x= tDuracion, y=0, xend= tDuracion, yend= yend9), colour= "orange", linetype=1,size=1)
#GmpzFull2


##### Calclulo de los errores del modelo #####

ErrorGmpzF   <- replicate(Diasiniciopand, 0)
ValorModF    <- replicate(Diasiniciopand, 0)
ValoresF   <- Datos$casos

for (i in 0:Diasiniciopand) 
{
  ValorModF[i] <- alpha*exp(-beta*exp(-kappa*i))*beta*kappa*exp(-kappa*i)
}

ErrorGmpzF   <- ValorModF - ValoresF

# MSE
MSEGF = 0
for (i in 1:length(ErrorGmpzF)) 
{
  MSEGF = MSEGF + (ErrorGmpzF[i])^2
}
MSEGF     <- MSEGF/Diasiniciopand
MSEGF

# RMSE
RMSEGF      <- sqrt(MSEGF)
RMSEGF

# MAE
MAEGF = 0
for (i in 1:length(ErrorGmpzF)) 
{
  MAEGF = MAEGF + abs(ErrorGmpzF[i])
}
MAEGF     <- MAEGF/Diasiniciopand
MAEGF

####################    Modelo Logistico Acumulado    ####################

#Valores iniciales para las estimaciones

gam0  <-  100
gam1  <-  10
gam2  <-  0.2

esti.spl  <- nlsLM(casosa ~ gam0/(1+gam1*exp(-gam2*Time)), start=list(gam0=gam0,gam1=gam1,gam2=gam2),data=Datos)
summary(esti.spl)

gamma0   <- coefficients(esti.spl)[1]
gamma1   <- coefficients(esti.spl)[2]
gamma2   <- coefficients(esti.spl)[3]

Log    <-function(x)
{
  gamma0/(1+gamma1*exp(-gamma2*x))
}

##### Datos interesantes del modelo Logistico #####

# Cuantas personas infectadas habra al final de la pandemia?
CMaxL  <- round(gamma0)

# Cuando va a ser el pico?
PMaxLsr  <- log(gamma1)/gamma2
PMaxL    <- round(PMaxLsr)

#Que fecha?
FMaxL  <- tiempo[1] + PMaxL

#cuantos casos se esperaban en tal fecha?
CPMaxL  <- gamma0/(1+gamma1*exp(-gamma2*PMaxL))
CaMaxiL <- round(CPMaxL)

#En que fecha llegara al 95% del maximo PMaxL?
tDuracionL   <- round(-log(((1/cota) - 1)/gamma1)/gamma2)     
fDuracionL   <- tiempo[1] + tDuracionL

#AcuLog   <- ggplot(Datos, aes(tt, casosa)) + geom_point() + xlab("Tiempo") + ylab("Numero de casos acumulados") + 
#            ggtitle(paste("Modelo Logistico,",paisn,":",fechau,"Casos acumulados = ",casosf," \n ( Total Predichos ",CMaxL,")")) 
#AcuLog   <- AcuLog + stat_function(fun = Gmpz, colour = "blue")
#AcuLog

AcuLog2  <- ggplot(Datos) + geom_step(mapping=aes(x=tt,y=casosa), linetype=1,size=.7,direction="hv") + xlab("Tiempo") + 
            ylab("Casos") + ggtitle(paste("Modelo Logistico,",paisn,":",fechau,"Casos acumulados = ",casosf," 
            ( Total Predichos ",CMaxL,")")) + stat_function(fun = Gmpz, colour = "blue", linetype=1,size=.9)
AcuLog2

##### Calclulo de los errores del modelo y prueba de bondad de ajuste #####

#Errores del modelo acumulado
ErrorLogA   <- replicate(Diasiniciopand, 0)
ValorModLA  <- replicate(Diasiniciopand, 0)

for (i in 0:Diasiniciopand) 
{
  ValorModLA[i] <- gamma0/(1+gamma1*exp(-gamma2*i))
}
ErrorLogA   <- ValorModLA - ValoresA
Datos        <- cbind(Datos, ErrorL = ErrorLogA)

#Prueba de bondad de ajuste
hist(Datos$ErrorG, main = "Histograma", xlab = "Error Modelo Logistico")

ajuste <- fitdist(Datos$ErrorL, "norm")
ajuste$estimate

plot(ajuste)
prueba <- gofstat(ajuste)
prueba$kstest
prueba$chisqpvalue

# MSE
MSELA = 0
for (i in 1:length(ErrorLogA)) 
{
  MSELA = MSELA + (ErrorLogA[i])^2
}
MSELA     <- MSELA/Diasiniciopand
MSELA

# RMSE
RMSELA    <- sqrt(MSELA)
RMSELA

# MAE
MAELA = 0
for (i in 1:length(ErrorLogA)) 
{
  MAELA = MAELA + abs(ErrorLogA[i])
}
MAELA     <- MAELA/Diasiniciopand
MAELA

####################    Modelo Logistico    ####################

LogF    <-function(x)
{
  gamma0*gamma1*gamma2*exp(-gamma2*x)/((1+gamma1*exp(-gamma2*x))^2)
}

yendpL      <- gamma0*gamma1*gamma2*exp(-gamma2*PMaxL)/((1+gamma1*exp(-gamma2*PMaxL))^2)
yend9L      <- gamma0*gamma1*gamma2*exp(-gamma2*tDuracionL)/((1+gamma1*exp(-gamma2*tDuracionL))^2)

LogFull  <- ggplot(Datos) + geom_step(mapping=aes(x=tt,y=casos), linetype=1,size=.5,direction="hv") + xlab("Tiempo") + 
            ylab("Casos") + ggtitle(paste("Modelo Logistico: Casos diarios COVID en" ,paisn, "\nCasos Acum Pico" ,CaMaxiL))
LogFull  <- LogFull + stat_function(fun = LogF, colour = "blue", linetype=1,size=1) + geom_text(x=PMaxL,y=0,label=FMaxL) + 
            geom_text(x=tDuracionL,y=0,label=fDuracionL) + geom_segment(aes(x= PMaxL, y=0, xend= PMaxL, yend= yendpL), colour= "hotpink4",
            linetype=1,size=1) + geom_segment(aes(x= tDuracionL, y=0, xend= tDuracionL, yend= yend9L), colour= "hotpink4", linetype=1,size=1)
LogFull

#LogFull2  <- ggplot(Datos) + geom_linerange(mapping=aes(x=tt,y=casos), ymin=0, ymax=casos, linetype=1,size=.5) + xlab("Tiempo") +
#             ylab("Casos diarios") + ggtitle(paste("Modelo Gompertz: Casos diarios COVID en" ,paisn, "\nCasos Acum Pico" ,CaMaxi))
#LogFull2  <- LogFull2 + stat_function(fun = LogF, colour = "blue", linetype=1,size=1) + geom_text(x=PMaxL,y=0,label=FMaxL) + 
#             geom_text(x=tDuracionL,y=0,label=fDuracionL) + geom_segment(aes(x= PMaxL, y=0, xend= PMaxL, yend= yendpL), colour= "hotpink4",
#             linetype=1,size=1) + geom_segment(aes(x= tDuracionL, y=0, xend= tDuracionL, yend= yend9L), colour= "hotpink4", linetype=1,size=1)

##### Calclulo de los errores del modelo #####

ErrorLogF   <- replicate(Diasiniciopand, 0)
ValorLogF    <- replicate(Diasiniciopand, 0)

for (i in 0:Diasiniciopand) 
{
  ValorLogF[i] <- gamma0*gamma1*gamma2*exp(-gamma2*i)/((1+gamma1*exp(-gamma2*i))^2)
}

ErrorLogF   <- ValorLogF - ValoresF

# MSE
MSELF = 0
for (i in 1:length(ErrorLogF)) 
{
  MSELF = MSELF + (ErrorLogF[i])^2
}
MSELF     <- MSELF/Diasiniciopand
MSELF

# RMSE
RMSELF    <- sqrt(MSELF)
RMSELF

# MAE
MAELF = 0
for (i in 1:length(ErrorLogF)) 
{
  MAELF = MAELF + abs(ErrorLogF[i])
}
MAELF     <- MAELF/Diasiniciopand
MAELF

####################    Modelo Normal Acumulado    ####################

#Valores iniciales para las estimaciones

a <-  1
b <-  50
p <-  10000

esti.spn  <- nlsLM(casosa ~ (p/2)*pnorm(a*(Time-b)), start=list(a=a,b=b, p=p),data=Datos)
summary(esti.spn)

af  <- coefficients(esti.spn)[1]
bf  <- coefficients(esti.spn)[2]
pf  <- coefficients(esti.spn)[3]

Nom    <-function(x)
{
  (pf/2)*pnorm(af*(x-bf))
}

##### Datos interesantes del modelo Normal #####

# Cuantas personas infectadas habra al final de la pandemia?
CMaxN  <- round(pf/2)

# Cuando va a ser el pico?
PMaxNsr  <- bf
PMaxN    <- round(PMaxNsr)

#Que fecha?
FMaxN  <- tiempo[1] + PMaxN

#cuantos casos se esperaban en tal fecha?
CPMaxN  <- (pf/2)*pnorm(af*(PMaxN-bf))
CaMaxiN <- round(CPMaxN)

#En que fecha llegara al 95% del maximo PMaxL?
tDuracionN   <- round(bf+ qnorm(.95)/af)
fDuracionN   <- tiempo[1] + tDuracionN

#AcuNom   <- ggplot(Datos, aes(tt, casosa)) + geom_point() + xlab("Tiempo") + ylab("Numero de casos acumulados") + 
#            ggtitle(paste("Modelo Normal,",paisn,":",fechau,"Casos acumulados = ",casosf," \n ( Total Predichos ",CMaxN,")")) 
#AcuNom   <- AcuNom + stat_function(fun = Nom, colour = "gold")
#AcuNom

AcuNom2  <- ggplot(Datos) + geom_step(mapping=aes(x=tt,y=casosa), linetype=1,size=.7,) + xlab("Tiempo") + 
            ylab("Casos") + ggtitle(paste("Modelo Normal,",paisn,":",fechau,"Casos acumulados = ",casosf," 
            ( Total Predichos ",CMaxN,")")) + stat_function(fun = Nom, colour = "gold", linetype=1,size=.9)
AcuNom2

##### Calclulo de los errores del modelo y prueba de bondad de ajuste #####

#Errores del modelo acumulado
ErrorNomA   <- replicate(Diasiniciopand, 0)
ValorModNA  <- replicate(Diasiniciopand, 0)

for (i in 0:Diasiniciopand) 
{
  ValorModNA[i] <- (pf/2)*pnorm(af*(i-bf))
}
ErrorNomA   <- ValorModNA - ValoresA
Datos        <- cbind(Datos, ErrorN = ErrorNomA)

#Prueba de bondad de ajuste
hist(Datos$ErrorG, main = "Histograma", xlab = "Error Modelo Normal")

ajuste <- fitdist(Datos$ErrorN, "norm")
ajuste$estimate

plot(ajuste)
prueba <- gofstat(ajuste)
prueba$kstest
prueba$chisqpvalue

# MSE
MSENA = 0
for (i in 1:length(ErrorNomA)) 
{
  MSENA = MSENA + (ErrorNomA[i])^2
}
MSENA     <- MSENA/Diasiniciopand
MSENA

# RMSE
RMSENA      <- sqrt(MSENA)
RMSENA

# MAE
MAENA = 0
for (i in 1:length(ErrorNomA)) 
{
  MAENA = MAENA + abs(ErrorNomA[i])
}
MAENA     <- MAENA/Diasiniciopand
MAENA

####################    Modelo Normal    ####################

#Valores iniciales para las estimaciones

NomF    <-function(x)
{
  (pf*af/(2*sqrt(2*pi))*exp(-0.5*af^2*(x-bf)^2))
}

yendpN      <- (pf*af/(2*sqrt(2*pi))*exp(-0.5*af^2*(PMaxN-bf)^2))
yend9N      <- (pf*af/(2*sqrt(2*pi))*exp(-0.5*af^2*(tDuracionN-bf)^2))

FullNom  <- ggplot(Datos) + geom_step(mapping=aes(x=tt,y=casos), linetype=1,size=.5,direction="hv") + xlab("Tiempo") + 
            ylab("Casos") + ggtitle(paste("Modelo Normal: Casos diarios COVID en" ,paisn, "\nCasos Acum Pico" ,CaMaxiN))
FullNom  <- FullNom + stat_function(fun = NomF, colour = "gold", linetype=1,size=1) + geom_text(x=PMaxN,y=0,label=FMaxN) + 
            geom_text(x=tDuracionN,y=0,label=fDuracionN) + geom_segment(aes(x= PMaxN, y=0, xend= PMaxN, yend= yendpN), colour= "orange4",
            linetype=1,size=1) + geom_segment(aes(x= tDuracionN, y=0, xend= tDuracionN, yend= yend9N), colour= "orange4", linetype=1,size=1)
FullNom

##### Calclulo de los errores del modelo #####

#Errores del modelo acumulado
ErrorNomF   <- replicate(Diasiniciopand, 0)
ValorModNF  <- replicate(Diasiniciopand, 0)

for (i in 0:Diasiniciopand) 
{
  ValorModNF[i] <- (pf*af/(2*sqrt(2*pi))*exp(-0.5*af^2*(i-bf)^2))
}
ErrorNomF   <- ValorModNF - ValoresA

# MSE
MSENF = 0
for (i in 1:length(ErrorNomF)) 
{
  MSENF = MSENF + (ErrorNomF[i])^2
}
MSENF     <- MSENF/Diasiniciopand
MSENF

# RMSE
RMSENF      <- sqrt(MSENF)
RMSENF

# MAE
MAENF = 0
for (i in 1:length(ErrorNomF)) 
{
  MAENF = MAENF + abs(ErrorNomF[i])
}
MAENF     <- MAENF/Diasiniciopand
MAENF

