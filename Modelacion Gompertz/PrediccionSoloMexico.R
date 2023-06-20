# Modelo Gompertz

library(ggplot2)
library(minpack.lm)
library(ggpubr)
library(MASS)
suppressPackageStartupMessages(library(tidyverse))

# https://community.rstudio.com/c/tidyverse
# https://r4ds.had.co.nz/

## Nota: En líneas 11 - 222 están las funciones Gomp1, Gomp2 y Gomp3.

########################### EMPIEZA ESTIMACIÓN MODELO ##################
# Función para el ajuste de modelo Gompertz
# Argumento dat2:
# dat2    <- data.frame(Unidad=unit, dateRep=tiempo, casos=casos, 
#                      casosa=casosa, pob=pob2018)
# dat2: data frame con 5 columnas. (pob2018 es población actual)

Gomp1 <- function(dat2){
  nacion <- (dat2$Unidad)[1]
  hab    <- ((dat2$pob)[1])/100000  # convertir a pob en ciento de miles
  tiempo <- dat2$dateRep
  casosd <- dat2$casos
  casosa <- dat2$casosa
  nn     <- length(tiempo)
  tt     <- 1:nn
  dfdat  <- data.frame(hab=hab,tt=tt,casosa=casosa)
  dfdato <- data.frame(Unidad=nacion,hab=hab,tiempo=tiempo,tt=tt,casosa=casosa)
  alfa   <- 10                    # Valores iniciales que han funcionado para
  beta   <- 5                     # casi todos los países. Tal vez requieran
  kapa   <- 0.1                   # "tuning" en caso de no convergencia
  outg   <- nlsLM(casosa ~ hab*alfa*exp(-beta*exp(-kapa*tt)),
                  start=list(alfa=alfa,beta=beta,kapa=kapa),
                  control=list(maxiter=1000),data=dfdat)
  aux    <- summary(outg)
  estd   <- (aux$coefficients)[,2]
  sigma  <- aux$sigma
  coe    <- (aux$coefficients)[,1] 
  alfa   <- coe[1]
  beta   <- coe[2]
  kapa   <- coe[3]
  # Cuántas personas infectadas (confirmadas), habrá al final de la epidemia?
  Nmax   <- round(alfa*hab)
  names(Nmax) <- NULL
  dfdato$CasosTotales <- Nmax
  # En que día se estima el pico de la epidemia?
  tmax   <- log(beta)/kapa
  names(tmax) <- NULL
  # En que fecha?
  fmax   <- tiempo[1] + round(tmax)
  names(fmax) <- NULL
  # Cuál será la incidencia de casos en el día del máximo (tmax)?
  aux    <- alfa*exp(-beta*exp(-kapa*tmax))
  ctmax  <- round( hab*beta*kapa*exp(-kapa*tmax)*aux )
  names(ctmax) <- NULL
  # Cuántos casos acumulados habrá hasta el pico? (aux previa)
  casosAP <- round( hab*alfa*exp(-beta*exp(-kapa*tmax)) )
  names(casosAP) <- NULL
  dfdato$CasosAcumPico <- casosAP
  # Cuántos casos ocurrirán entre el pico y una semana antes del pico?
  casos1S <- hab*alfa*( exp(-beta*exp(-kapa*tmax))-exp(-beta*exp(-kapa*(tmax-7))) )
  casos1S <- round(casos1S)  
  names(casos1S) <- NULL
  dfdato$CasosAcum1SyPico <- casos1S
  # Cuánto durará la epidemia? : == Cuándo se llega al 95% del máximo total? 
  cota   <- 0.95
  tdura  <- round(-log(-log(cota)/beta)/kapa)
  names(tdura) <- NULL
  fdura  <- tiempo[1] + tdura
  names(fdura) <- NULL
  # Tasa de casos confirmados por millón de hab
  casosPM <- round( tail(casosa,1)/(hab/10), 1 )
  dfdato$CasosPorMillon <- casosPM
  return(list(Unidad=nacion,
              Coeficientes=coe, Err_Std=estd, sigma=sigma, Casos_Totales=Nmax,
              Dia_Pico=round(tmax), Fecha_Pico=fmax, Casos_Fecha_Pico=ctmax,
              Casos_Acum_al_Pico=casosAP, Casos_Entre_1_Sem_antes_y_el_Pico=casos1S,
              Duracion_Epidemia=tdura, Fecha_Terminacion_Epidemia=fdura,
              Casos_Por_Millon=casosPM,
              DatosId = dfdato) ) }


########################### TERMINA ESTIMACIÓN MODELO ##################



########################### EMPIEZA PREDICCIÓN ##################
# Función para el cálculo de prediciones y límites de confianza
# Argumentos: h (número de días a predecir)
#             coe,sigma,dfdato (salida de Gomp1)
#             conf=0.99 (default)

Gomp2 <- function(h,coe,sigma,dfdato,conf=0.99){
  alfa   <- coe[1]
  beta   <- coe[2]
  kapa   <- coe[3]
  hab    <- (dfdato$hab)[1]
  tt     <- dfdato$tt
  nn     <- length(tt)
  x0     <- c(tt,tt[nn]+(1:h))
  predic <- hab*alfa*exp(-beta*exp(-kapa*x0))
  ff1    <- hab*alfa*exp(-beta*exp(-kapa*tt))
  Df     <- cbind(1/(hab*alfa),-exp(-kapa*tt),beta*tt*exp(-kapa*tt)) * ff1
  aux    <- t(Df)%*%Df
  uu     <- diag(c(1,1,1e-5))
  uu     <- diag(c(1,1e-4,1e-8)) # 24 abril modificado Mexico
  # vv es matriz de var y cov de est de max vero
  # es análoga a sigma^2 * (X^T X)^{-1} del caso lineal
  # los parámetros están en escalas muy disímiles, se "escala
  # y desescala" a la hora de invertr
  vv     <- sigma^2 * (uu %*% solve(uu%*%aux%*%uu) %*% uu)
  Dlam   <- cbind(1/(hab*alfa),-exp(-kapa*x0),beta*x0*exp(-kapa*x0)) * predic
  aux    <- qnorm(1-(1-conf)/2)
  # se añade sigma^2 para errores de predicción
  zv     <- aux*sqrt(sigma^2 + diag(Dlam %*% vv %*% t(Dlam)))
  lims   <- (predic + zv)
  limi   <- ifelse( (predic - zv)<0, 0, predic - zv ) 
  preds  <- as.data.frame(cbind(limi,predic,lims))
  return(list(preds=preds,Mvarianza=vv)) }

########################### TERMINA PREDICCIÓN ##################



########################### EMPIEZA GRÁFICAS ##################
# Función para la graficación del ajuste de modelo Gompertz

Gomp3 <- function(h, pred, coe=bb$Coeficientes, dfdato=bb$DatosId, VV,GHH,
                  casos=dat2$casos){
  fechau <- tail(dfdato$tiempo,1)
  casosu <- tail(dfdato$casosa,1)
  casosf <- format( casosu, big.mark = ",")
  Nmax   <- (dfdato$CasosTotales)[1]
  Nmaxf  <- format(Nmax,big.mark = ",")
  percen <- round(100*casosu/Nmax)
  casap  <- (dfdato$CasosAcumPico)[1]
  casapf <- format( casap, big.mark = ",")
  cas1S  <- (dfdato$CasosAcum1SyPico)[1]
  cas1Sf <- format( cas1S, big.mark = ",")
  caspm  <- (dfdato$CasosPorMillon)[1]
  caspmf <- format( caspm, big.mark = ",")
  hab    <- (dfdato$hab)[1]
  Nh     <- dim(pred)[1]
  tiempo <- seq.Date(from=dfdato$tiempo[1],to=dfdato$tiempo[1]+(Nh-1),by=1)
  pred$tiempo <- tiempo
  pred$tisD   <- 1:Nh
  paisn       <- (dfdato$Unidad)[1]
  labtiempo   <- function(x){
    as.Date( tiempo[round(x)+1], format="%y.%m.%d" ) }
  pred$limi <- ifelse(pred$tisD < Nh-h,pred$predic,pred$limi)
  pred$lims <- ifelse(pred$tisD < Nh-h,pred$predic,pred$lims)
  # graf1: Datos acumulados con curva estimada y predicciones con int conf
  graf1  <- ggplot( ) +
    xlab("Días") + ylab("casos confirmados acum.") +
    geom_line(data=pred, aes(x=tisD,y=limi), colour="red", linetype=2, size=.4) +
    geom_line(data=pred, aes(x=tisD,y=lims), colour="red", linetype=2, size=.4) +
    geom_line(data=pred, aes(x=tisD,y=predic), colour="red", size=1.05) +
    geom_step(data=dfdato, mapping=aes(x=tt,y=casosa), linetype=1,size=.8,
              direction="hv") +
    ggtitle(paste(paisn,":",fechau,"Casos = ",casosf,"  ( Total Pred.",Nmaxf,")")) +
    scale_x_continuous( labels = labtiempo ) 
  # graf2: Gráficas de derivada de Gompertz
  fgom  <- function(tt,alfa,beta,kapa){
    return(alfa*exp(-beta*exp(-kapa*tt))) }
  dfgom <- function(tt,alfa,beta,kapa){
    return(beta*kapa*exp(-kapa*tt)*fgom(tt,alfa,beta,kapa)) }
  ttt   <- seq(0,GHH,length.out = 301)
  alfa  <- coe[1]
  beta  <- coe[2]
  kapa  <- coe[3]
  tmax  <- log(beta)/kapa
  fmax  <- dfdato$tiempo[1] + round(tmax)
  cota  <- 0.95
  tdura <- -log(-log(cota)/beta)/kapa
  fdura <- tiempo[1] + round(tdura)
  yyf   <- hab*dfgom(tdura,alfa,beta,kapa)
  if(tdura > GHH){ tdura <- GHH-10 }
  datgraf <- data.frame(tt=ttt,ff=hab*dfgom(ttt,alfa,beta,kapa))
  df      <- data.frame(x1 = tmax, y1 = 0, x2 = tmax, 
                        y2 = hab*dfgom(tmax,alfa,beta,kapa),
                        x1f = tdura, y1f = 0, x2f = tdura,
                        y2f = yyf )
  df2     <- data.frame(casos=casos,xx=1:length(casos))
  graf2 <- ggplot( data=datgraf, aes(x=ttt,y=ff) ) + 
    xlab(paste(paisn," (Num. de días desde caso 1)")) + 
    ylab("Num. casos por día") +
    geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data = df) +
    geom_segment(aes(x=x1f, y=y1f, xend=x2f, yend=y2f), data = df) +
    geom_line( data=datgraf, aes(x=tt,y=ff), size=1, colour="blue" ) +
    geom_text(x=tmax,y=0,label=fmax) +
    geom_text(x=tdura,y=0,label=fdura) +
    ggtitle(paste("Acum Pico",casapf,
                  ". Acum 1 Sem y Pico", cas1Sf, ". Casos por Mill",caspmf)) 
  # graf3: Gráficas de derivada de Gompertz con simulaciones de la curva
  # usando la normalidad asintótica de los max vero, para dimensionar
  # el nivel de incertidumbre
  x    <- seq(0,GHH,length.out = 301)
  dddf <- data.frame(x=x,
                     y=hab*coe[1]*exp(-coe[2]*exp(-coe[3]*x))*coe[2]*coe[3]*exp(-coe[3]*x))
  aux  <- diag(c(1/hab,1,1))
  Vcoe <- aux %*% VV %*% aux
  N    <- 100
  tet  <- mvrnorm(n = N, mu=coe, Sigma=Vcoe)
  alf  <- ifelse(tet[,1]<0,0,tet[,1])
  bet  <- ifelse(tet[,2]<0,0,tet[,2])
  kap  <- ifelse(tet[,3]<0,0,tet[,3])
  gompertz_curves <- tibble(Group = 1:N, alf, bet, kap)
  plot_data <- pmap_df(gompertz_curves,
                       function(Group, alf,bet,kap) {
                         tibble(Group = Group,
                                x = seq(0,GHH,length.out = 301),
                                y = hab*alf*exp(-bet*exp(-kap*x))*bet*kap*exp(-kap*x) ) })
  graf3 <- ggplot(data = plot_data) +
    xlab(paste(paisn," (Num. de días desde caso 1)")) + 
    ylab("Num. casos por día") +
    geom_line(aes(group = Group, x = x, y = y), colour="red",alpha = .05) +
    geom_line(data=dddf, aes(x=x, y=y), colour="blue",size=1.2) +
    geom_text(x=tmax,y=0,label=fmax) +
    geom_text(x=tdura,y=0,label=fdura) +
    geom_linerange(mapping = aes(x=xx,y=casos),ymin=0,ymax=casos,data=df2) +
    geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2), colour="red", size=1.1,data = df) +
    geom_segment(aes(x=x1f,y=y1f,xend=x2f,yend=y2f),colour="red",size=1.1,data = df) +
    ggtitle(paste("Acum Pico",casapf,
                  ". Acum 1 Sem y Pico", cas1Sf, ". Casos por Mill",caspmf)) 
  return(list(graf1=graf1,graf2=graf2,graf3=graf3)) }

########################### TERMINA GRÁFICAS ##################





########################### EMPIEZA MEXICO ##################
# Lectura de incidencia acumulada de casos
ff   <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master"
ff   <- paste(ff,"/csse_covid_19_data/csse_covid_19_time_series/",sep="")
ff   <- paste(ff,"time_series_covid19_confirmed_global.csv",sep="")
dat  <- read.csv2(ff,header=TRUE,sep=",")
dim(dat) # 266 x 159
nc   <- dim(dat)[2] 

# Países en la base de datos
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
pob2018 <- 128932753               # Población actual de México
unit    <- pais[id]
zz      <- datu[id,]
zzc     <- diff(zz)
nonn    <- min( which(zzc > 0) )
tiempo  <- dd[-(1:nonn)]
casos   <- zzc[-(1:(nonn-1))]    # no se usa esta variable
casosa  <- datu[id,-(1:nonn)]
dat2    <- data.frame(Unidad=unit, dateRep=tiempo, casos=casos, 
                        casosa=casosa, pob=pob2018)
bb      <- Gomp1(dat2)
h       <- 3
aux     <- Gomp2(h, coe=bb$Coeficientes, sigma=bb$sigma, dfdato=bb$DatosId)
pred    <- aux$preds
VVar    <- aux$Mvarianza
GHH     <- 240
ggra    <- Gomp3(h, pred, coe=bb$Coeficientes, dfdato=bb$DatosId, VV=VVar, GHH=GHH)
gg      <- ggarrange(ggra$graf1,ggra$graf3, ncol = 1, nrow = 2)

plot(gg)





































########################### EMPIEZA MEXICO VERSION 2 ##################
# Lectura de incidencia acumulada de casos
ff   <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master"
ff   <- paste(ff,"/csse_covid_19_data/csse_covid_19_time_series/",sep="")
ff   <- paste(ff,"time_series_covid19_confirmed_global.csv",sep="")
dat  <- read.csv2(ff,header=TRUE,sep=",")
dim(dat) # 266 x 151
nc   <- dim(dat)[2] 

# Países en la base de datos
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

id       <- which(pais=="Mexico")
pob2018  <- 128932753               # Población actual de México
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
for(MM in 61:N){
  i       <- i+1
  tiempo  <- tiempoA[1:MM]
  casos   <- casosA[1:MM]     
  casosa  <- casosaA[1:MM] 
  dat2    <- data.frame(Unidad=unit, dateRep=tiempo, casos=casos, 
                      casosa=casosa, pob=pob2018)
  bb      <- Gomp1(dat2)   # Gomp1 es la misma
  h       <- 300 - i
  aux     <- Gomp2(h, coe=bb$Coeficientes, sigma=bb$sigma, dfdato=bb$DatosId)
  pred    <- aux$preds
  VVar    <- aux$Mvarianza
  GHH     <- 360
  ggra    <- Gomp3X(h, pred, coe=bb$Coeficientes, dfdato=bb$DatosId, VV=VVar, GHH=GHH)
  grafejesfijos[[i]] <- ggarrange(ggra$graf1,ggra$graf3, ncol = 1, nrow = 2)
}

pdf("G_MXEjesFijos.pdf")
for(i in 1:(N-60)){ print(grafejesfijos[[i]]) }
dev.off()




########################### EMPIEZA GRÁFICAS VERSION 2 ##################
# Función para la graficación del ajuste de modelo Gompertz

Gomp3X <- function(h, pred, coe=bb$Coeficientes, dfdato=bb$DatosId, VV,GHH,
                  casos=dat2$casos){
  fechau <- tail(dfdato$tiempo,1)
  casosu <- tail(dfdato$casosa,1)
  casosf <- format( casosu, big.mark = ",")
  Nmax   <- (dfdato$CasosTotales)[1]
  Nmaxf  <- format(Nmax,big.mark = ",")
  percen <- round(100*casosu/Nmax)
  casap  <- (dfdato$CasosAcumPico)[1]
  casapf <- format( casap, big.mark = ",")
  cas1S  <- (dfdato$CasosAcum1SyPico)[1]
  cas1Sf <- format( cas1S, big.mark = ",")
  caspm  <- (dfdato$CasosPorMillon)[1]
  caspmf <- format( caspm, big.mark = ",")
  hab    <- (dfdato$hab)[1]
  Nh     <- dim(pred)[1]
  tiempo <- seq.Date(from=dfdato$tiempo[1],to=dfdato$tiempo[1]+(Nh-1),by=1)
  pred$tiempo <- tiempo
  pred$tisD   <- 1:Nh
  paisn       <- (dfdato$Unidad)[1]
  labtiempo   <- function(x){
    as.Date( tiempo[round(x)+1], format="%y.%m.%d" ) }
  pred$limi <- ifelse(pred$tisD < Nh-h,pred$predic,pred$limi)
  pred$lims <- ifelse(pred$tisD < Nh-h,pred$predic,pred$lims)
  # graf1: Datos acumulados con curva estimada y predicciones con int conf
  graf1  <- ggplot( ) +
    xlab("Días") + ylab("casos confirmados acum.") +
    ylim(0,800000) +
    geom_line(data=pred, aes(x=tisD,y=limi), colour="red", linetype=2, size=.4) +
    geom_line(data=pred, aes(x=tisD,y=lims), colour="red", linetype=2, size=.4) +
    geom_line(data=pred, aes(x=tisD,y=predic), colour="red", size=1.05) +
    geom_step(data=dfdato, mapping=aes(x=tt,y=casosa), linetype=1,size=1.1,
              direction="hv") +
    ggtitle(paste(paisn,":",fechau,"Casos = ",casosf,"  ( Total Pred.",Nmaxf,")")) +
    scale_x_continuous( labels = labtiempo ) 
  # graf2: Gráficas de derivada de Gompertz
  fgom  <- function(tt,alfa,beta,kapa){
    return(alfa*exp(-beta*exp(-kapa*tt))) }
  dfgom <- function(tt,alfa,beta,kapa){
    return(beta*kapa*exp(-kapa*tt)*fgom(tt,alfa,beta,kapa)) }
  ttt   <- seq(0,GHH,length.out = 301)
  alfa  <- coe[1]
  beta  <- coe[2]
  kapa  <- coe[3]
  tmax  <- log(beta)/kapa
  fmax  <- dfdato$tiempo[1] + round(tmax)
  cota  <- 0.95
  tdura <- -log(-log(cota)/beta)/kapa
  fdura <- tiempo[1] + round(tdura)
  yyf   <- hab*dfgom(tdura,alfa,beta,kapa)
  if(tdura > GHH){ tdura <- GHH-10 }
  datgraf <- data.frame(tt=ttt,ff=hab*dfgom(ttt,alfa,beta,kapa))
  df      <- data.frame(x1 = tmax, y1 = 0, x2 = tmax, 
                        y2 = hab*dfgom(tmax,alfa,beta,kapa),
                        x1f = tdura, y1f = 0, x2f = tdura,
                        y2f = yyf )
  df2     <- data.frame(casos=casos,xx=1:length(casos))
  graf2 <- ggplot( data=datgraf, aes(x=ttt,y=ff) ) + 
    xlab(paste(paisn," (Num. de días desde caso 1)")) + 
    ylab("Num. casos por día") +
    geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), data = df) +
    geom_segment(aes(x=x1f, y=y1f, xend=x2f, yend=y2f), data = df) +
    geom_line( data=datgraf, aes(x=tt,y=ff), size=1, colour="blue" ) +
    geom_text(x=tmax,y=0,label=fmax) +
    geom_text(x=tdura,y=0,label=fdura) +
    ggtitle(paste("Acum Pico",casapf,
                  ". Acum 1 Sem y Pico", cas1Sf, ". Casos por Mill",caspmf)) 
  # graf3: Gráficas de derivada de Gompertz con simulaciones de la curva
  # usando la normalidad asintótica de los max vero, para dimensionar
  # el nivel de incertidumbre
  x    <- seq(0,GHH,length.out = 301)
  dddf <- data.frame(x=x,
                     y=hab*coe[1]*exp(-coe[2]*exp(-coe[3]*x))*coe[2]*coe[3]*exp(-coe[3]*x))
  aux  <- diag(c(1/hab,1,1))
  Vcoe <- aux %*% VV %*% aux
  N    <- 100
  tet  <- mvrnorm(n = N, mu=coe, Sigma=Vcoe)
  alf  <- ifelse(tet[,1]<0,0,tet[,1])
  bet  <- ifelse(tet[,2]<0,0,tet[,2])
  kap  <- ifelse(tet[,3]<0,0,tet[,3])
  gompertz_curves <- tibble(Group = 1:N, alf, bet, kap)
  plot_data <- pmap_df(gompertz_curves,
                       function(Group, alf,bet,kap) {
                         tibble(Group = Group,
                                x = seq(0,GHH,length.out = 301),
                                y = hab*alf*exp(-bet*exp(-kap*x))*bet*kap*exp(-kap*x) ) })
  graf3 <- ggplot(data = plot_data) +
    xlab(paste(paisn," (Num. de días desde caso 1)")) + 
    ylab("Num. casos por día") + ylim(0,5500) +
    geom_linerange(mapping = aes(x=xx,y=casos),ymin=0,ymax=casos,data=df2) +
    geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2),size=1.2,colour="red",data = df) +
    geom_segment(aes(x=x1f,y=y1f,xend=x2f,yend=y2f),size=1.2,colour="red",data = df) +
    geom_line(aes(group = Group, x = x, y = y), colour="red",alpha = .05) +
    geom_line(data=dddf, aes(x=x, y=y), colour="blue",size=1.2) +
    geom_text(x=tmax,y=0,label=fmax) +
    geom_text(x=tdura,y=0,label=fdura) +
    geom_linerange(mapping = aes(x=xx,y=casos),ymin=0,ymax=casos,data=df2) +
    ggtitle(paste("Acum Pico",casapf,
                  ". Acum 1 Sem y Pico", cas1Sf, ". Casos por Mill",caspmf)) 
  return(list(graf1=graf1,graf2=graf2,graf3=graf3)) }

########################### TERMINA GRÁFICAS VERSION 2 ##################

