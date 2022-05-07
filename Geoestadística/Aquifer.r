# Aquifer

###############################################################################

rm(list = ls())

################################################################################
################################ paquetes ######################################
################################################################################

# install.packages("stars")

################################################################################
############################### librerias ######################################
################################################################################

library(geoR)
library(fields)
library(akima)
library(dplyr)


################################################################################
#################################### path ######################################
################################################################################

getCurrentFileLocation <-  function()
{
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}


path = getCurrentFileLocation()
setwd(path)
getwd()

################################################################################
############################ Lectura de datos ##################################
################################################################################

aquifer <- read.table("aquifer.txt", head = TRUE, dec = ",")
head(aquifer)
summary(aquifer)

names(aquifer)[1] <- "Este"
names(aquifer)[2] <- "Norte"
names(aquifer)[3] <- "Profundidad"


################################################################################
############### Convertir aquifer a un objeto geodata (geoR obj) ###############
################################################################################

aquiferg <- as.geodata(aquifer)
summary(aquiferg)

#-----------------------------------------------------------------------------#
# Gráfico de objeto geodata
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Gráfico del objeto geodata

windows();
plot(aquiferg, qt.col = c("purple",
                         "pink",
                         "green",
                         "yellow"))

#-----------------------------------------------------------------------------#
# Gráfico con el parametro 3d

windows();
plot(aquiferg, scatter3d = T)

#-----------------------------------------------------------------------------#
# Gráfico removiendo la tendencia (trend )

windows();
plot(aquiferg, trend = "1st")


#-----------------------------------------------------------------------------#
# Gráficos descriptivos interpolación
#-----------------------------------------------------------------------------#

windows();
par(mfrow = c(2, 2),
    mar = c(3, 3, 1, 1),
    mgp = c(2, 1, 0))
# Esta función agrupa los siguientes gráficos en
# una matrix 2x2

grillas <- interp(aquifer$Este,
                  aquifer$Norte,
                  aquifer$Profundidad)

persp(grillas$x,
      grillas$y,
      grillas$z,
      xlab = "Este",
      ylab = "Norte",
      zlab = "Nivel freatico",
      phi = 30,
      theta = 20,
      col = "lightblue",
      expand = .5,
      ticktype = "detailed")

drape.plot(grillas$x,
           grillas$y,
           grillas$z,
           xlab = "Este",
           ylab = "Norte",
           zlab = "z",
           theta = 45,
           col = topo.colors(64),
           expand = .5,
           ticktype = "detailed")


drape.plot(grillas$x,
           grillas$y,
           grillas$z,
           xlab = "Este",
           ylab = "Norte",
           zlab = "z",
           theta = -10,
           col = topo.colors(64),
           expand = .5,
           ticktype = "detailed")


drape.plot(grillas$x,
           grillas$y,
           grillas$z,
           xlab = "Este",
           ylab = "Norte",
           zlab = "z",
           theta = 60,
           col = topo.colors(64),
           expand = .5,
           ticktype = "detailed")


#-----------------------------------------------------------------------------#
# Gráficos de contorno
#-----------------------------------------------------------------------------#

windows();
par(mfrow = c(2, 1),
    mar = c(1,1,1,1))

contour(grillas, nlevels = 10, main = "Contorno")
image(grillas$z, main =  "Grilla")

filled.contour(grillas, levels = seq(1000,
                                     5000,
                                     len = 10),
               col = heat.colors(10),
                main = "grilla niveles")

#-----------------------------------------------------------------------------#
# Funciones y gráficas a partir de la función outer
#-----------------------------------------------------------------------------#

h <- seq(0, 1, len = 50)
u <- seq(0, 1, len = 50)

ejemplo1CH  <- function(h, u, sigma, a, b, c, d, delta) {
    (sigma^2/((a^2*u^2+c)^(d/2)))*exp(-(b^2*h^2)/(a^2*u^2+c))*exp(-delta*u^2)
    }
h <- seq(0, 1, len = 20)
u <- seq(1, 10, len = 20)
f <- outer(h, u, ejemplo1CH, sigma=3, a=1, b=3, c=1, d=2, delta=0)

windows();
par(mfrow = c(2, 2),
    mar = c(3, 3, 1, 1),
    mgp = c(2, 1, 0))

drape.plot(h,
           u,
           f,
           main = "Cressie-Huang; 1 (25,1,0.6)",
           xlab = "h",
           ylab = "u",
           zlab = "Covarianza",
           ltheta = 75,
           col = terrain.colors(64))

drape.plot(h,
           u,
           f,
           main = "Cressie-Huang; 1 (25,1,0.6)",
           xlab = "h",
           ylab = "u",
           zlab = "Covarianza",
           theta = -150,
           col = terrain.colors(64))
persp(h,
      u,
      f,
      main = "Cressie-Huang; 1 (25,1,0.6)",
      xlab = "h",
      ylab = "u",
      zlab = "Covarianza",
      ltheta = 75)

contour(h,
        u,
        f,
        col = topo.colors(10),
        xlim = c(0,0.6))


#-----------------------------------------------------------------------------#
# Modelando la media con regresión polinomial
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Primer modelo

reg1 <- lm(Profundidad ~ Este + Norte, data = aquifer)
residuales1  <-  residuals(reg1)
summary(reg1)
anova(reg1)


#-----------------------------------------------------------------------------#
# Segundo modelo

reg2 <- lm(Profundidad ~ Este + Norte +
           I(Este^2) + I(Norte^2) +
           I(Este * Norte),
           data = aquifer)
residuales2  <-  residuals(reg2)
summary(reg2)
anova(reg2)


#-----------------------------------------------------------------------------#
# Tercer modelo

reg3 <- lm(Profundidad ~ Este * Norte,
           data = aquifer)
residuales3  <-  residuals(reg3)
summary(reg3)
anova(reg3)


#-----------------------------------------------------------------------------#
# Estimación del semivariograma empírico
#-----------------------------------------------------------------------------#

vari2 <- variog(aquiferg, trend = "1st")

vari2Cloud <- variog(aquiferg, op = "cloud", trend = "1st")

vari2BinCloud <- variog(aquiferg,
                       max.dist = 200,
                       op = "cloud",
                       bin.cloud = TRUE)

vari2Sm <- variog(aquiferg,
                  trend = "1st",
                  op = "sm",
                  band=11)

windows();
par(mfrow = c(2, 2), mar = c(3, 3, 1, 1), mgp = c(2, 1, 0))
     plot(vari2, main = "binned variogram")
     plot(vari2Cloud, main = "variogram cloud")
     plot(vari2BinCloud,main = "clouds for binned variogram")
     plot(vari2Sm, main = "smoothed variogram")

#-----------------------------------------------------------------------------#
# Explorando estimación clásica, removiendo tendencia
#-----------------------------------------------------------------------------#

vari1 <- variog(aquiferg)
vari2 <- variog(aquiferg, trend = "1st")
vari3 <- variog(aquiferg, trend = "2nd")

#-----------------------------------------------------------------------------#
# sin remover tendencia 

windows();
plot(vari1, main = "Sin remover tendencia")

#-----------------------------------------------------------------------------#
# removiendo la primera tendencia

windows();
plot(vari2, main  = "Trend 1 ")

#-----------------------------------------------------------------------------#
# removiendo la segunda tendencia

windows();
plot(vari3, main  = "Trend 2 ")

#-----------------------------------------------------------------------------#
# Explorando estimación resistente a datos atípicos y removiendo tendencia
#-----------------------------------------------------------------------------#

vari1 <- variog(aquiferg, estimator.type = "modulus")
vari2 <- variog(aquiferg, trend = "1st", estimator.type = "modulus")
vari3 <- variog(aquiferg, trend = "2nd", estimator.type = "modulus")

#-----------------------------------------------------------------------------#
# sin remover tendencia 

windows();
plot(vari1, main = "Sin remover tendencia")

#-----------------------------------------------------------------------------#
# removiendo la primera tendencia

windows();
plot(vari2, main  = "Trend 1 ")

#-----------------------------------------------------------------------------#
# removiendo la segunda tendencia

windows();
plot(vari3, main  = "Trend 2 ")


#-----------------------------------------------------------------------------#
# Explorando anisotropía
#-----------------------------------------------------------------------------#

vari_0 <- variog(aquiferg,
                 trend = "1st",
                 max.dist = 200,
                 dir = 0)

vari_45 <- variog(aquiferg,
                  trend = "1st",
                  max.dist = 200,
                  dir = pi / 4)
vari_90 <- variog(aquiferg,
                  trend = "1st",
                  max.dist = 200,
                  dir = pi / 2)
vari_135 <- variog(aquiferg,
                   trend = "1st",
                   max.dist = 200,
                   dir = 3 * pi / 4)


windows();
par(mfrow = c(2, 2),
    mar = c(3, 3, 1, 1),
    mgp = c(2, 1, 0))

plot(vari_0, main = "vari 0")
plot(vari_45, main = "vari 45")
plot(vari_90, main = "vari 90")
plot(vari_135, main = "vari 195")


#-----------------------------------------------------------------------------#
# Estimación teórica del semivariograma.
#-----------------------------------------------------------------------------#

var1 <- variog(aquiferg,trend="1st",max.dist=200)


#ini1 <- eyefit(var1)
#cov.model  sigmasq phi   tausq kappa kappa2   practicalRange
#1      wave 30805.52  13 8984.94  <NA>   <NA> 38.8889336320589
ini1 <- c(30805.52, 13)
fitvar1 <- variofit(var1,
                    cov.model = "wave",
                    ini1,
                    fix.nugget = TRUE,
                    nugget = 8984.94,
                    wei = "equal")

fitvar2 <- variofit(var1,
                    cov.model = "wave",
                    ini1,
                    fix.nugget = TRUE,
                    nugget = 8984.94,
                    wei = "npairs")

fitvar3 <- variofit(var1,
                    ini1,
                    fix.nugget = TRUE,
                    nugget = 8984.94,
                    wei = "cressie")


fitvar4 <- likfit(aquiferg,
                  coords = aquiferg$coords,
                  data = aquiferg$data,
                  trend = "1st",
                  ini.cov.pars = ini1,
                  fix.nugget = T,
                  nugget = 8984.94,
                  cov.model = "wave",
                  lik.method = "ML")

fitvar5 <- likfit(aquiferg,
                  coords = aquiferg$coords,
                  data = aquiferg$data,
                  trend = "1st",
                  ini.cov.pars = ini1,
                  fix.nugget = T,
                  nugget = 8984.94,
                  cov.model = "wave",
                  lik.method = "REML")

windows();
plot(var1,
     xlab = "h",
     ylab = "semivarianza",
     cex.lab = 1.3,
     cex.axis = 1.2,
     main = "Estimación teórica del modelo de semivariograma",
     col.main = 4, cex.main =1.3)
lines(fitvar1, col = 1)
lines(fitvar2, col = 2)
lines(fitvar3, col = 3)
lines(fitvar4, col = 4)
lines(fitvar5, col = 5)
legend(130, 18000,
       c("MCO", "MCPnpairs", "MCPcressie", "ML", "REML"),
       lwd = 2,
       lty = 2:7,
       col = 2:7,
       box.col = 9,
       text.col = 2:7)


################################################################################
#################################### Kriging ###################################
################################################################################

library(sf)
library(gstat)
library(stars)

head(aquifer)

#-----------------------------------------------------------------------------#
# empleando KO (kriging ordinario)

# aquifer$Profundidad <- aquifer$Profundidad/100 # en cientos de pies...
aquifer_sf <- st_as_sf(aquifer, coords = c("Este", "Norte"), remove = FALSE, agr = "constant")

vario <- variogram(Profundidad ~ Este + Norte, aquifer_sf, cutoff = 200)
fit <- fit.variogram(vario, vgm(model = "Sph", nugget = NA), fit.method = 2)

windows();
plot(vario)


#-----------------------------------------------------------------------------#
# Grilla

buffer <- aquifer_sf %>% st_geometry() %>% st_buffer(40)

grid <- buffer %>%  st_as_stars(nx = 50, ny = 50)


#-----------------------------------------------------------------------------#
# Como suponemos un modelo (no constante) para la tendencia, 
# es necesario añadir los valores de las variables explicativas
# a la rejilla de predicción:

coord <- st_coordinates(grid)
grid$Este <- coord$x
grid$Norte <- coord$y


#-----------------------------------------------------------------------------#
# recortamos la rejeilla

grid <- grid %>% st_crop(buffer)


#-----------------------------------------------------------------------------#
# predicciones mediante kriging universal

pred <- krige(Profundidad ~ Este + Norte, aquifer_sf, model = fit,
              newdata = grid)

#-----------------------------------------------------------------------------#
# coordenadas del objeto  

summary(st_coordinates(grid))

summary(st_coordinates(pred))

grid$var1.pred <- pred$var1.pred
grid$var1.var <- pred$var1.var


#-----------------------------------------------------------------------------#
# grafico de las predicciones y las varianzas kriging 

windows();
plot(grid["var1.pred"], breaks = "equal", col = sf.colors(64), key.pos = 4,
     main = "Predicciones kriging")

windows();
plot(grid["var1.var"], breaks = "equal", col = sf.colors(64), key.pos = 4,
     main = "Varianzas kriging")


#-----------------------------------------------------------------------------#
# con el paquete ggplot2

library(ggplot2)
library(gridExtra)


p1 <- ggplot() + geom_stars(data = grid, aes(fill = var1.pred, x = x, y = y)) +
    scale_fill_viridis_c() + geom_sf(data = aquifer_sf) 
    # + coord_sf(lims_method = "geometry_bbox")


p2 <- ggplot() + geom_stars(data = grid, aes(fill = var1.var, x = x, y = y)) +
    scale_fill_viridis_c() + geom_sf(data = aquifer_sf) 
    #+ coord_sf(lims_method = "geometry_bbox")

windows();
grid.arrange(p1, p2, ncol = 2)


################################################################################
########################### validacion cruzada #################################
################################################################################

system.time(cv <- krige.cv(formula = Profundidad ~ Este + Norte, locations = aquifer_sf,
                           model = fit))

str(cv)

#-----------------------------------------------------------------------------#
#  partición en grupos es aleatoria
# Un parámetro opcional que establece el número de pliegues para la validación cruzada.

set.seed(1)
system.time(cv <- krige.cv(formula = Profundidad ~ Este + Norte, locations = aquifer_sf,
                           model = fit, nfold = 10))

#-----------------------------------------------------------------------------#
# estadísticos

summary_cv <- function(cv.data, na.rm = FALSE,
                       tol = sqrt(.Machine$double.eps)) {
  err <- cv.data$residual      # Errores
  obs <- cv.data$observed
  z <- cv.data$zscore
  w <- 1/pmax(cv.data$var1.var, tol) # Ponderación según varianza kriging
  if(na.rm) {
    is.a <- !is.na(err)
    err <- err[is.a]
    obs <- obs[is.a]
    z <- z[is.a]
    w <- w[is.a]
  }
  perr <- 100*err/pmax(obs, tol)  # Errores porcentuales
  return(c(
    # Medidas de error tradicionales
    me = mean(err),           # Error medio
    rmse = sqrt(mean(err^2)), # Raíz del error cuadrático medio
    mae = mean(abs(err)),     # Error absoluto medio
    mpe = mean(perr),         # Error porcentual medio
    mape = mean(abs(perr)),   # Error porcentual absoluto medio
    r.squared = 1 - sum(err^2)/sum((obs - mean(obs))^2), # Pseudo R-cuadrado
    # Medidas de error que tienen en cuenta la varianza kriging
    dme = mean(z),            # Error estandarizado medio
    dmse = sqrt(mean(z^2)),    # Error cuadrático medio adimensional
    rwmse = sqrt(weighted.mean(err^2, w)) # Raíz del ECM ponderado
  ))
}

summary_cv(cv)
