##################################################################
##                                                              ##
## Carlos Serrano Fajardo                                       ##
##                                                              ##
## Grado en Bioquímica                                          ##
## Dpto. Nutrición y Bromatología, Toxicología y Medicina Legal ##
## Grupo Color y Calidad de Alimentos                           ##
##                                                              ##
## ANÁLISIS TEXTUROMÉTRICO Y DE IMAGEN DE SEMILLAS DE UVA       ##
## TRABAJO DE FIN DE GRADO - SCRIPT III                         ##
##                                                              ##
## Esta obra está sujeta a la licencia                          ##
## Reconocimiento-NoComercial-CompartirIgual                    ##
## 4.0 Internacional de Creative Commons. Para                  ##
## ver una copia de esta licencia, visite                       ##
## http://creativecommons.org/licenses/by-nc-sa/4.0/.           ##
##                                                              ##
##################################################################

## En este archivo se realiza el estudio de la variación de cada una de las propiedades
## texturométricas, de imagen o hiperespectrales de las muestras a lo largo de un período
## de maduración en planta. También se realiza un estudio similar considerando semillas que
## han sufrido un proceso de desecación en uvas colocadas en paseras. Posteriormente, se
## estudia cómo las diferentes variables pueden estar correlacionadas entre sí, tanto de forma
## global como categorizando las muestras por condiciones.

##########################
## PAQUETES Y FUNCIONES ##
##########################

library(reshape)
library(scales)
library(plyr)
library(gdata)
library(Hmisc)
library(corrplot)

setwd("~/Académico/Universidad/04_TFG/03_DATOS/output1/")

####################################
####################################
####                            ####
#### MADURACIÓN TEXTURA-TIEMPO  ####
####                            ####
####################################
####################################

## PEDRO XIMÉNEZ; 2 PASERAS CONJUNTAS
## Para esto sólo seleccionamos las muestras LP01X - LP17X. Colocamos las muestras a analizar y
## paralelamente el tiempo que ha transcurrido en su toma. Se le asigna día 0 a la primera muestra.

texture_data <- c("output_LP01X.txt", "output_LP03X.txt", "output_LP05X.txt", "output_LP06X.txt",
                  "output_LP07X.txt", "output_LP10X.txt", "output_LP11X.txt", "output_LP17X.txt",
                  "output_LP12X.txt", "output_LP14X.txt", "output_LP15X.txt")
time_data <- c(0, 1, 3, 4, 6, 15, 35, 36, 37, 41, 45)

texture_data_list <- list()
texture_means_list <- list()
texture_sds_list <- list()
for(i in 1:length(texture_data))
{
   texture_data_list[[i]] <- read.csv(file=texture_data[i], sep="", header=TRUE)
   texture_means_list[[i]] <- colMeans(texture_data_list[[i]])
   texture_sds_list[[i]] <- sapply(texture_data_list[[i]], sd)
}
names(texture_means_list) <- texture_data
names(texture_sds_list) <- texture_data

## Para parámetros de los que sí se tienen las desviaciones estándar.
variation_in_time <- function(texture_means_list, texture_sds_list, time_data, parameter_name, main=main, ylab=ylab)
{
   means <- numeric()
   sds <- numeric()
   for(i in 1:length(texture_means_list))
   {
      means[[i]] <- texture_means_list[[i]][[parameter_name]]
      sds[[i]] <- texture_sds_list[[i]][[parameter_name]]
   }
   data <- data.frame(time_data, means, sds)
   with (data = data, expr = errbar(time_data, means,
                                    means+sds, means-sds, col="blue", cex=1.2, pch=16, cap=0.015,
                                    add=F, ylab=ylab, xlab="Tiempo (días)"))
   title(main=main)
   
   #abline(lm(means ~ time_data), col="red")
   #text(text_x, text_y, paste("R-squared:", round(summary(lm(means ~ time_data))[["r.squared"]], 4)))
   return(paste("R-squared:", round(summary(lm(means ~ time_data))[["r.squared"]], 4)))
}

## Para parámetros de los que no se tienen las desviaciones estándar.
variation_in_time2 <- function(texture_means_list, time_data, parameter_name, main=main, ylab=ylab,text_x, text_y)
{
   means <- numeric()
   for(i in 1:length(texture_means_list))
   {
      means[[i]] <- texture_means_list[[i]][[parameter_name]]
   }
   
   plot(time_data, means, col="blue", cex=1.2, pch=16,
        ylab=ylab, xlab="Tiempo (días)", main=main)
   return(paste("R-squared:", round(summary(lm(means ~ time_data))[["r.squared"]], 4)))
}

variation_in_time(texture_means_list, texture_sds_list, time_data,
                  "break_force", main="Fuerza de rotura", ylab="Fuerza (N)")

variation_in_time(texture_means_list, texture_sds_list, time_data,
                  "break_step_force", main="Fuerza de caída", ylab="Fuerza (N)")

variation_in_time(texture_means_list, texture_sds_list, time_data,
                  "break_deformation", main="Deformación de rotura", ylab="Deformación (%)")

variation_in_time(texture_means_list, texture_sds_list, time_data,
                  "break_energy", main="Energía de rotura", ylab="Energía (mJ)")

variation_in_time(texture_means_list, texture_sds_list, time_data,
                  "slope", main="Pendiente", ylab="Pendiente (N/mm)")

variation_in_time(texture_means_list, texture_sds_list, time_data,
                  "total_energy", main="Energía de compresión 50%", ylab="Energía (mJ)")

variation_in_time(texture_means_list, texture_sds_list, time_data,
                  "seed_width", main="Anchura de semilla", ylab="Anchura (mm)")

######################
## ESTUDIO PASERA 1 ##
## PEDRO XIMÉNEZ    ##
######################

texture_data_pas1 <- c("output_LP01X.txt", "output_LP03X.txt", "output_LP05X.txt", "output_LP06X.txt",
                  "output_LP07X.txt")
time_data_pas1 <- c(0, 1, 3, 4, 6)

texture_data_pas1_list <- list()
texture_means_pas1_list <- list()
texture_sds_pas1_list <- list()
for(i in 1:length(texture_data_pas1))
{
   texture_data_pas1_list[[i]] <- read.csv(file=texture_data_pas1[i], sep="", header=TRUE)
   texture_means_pas1_list[[i]] <- colMeans(texture_data_pas1_list[[i]])
   texture_sds_pas1_list[[i]] <- sapply(texture_data_pas1_list[[i]], sd)
}
names(texture_means_pas1_list) <- texture_data_pas1
names(texture_sds_pas1_list) <- texture_data_pas1

## Sobre esta instrucción, variar el parámetro de entrada para seleccionar otras visualizaciones.
variation_in_time(texture_means_pas1_list, texture_sds_pas1_list, time_data_pas1,
                  "break_force", main="Fuerza de rotura", ylab="Fuerza (N)")

######################
## ESTUDIO PASERA 2 ##
## PEDRO XIMÉNEZ    ##
######################

texture_data_pas2 <- c("output_LP10X.txt", "output_LP11X.txt", "output_LP17X.txt",
                  "output_LP12X.txt", "output_LP14X.txt", "output_LP15X.txt")
time_data_pas2 <- c(15, 35, 36, 37, 41, 45)

texture_data_pas2_list <- list()
texture_means_pas2_list <- list()
texture_sds_pas2_list <- list()
for(i in 1:length(texture_data_pas2))
{
   texture_data_pas2_list[[i]] <- read.csv(file=texture_data_pas2[i], sep="", header=TRUE)
   texture_means_pas2_list[[i]] <- colMeans(texture_data_pas2_list[[i]])
   texture_sds_pas2_list[[i]] <- sapply(texture_data_pas2_list[[i]], sd)
}
names(texture_means_pas2_list) <- texture_data_pas2
names(texture_sds_pas2_list) <- texture_data_pas2

## Sobre esta instrucción, variar el parámetro de entrada para seleccionar otras visualizaciones.
variation_in_time(texture_means_pas2_list, texture_sds_pas2_list, time_data_pas2,
                  "break_force", main="Fuerza de rotura", ylab="Fuerza (N)")

############################################
## ESTUDIO MADURACIÓN EN PLANTA           ##
## PEDRO XIMÉNEZ, 31/08/2015 - 22/09/2015 ##
############################################

texture_data <- c("output_LM03X.txt", "output_LM07X.txt","output_LM09X.txt")
time_data <- c(0, 14, 22)

texture_data_list <- list()
texture_means_list <- list()
texture_sds_list <- list()

for(i in 1:length(texture_data))
{
   texture_data_list[[i]] <- read.csv(file=texture_data[i], sep="", header=TRUE)
   texture_means_list[[i]] <- colMeans(texture_data_list[[i]])
   texture_sds_list[[i]] <- sapply(texture_data_list[[i]], sd)
}
names(texture_means_list) <- texture_data
names(texture_sds_list) <- texture_data

variation_in_time(texture_means_list, texture_sds_list,  time_data,
                  "break_force", main="Fuerza de rotura", ylab="Fuerza (N)")
variation_in_time(texture_means_list, texture_sds_list,  time_data,
                  "break_step_force", main="Fuerza de caída", ylab="Fuerza (N)")
variation_in_time(texture_means_list, texture_sds_list,  time_data,
                  "break_deformation", main="Deformación de rotura", ylab="Deformación (%)")
variation_in_time(texture_means_list, texture_sds_list,  time_data,
                  "break_energy", main="Energía de rotura", ylab="Energía (mJ)")
variation_in_time(texture_means_list, texture_sds_list,  time_data,
                  "slope", main="Pendiente", ylab="Pendiente (N/mm)")
variation_in_time(texture_means_list, texture_sds_list,  time_data,
                  "total_energy", main="Energía de compresión", ylab="Energía (mJ)")
variation_in_time(texture_means_list, texture_sds_list,  time_data,
                  "seed_width", main="Anchura de semilla", ylab="Anchura (mm)")


##############################################
## ESTUDIO SOBREMADURAS                     ##
## MOSCATEL MORISCO; 2/08/2016 - 22/08/2015 ##
##############################################

texture_data <- c("output_LM04M.txt", "output_LM08M.txt","output_LM10M.txt")
time_data <- c(0, 12, 20)

texture_data_list <- list()
texture_means_list <- list()
texture_sds_list <- list()
for(i in 1:length(texture_data))
{
   texture_data_list[[i]] <- read.csv(file=texture_data[i], sep="", header=TRUE)
   texture_means_list[[i]] <- colMeans(texture_data_list[[i]])
   texture_sds_list[[i]] <- sapply(texture_data_list[[i]], sd)
}
names(texture_means_list) <- texture_data
names(texture_sds_list) <- texture_data

variation_in_time2(texture_means_list, time_data,
                   "break_force", main="Fuerza de rotura", ylab="Fuerza (N)")
variation_in_time2(texture_means_list, time_data,
                   "break_step_force", main="Fuerza de caída", ylab="Fuerza (N)")
variation_in_time2(texture_means_list, time_data,
                   "break_deformation", main="Deformación de rotura", ylab="Deformación (%)")
variation_in_time2(texture_means_list, time_data,
                   "break_energy", main="Energía de rotura", ylab="Energía (mJ)")
variation_in_time2(texture_means_list, time_data,
                   "slope", main="Pendiente", ylab="Pendiente (N/mm)")
variation_in_time2(texture_means_list, time_data,
                   "total_energy", main="Energía de compresión 50%", ylab="Energía (mJ)")
variation_in_time2(texture_means_list, time_data,
                   "seed_width", main="Anchura de semilla", ylab="Anchura (mm)")

####################################
####################################
####                            ####
#### MADURACIÓN IMAGEN-TIEMPO   ####
####                            ####
####################################
####################################

## PEDRO XIMÉNEZ; 2 PASERAS CONJUNTAS

maduration_samples <- c("LP01X", "LP03X", "LP05X", "LP06X",
                  "LP07X", "LP10X", "LP11X", "LP17X",
                  "LP12X", "LP14X", "LP15X")

image_data <- read.csv(file="imagen_semillas_puntos.txt", sep="", header=TRUE)
rownames(image_data) <- image_data[,1]

image_selected_data <- image_data[maduration_samples,]
time_data <- c(0, 1, 3, 4, 6, 15, 35, 36, 37, 41, 45)

image_variation_in_time <- function(image_selected_data, time_data, parameter, x_text=35, y_text, main)
{
   plot(time_data, image_selected_data[,parameter], type="p", cex=1.2,
        pch=16, col="blue", xlab="Tiempo (días)", ylab=parameter, main=main)
   #abline(lm(image_selected_data[,parameter] ~ time_data), col="red")
   #text(x_text, y_text, paste("R-squared:", round(summary(lm(image_selected_data[,parameter] ~ time_data))[["r.squared"]], 4)))
   return(paste("R-squared:", round(summary(lm(image_selected_data[,parameter] ~ time_data))[["r.squared"]], 3)))
}


## De entre todos los pa´rametros posibles, descartamos las desviaciónes (que se incluyen
## en el gráfico correspondiente en forma de barras de error) y los parámetros colorimétricos
## en formato L.a.b (lo reemplazamos por el L.h.c).

image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="MajorAxis", x_text=30, y_text=3.0)
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="MinorAxis", x_text=25, y_text=11.3)
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="Area", x_text=25, y_text=0.80)
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="Eccentricity", x_text=25, y_text=3.78)
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="EquivDiameter", x_text=25, y_text=13.4)
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="Perimeter", x_text=25, y_text=1.29)
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="Roundness", x_text=25, y_text=120)
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="L.", x_text=25, y_text=6.1)
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="desv_L.", x_text=25, y_text=11)
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="a.", x_text=25, y_text=2.6)
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="desv_a.", x_text=25, y_text=15.6)
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="b.", x_text=25, y_text=4.6)
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="desv_b.", x_text=25, y_text=19)
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="C.ab", x_text=25, y_text=4.6)
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="desv_C.ab", x_text=25, y_text=56)
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="hab", x_text=30, y_text=8.5)
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="desv_hab", x_text=25, y_text=6.9)

############################################
## ESTUDIO MADURACIÓN EN PLANTA           ##
## PEDRO XIMÉNEZ, 31/08/2015 - 22/09/2015 ##
############################################

maduration_samples <- c("LM03X", "LM07X", "LM09X")

image_data <- read.csv(file="imagen_semillas_puntos.txt", sep="", header=TRUE)

help(read.csv)
rownames(image_data) <- image_data[,1]

image_selected_data <- image_data[maduration_samples,]
time_data <- c(0, 14, 22)

image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="MajorAxis", x_text=30, y_text=3.0, main="Eje mayor")
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="MinorAxis", x_text=25, y_text=11.3, main="Eje menor")
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="Area", x_text=25, y_text=0.80, main="Área")
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="Eccentricity", x_text=25, y_text=3.78, main="Excentricidad")
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="EquivDiameter", x_text=25, y_text=13.4, main="Eje mayor")
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="Perimeter", x_text=25, y_text=1.29, main="Perímetro")
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="Roundness", x_text=25, y_text=120, main="Redondez")
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="L.", x_text=25, y_text=6.1, main="L (claridad)")
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="desv_L.", x_text=25, y_text=11, main="Eje mayor")
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="a.", x_text=25, y_text=2.6)
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="desv_a.", x_text=25, y_text=15.6)
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="b.", x_text=25, y_text=4.6)
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="desv_b.", x_text=25, y_text=19)
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="C.ab", x_text=25, y_text=4.6, main="C (croma)")
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="desv_C.ab", x_text=25, y_text=56, main="Eje mayor")
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="hab", x_text=30, y_text=8.5, main="H (tono)")
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="desv_hab", x_text=25, y_text=6.9, main="Eje mayor")


## Realizamos la representación gráfica con barras de error de los parámetros que lo permiten,
# los parámetros colorimétricos L, C y H.

## L.

means <- c(42.20206, 42.53953, 41.62059)
sds <- c(4.274364, 5.360626, 5.136029)
data <- data.frame(time_data, means, sds)
with (data = data, expr = errbar(time_data, means, means+sds, means-sds, col="blue",
                                 cex=1.2, pch=16, cap=0.015, ylab="L.",
                                 xlab="Tiempo (días)"))
title(main="L (claridad)")

## C.

means <- c(18.99770, 18.75979, 16.94469)
sds <- c(3.243087, 4.079245, 2.796113)
data <- data.frame(time_data, means, sds)
with(data = data, expr = errbar(time_data, means, means+sds, means-sds, col="blue",
                                cex=1.2, pch=16, cap=0.015, ylab="C.ab",
                                xlab="Tiempo (días)"))
title(main="C (croma)")

## H.

means <- c(62.91745, 60.37192, 55.31861)
sds <- c(5.984196, 6.812011, 6.943222)
data <- data.frame(time_data, means, sds)
with (data = data, expr = errbar(time_data, means, means+sds, means-sds, col="blue",
                                 cex=1.2, pch=16, cap=0.015, ylab="h.ab",
                                 xlab="Tiempo (días)"))
title(main="H (tono)")


##############################################
## ESTUDIO SOBREMADURAS                     ##
## MOSCATEL MORISCO; 2/08/2016 - 22/08/2015 ##
##############################################

maduration_samples <- c("LM04M", "LM08M", "LM10M")

image_data <- read.csv(file="imagen_semillas_puntos.txt", sep="", header=TRUE)
rownames(image_data) <- image_data[,1]

image_selected_data <- image_data[maduration_samples,]
time_data <- c(0, 12, 20)

image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="MajorAxis", x_text=30, y_text=3.0, main="Eje mayor")
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="MinorAxis", x_text=25, y_text=11.3, main="Eje menor")
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="Area", x_text=25, y_text=0.80, main="Área")
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="Eccentricity", x_text=25, y_text=3.78, main="Excentricidad")
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="EquivDiameter", x_text=25, y_text=13.4, main="Eje mayor")
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="Perimeter", x_text=25, y_text=1.29, main="Perímetro")
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="Roundness", x_text=25, y_text=120, main="Redondez")
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="L.", x_text=25, y_text=6.1, main="L (claridad)")
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="desv_L.", x_text=25, y_text=11, main="Eje mayor")
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="a.", x_text=25, y_text=2.6)
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="desv_a.", x_text=25, y_text=15.6)
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="b.", x_text=25, y_text=4.6)
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="desv_b.", x_text=25, y_text=19)
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="C.ab", x_text=25, y_text=4.6, main="C (croma)")
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="desv_C.ab", x_text=25, y_text=56, main="Eje mayor")
image_variation_in_time(image_selected_data, time_data=time_data,
                        parameter="hab", x_text=30, y_text=8.5, main="H (tono)")
#image_variation_in_time(image_selected_data, time_data=time_data,
#                        parameter="desv_hab", x_text=25, y_text=6.9, main="Eje mayor")


#####################################
#####################################
####                             ####
#### CORRELACIÓN TEXTURA-IMAGEN  ####
####                             ####
#####################################
#####################################
## Hay que sustituir antes las comas por puntos. Vamos a crear un correlograma en el que
## se muestran todas las condiciones frente a la media de sus parámetros de textura e imagen.

## De image_data hay que eliminar LM12M
## De los archivos input hay que eliminar IP01X y LP04X
## De los siguientes grupos de muestras que configurar texture_data, seleccionar ejecutando
## uno de ellos (A, B, C ó D) y después continuar con el apartado ANÁLISIS común a cualquier
## grupo de muestras escogido.

## A. TODAS LAS MUESTRAS
texture_data <- c("LM01A", "LM02A", "LM03X", "LM04M",
                "LM05A", "LM06A", "LM07X", "LM08M", "LM09X", "LM10M", "LP01X", 
                "LP03X", "LP05X", "LP06X", "LP07X", "LP10X",
                "LP11X", "LP12X", "LP14X", "LP15X", "LP17X")

## B. MUESTRAS DE MADURACIÓN EN PLANTA
texture_data <- c("LM01A", "LM03X", "LM07X", "LM09X", "LP01X")

## C. MUESTRAS DE PASERA 1
texture_data <- c("LP01X", "LP03X", "LP05X", "LP06X", "LP07X")

## D. MUESTRAS DE PASERA 2
texture_data <- c("LP11X", "LP12X", "LP14X", "LP15X", "LP17X")

## ANÁLISIS

texture_data_frame <- data.frame(stringsAsFactors=FALSE)
for(i in 1:length(texture_data))
{
   current_data <- colMeans(read.csv(file=paste("output_", texture_data[i],".txt", sep=""), sep="", header=TRUE)[-1])
   texture_data_frame <- rbind(texture_data_frame, current_data)
   rownames(texture_data_frame)[i] <- texture_data[i]
}

data1 <- colMeans(read.csv(file=paste("output_", "LP01X", ".txt", sep=""), sep="", header=TRUE)[-1])
colnames(texture_data_frame) <- names(data1)

## Leer los datos de imagen de DigiEye.
image_data <- read.csv(file="imagen_semillas_puntos.txt", sep="", header=TRUE)
rownames(image_data) <- image_data[[1]]
image_data <- image_data[texture_data,]
image_numeric_data <- image_data[-c(1:5, 10, 13, 15, 16, 17, 18, 19, 21, 23, 24)]

## Leer datos de imagen hiperespectral.
spectrum_data <- read.csv("espectro_semillas_puntos.txt", sep="", header=TRUE)
rownames(spectrum_data) <- spectrum_data[[1]]
head(spectrum_data)
spectrum_data <- spectrum_data[texture_data,]
spectrum_numeric_data <- spectrum_data[-c(1:5)]

plot(as.numeric(spectrum_numeric_data["LM01A",]), type="l", ylim=c(0,0.75),
     ylab="Absorción", xlab="Longitud de onda (nm)")
spectrum_pca <- prcomp(spectrum_numeric_data)

## Aquí se muestra la importancia de las nuevas variables PCA definidas. Entre PCA1, PCA2
## y PCA3 se determina el 98% de la varianza total.
summary(spectrum_pca)
plot(spectrum_pca)
#biplot(spectrum_pca)
#spectrum_pca$rotation[,1]

full_data_frame1 <- cbind(texture_data_frame, image_numeric_data, spectrum_pca$x[,1:3])
mcor1 <- cor(full_data_frame1, method="pearson")
head(round(mcor1,2))
corr_coef <- rcorr(as.matrix(full_data_frame1), type="pearson")

cor.mtest <- function(mat, ...) {
   mat <- as.matrix(mat)
   n <- ncol(mat)
   p.mat<- matrix(NA, n, n)
   diag(p.mat) <- 0
   for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
         tmp <- cor.test(mat[, i], mat[, j], ...)
         p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
   }
   colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
   p.mat
}

p.mat <- round(cor.mtest(full_data_frame1),2)

#################################
## CORRELACIÓN TEXTURA-TEXTURA ##
#################################

## A este apartado no se le presta especial atención en el TFG porque se sale de los
## objettivos del estudio.
## Por defecto se crea la matriz de correlaciones segun el coeficiente de
## Pearson, para correlaciones lineales.

mcor1_text_text <- mcor1[1:7,1:7]
colnames(mcor1_text_text) <- c("Fuerza de rotura", "Fuerza de caída", "Energía de rotura", "Deformación de rotura",
                               "Pendiente", "Energía de compresión 50%", "Anchura de semilla")
rownames(mcor1_text_text) <- colnames(mcor1_text_text)
corrplot(mcor1_text_text, method="color", addgrid.col = "grey", addCoef.col = "black",
         tl.cex=0.8, tl.srt=45, tl.col="black", type="lower",
         p.mat = p.mat[1:7,1:7], sig.level = 0.05, insig = "blank",
         diag=TRUE, number.cex = 0.8)

################################
## CORRELACIÓN TEXTURA-IMAGEN ##
################################

mcor1_text_img <- mcor1[c(1,2,3,4,5,6,7, 17, 18, 19), 8:19]
dim(mcor1)

colnames(mcor1_text_img)[4] <- "Excentricidad"
colnames(mcor1_text_img)[5] <- "Perímetro"
colnames(mcor1_text_img)[6] <- "Redondez"
colnames(mcor1_text_img)[8] <- "C."
colnames(mcor1_text_img)[9] <- "H."

rownames(mcor1_text_img) <- c(colnames(mcor1_text_text), "PC1", "PC2", "PC3")
corrplot(mcor1_text_img, method="color", addgrid.col = "grey", addCoef.col = "black",
         tl.cex=0.8, tl.srt=45, tl.col="black", 
         p.mat = p.mat[c(1,2,3,4,5,6,7, 17, 18, 19), 8:19], sig.level = 0.05, insig = "blank",
         diag=TRUE, number.cex = 0.8)




## ANEXO
# http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram
## Para modificar como quieras todo el correlograma.
#
#corrplot(M, method="number")
#col<- colorRampPalette(c("red", "white", "blue"))(20)
#corrplot(M, type="lower", order="hclust", col=col, tl.cex=0.5,
#         tl.col="black", tl.srt=45, sig.level=0.1, insig="pch", pch=1)
