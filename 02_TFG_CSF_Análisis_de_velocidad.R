##################################################################
##                                                              ##
## Carlos Serrano Fajardo                                       ##
##                                                              ##
## Grado en Bioquímica                                          ##
## Dpto. Nutrición y Bromatología, Toxicología y Medicina Legal ##
## Grupo Color y Calidad de Alimentos                           ##
##                                                              ##
## ANÁLISIS TEXTUROMÉTRICO Y DE IMAGEN DE SEMILLAS DE UVA       ##
## TRABAJO DE FIN DE GRADO - SCRIPT II                          ##
##                                                              ##
## Esta obra está sujeta a la licencia                          ##
## Reconocimiento-NoComercial-CompartirIgual                    ##
## 4.0 Internacional de Creative Commons. Para                  ##
## ver una copia de esta licencia, visite                       ##
## http://creativecommons.org/licenses/by-nc-sa/4.0/.           ##
##                                                              ##
##################################################################

## En este script se realiza el análisis de los parámetros físicos obtenidos y se observa
## cómo varían en función de la velocidad del ensayo para escoger la mejor velocidad
## de movimiento de la sonda en la texturometría.

##########################
## PAQUETES Y FUNCIONES ##
##########################

library(stats)

###########################
## PARÁMETROS DE ENTRADA ##
###########################

## RUTA DEL ESPACIO DE TRABAJO (DONDE SE ENCUENTRA EL ARCHIVO .TAB DE ENTRADA)
setwd("~/Académico/Universidad/04_TFG/03_DATOS/parametros_textura")

## NOMBRE DEL ARCHIVO .TAB DE ENTRADA
data_vel1 <- read.csv(file="output_enteraseca_vel1.txt", sep="", header=TRUE)
data_vel2 <- read.csv(file="output_enteraseca_vel2.txt", sep="", header=TRUE)
data_vel3 <- read.csv(file="output_enteraseca_vel3.txt", sep="", header=TRUE)
data_vel4 <- read.csv(file="output_enteraseca_vel4.txt", sep="", header=TRUE)

## Para llevar a cabo la comparación de varios grupos de muestras entre sí se realiza por ANOVA.
## Previamente hay que asegurar la homogeneidad de varianzas (homoskedasticity), para ello se emplean
## los tests de Bartlett y Fligner-Killeen. Aquí se realiza el test de Barlett.
## El grado de libertad es n-1
anova_4comparisons <- function(input_vel1, input_vel2, input_vel3, input_vel4, alpha=0.05, df=3)
{
   all_data <- c(input_vel4,input_vel1,input_vel2,input_vel3)
   
   groups <- factor(c(rep("vel4", length(input_vel4)),
                      rep("vel1", length(input_vel1)),
                      rep("vel2", length(input_vel2)),
                      rep("vel3", length(input_vel3))))
   fit <- lm(formula = all_data ~ groups)
   ## Si p-value del test de Bartlett es mayor que 0.05, se rechaza que las varianzas
   ## difieran. Si chi-squared > Barlett's K squared, se acepta la hipótesis
   ## nula de que las varianzas son iguales.
   ## Se puede verificar la homoskedasticity también con el test de Fligner.
   ## Una vez verificada la homogeneidad se realiza el test de ANOVA
   ## La salida de la función es: Df (degree of freedom), Sum Sq (desviación(en los grupos
   ## y residual),Mean Sq (varianza(en los grupos y residual), F value (estadistico de Fisher
   ## de la comparación (variance within groups)/(variance residual)), Pr(>F) = p-value
   ## Puesto que p-value > 0.05, aceptamos la hipótesis nula H0.
   ## Se puede comparar tambien el F-value obtenido con el F-value tabulado.
   ## Si Tabulated F-value > computed F-value, aceptamos
   ## la hipótesis nula H0 > Los 4 grupos son iguales.
   return(list(
            Test_Bartlett=bartlett.test(all_data, groups),
            Chi_squared=qchisq((1-alpha), df),
            Test_ANOVA=anova(fit),
            F_Fisher=qf((1-alpha), 20, df)
         ))
}

## Esta función sirve para generar los gráficos para las 4 condiciones analizadas
barplot_4graph <- function(param_vel4, param_vel1, param_vel2, param_vel3,                        ylab ,main ,ylim=c(0,60))
{
   medias <- c(mean(param_vel4),
               mean(param_vel1),
               mean(param_vel2),
               mean(param_vel3))
   sds <-  c(sd(param_vel4),
             sd(param_vel1),
             sd(param_vel2),
             sd(param_vel3))
   
   names(medias) <- c("Fuerza rotura vel4", "Fuerza Rotura vel1","Fuerza Rotura vel2","Fuerza Rotura vel3")
   xpos <- barplot(medias, names.arg=c("0.2 mm/s", "0.5 mm/s","1 mm/s","2 mm/s"),
                   ylim=ylim, xlab="Velocidad de ensayo", col=c("tomato", "slateblue", "limegreen", "goldenrod"),
                   ylab=ylab, main=main, border="white")
   arrows(xpos, medias-sds, xpos, medias+sds, angle=90, code=3, length=0.1)
}

# T-STUDENT
# H0 : 0 = 1
# H1 : 0 < 1 ó 0 > 1
# La función t.test recibe como entrada dos vectores con los datos correspondientes y la hipótesis
# alternativa alt="less" o "greater". Por defecto asume un nivel de confianza del 95% que puede
# tunearse con el parámetro conf.level.


#################
## BREAK FORCE ##
#################

break_force_vel1 <- data_vel1[["break_force"]]
break_force_vel2 <- data_vel2[["break_force"]]
break_force_vel3 <- data_vel3[["break_force"]]
break_force_vel4 <- data_vel4[["break_force"]]
anova_4comparisons(break_force_vel1,break_force_vel2,break_force_vel3,break_force_vel4)

barplot_4graph(break_force_vel4, break_force_vel1, break_force_vel2,break_force_vel3,
               ylab="Fuerza (N)", main="Fuerza de rotura", ylim=c(0,65))

t.test(break_force_vel1,break_force_vel2,alt="greater")
t.test(break_force_vel1,break_force_vel3,alt="greater")
t.test(break_force_vel1,break_force_vel4,alt="greater")
t.test(break_force_vel2,break_force_vel3,alt="greater")
t.test(break_force_vel2,break_force_vel4,alt="greater")
t.test(break_force_vel3,break_force_vel4,alt="greater")

######################
## BREAK STEP FORCE ##
######################
break_step_force_vel1 <- data_vel1[["break_step_force"]]
break_step_force_vel2 <- data_vel2[["break_step_force"]]
break_step_force_vel3 <- data_vel3[["break_step_force"]]
break_step_force_vel4 <- data_vel4[["break_step_force"]]
anova_4comparisons(break_step_force_vel1,break_step_force_vel2,break_step_force_vel3,break_step_force_vel4)

barplot_4graph(break_step_force_vel4, break_step_force_vel1, break_step_force_vel2,break_step_force_vel3,
               ylab="Fuerza (N)", main="Fuerza de caída", ylim=c(0,65))

t.test(break_step_force_vel1,break_step_force_vel2,alt="greater")
t.test(break_step_force_vel1,break_step_force_vel3,alt="greater")
t.test(break_step_force_vel1,break_step_force_vel4,alt="greater")
t.test(break_step_force_vel2,break_step_force_vel3,alt="greater")
t.test(break_step_force_vel2,break_step_force_vel4,alt="greater")
t.test(break_step_force_vel3,break_step_force_vel4,alt="greater")

##################
## BREAK ENERGY ##
##################

break_energy_vel1 <- data_vel1[["break_energy"]]
break_energy_vel2 <- data_vel2[["break_energy"]]
break_energy_vel3 <- data_vel3[["break_energy"]]
break_energy_vel4 <- data_vel4[["break_energy"]]
anova_4comparisons(break_energy_vel1,break_energy_vel2,break_energy_vel3,break_energy_vel4)

barplot_4graph(break_energy_vel4, break_energy_vel1, break_energy_vel2,break_energy_vel3,
               ylab="Energía (mJ)", main="Energía de rotura", ylim=c(0,15))

t.test(break_energy_vel1,break_energy_vel2,alt="greater")
t.test(break_energy_vel1,break_energy_vel3,alt="greater")
t.test(break_energy_vel1,break_energy_vel4,alt="greater")
t.test(break_energy_vel2,break_energy_vel3,alt="greater")
t.test(break_energy_vel2,break_energy_vel4,alt="greater")
t.test(break_energy_vel3,break_energy_vel4,alt="greater")

#######################
## BREAK DEFORMATION ##
#######################

break_deformation_vel1 <- data_vel1[["break_deformation"]]
break_deformation_vel2 <- data_vel2[["break_deformation"]]
break_deformation_vel3 <- data_vel3[["break_deformation"]]
break_deformation_vel4 <- data_vel4[["break_deformation"]]
anova_4comparisons(break_deformation_vel1,break_deformation_vel2,break_deformation_vel3,break_deformation_vel4)

barplot_4graph(break_deformation_vel4, break_deformation_vel1, break_deformation_vel2, break_deformation_vel3,
               ylab="Deformación (%)", main="Deformación de rotura", ylim=c(0,0.3))

t.test(break_deformation_vel1,break_deformation_vel2,alt="greater")
t.test(break_deformation_vel1,break_deformation_vel3,alt="greater")
t.test(break_deformation_vel1,break_deformation_vel4,alt="greater")
t.test(break_deformation_vel2,break_deformation_vel3,alt="greater")
t.test(break_deformation_vel2,break_deformation_vel4,alt="greater")
t.test(break_deformation_vel3,break_deformation_vel4,alt="greater")

###########
## SLOPE ##
###########

## La pendiente de subida debería variar en función de la velocidad. Por
## normalizamos la misma dividiendo entre la velocidad de cada condición.
slope_vel1 <- data_vel1[["slope"]]
slope_vel2 <- data_vel2[["slope"]]
slope_vel3 <- data_vel3[["slope"]]
slope_vel4 <- data_vel4[["slope"]]
anova_4comparisons(slope_vel1,slope_vel2,slope_vel3,slope_vel4)

barplot_4graph(slope_vel4, slope_vel1, slope_vel2,slope_vel3,
               ylab="Pendiente (N/mm)", main="Pendiente de subida", ylim=c(0,250))

t.test(slope_vel1,slope_vel2,alt="greater")
t.test(slope_vel1,slope_vel3,alt="greater")
t.test(slope_vel1,slope_vel4,alt="greater")
t.test(slope_vel2,slope_vel3,alt="greater")
t.test(slope_vel2,slope_vel4,alt="greater")
t.test(slope_vel3,slope_vel4,alt="greater")

##################
## TOTAL ENERGY ##
##################

total_energy_vel1 <- data_vel1[["total_energy"]]
total_energy_vel2 <- data_vel2[["total_energy"]]
total_energy_vel3 <- data_vel3[["total_energy"]]
total_energy_vel4 <- data_vel4[["total_energy"]]
anova_4comparisons(total_energy_vel1,total_energy_vel2,total_energy_vel3,total_energy_vel4)

barplot_4graph(total_energy_vel4, total_energy_vel1, total_energy_vel2,total_energy_vel3,
               ylab="Energía (mJ)", main="Energía total de deformación", ylim=c(0,35))

t.test(total_energy_vel1,total_energy_vel2,alt="greater")
t.test(total_energy_vel1,total_energy_vel3,alt="greater")
t.test(total_energy_vel1,total_energy_vel4,alt="greater")
t.test(total_energy_vel2,total_energy_vel3,alt="greater")
t.test(total_energy_vel2,total_energy_vel4,alt="greater")
t.test(total_energy_vel3,total_energy_vel4,alt="greater")

#################### ANÁLISIS COMPLEMENTARIO ##########
##
## También se pueden evaluar los parámetros por parejas y representando gráficamente
## el resultado con la siguiente función.

barplot_2graph <- function(param_vel1, param_vel2, ylab ,main ,ylim=c(0,60))
{
   medias <- c(mean(param_vel1),
               mean(param_vel2))
   sds <-  c(sd(param_vel1),
             sd(param_vel2))
   
   names(medias) <- c("Fuerza rotura vel1", "Fuerza Rotura vel2")
   xpos <- barplot(medias, names.arg=c("0.5 mm/s","1 mm/s"),
                   ylim=ylim, xlab="Medias", col=c("tomato", "slateblue"),
                   ylab=ylab, main=main, border="white")
   arrows(xpos, medias-sds, xpos, medias+sds, angle=90, code=3, length=0.1)
}

t.test(break_force_vel1,break_force_vel2,alt="greater")
barplot_2graph(break_force_vel1, break_force_vel2,
               ylab="Fuerza (N)", main="Fuerza de rotura", ylim=c(0,65))

t.test(break_step_force_vel1,break_step_force_vel2,alt="greater")
barplot_2graph(break_step_force_vel1, break_step_force_vel2,
               ylab="Fuerza (N)", main="Fuerza de caída", ylim=c(0,65))

t.test(break_energy_vel1,break_energy_vel2,alt="greater")
barplot_2graph(break_energy_vel1, break_energy_vel2,
               ylab="Energía (mJ)", main="Energía de rotura", ylim=c(0,20))

t.test(break_deformation_vel1,break_deformation_vel2,alt="greater")
barplot_2graph(break_deformation_vel1, break_deformation_vel2,
               ylab="Deformación (%)", main="Deformación de rotura", ylim=c(0,0.3))

t.test(total_energy_vel1,total_energy_vel2,alt="greater")
barplot_2graph(total_energy_vel1, total_energy_vel2,
               ylab="Energía (mJ)", main="Energía total de deformación (mJ)", ylim=c(0,35))

t.test(slope_vel1,slope_vel2,alt="greater")
barplot_2graph(slope_vel1, slope_vel2,
               ylab="Pendiente (N/mm)", main="Pendiente de subida", ylim=c(0,250))


