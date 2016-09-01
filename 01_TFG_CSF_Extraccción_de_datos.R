##################################################################
##                                                              ##
## Carlos Serrano Fajardo                                       ##
##                                                              ##
## Grado en Bioquímica                                          ##
## Dpto. Nutrición y Bromatología, Toxicología y Medicina Legal ##
## Grupo Color y Calidad de Alimentos                           ##
##                                                              ##
## ANÁLISIS TEXTUROMÉTRICO Y DE IMAGEN DE SEMILLAS DE UVA       ##
## TRABAJO DE FIN DE GRADO - SCRIPT I                           ##
##                                                              ##    
## Esta obra está sujeta a la licencia                          ##
## Reconocimiento-NoComercial-CompartirIgual                    ##
## 4.0 Internacional de Creative Commons. Para                  ##
## ver una copia de esta licencia, visite                       ##
## http://creativecommons.org/licenses/by-nc-sa/4.0/.           ##
##                                                              ##
##################################################################

## En este código se realiza la extracción de los parámetros texturométricos de interés
## para cada una de las semillas procesadas.

##########################
## PAQUETES Y FUNCIONES ##
##########################

library(zoo) #Para ejecutar la función rollmean.
library(polynom) # Para resolver derivadas de un ajuste polinómico.
require(pracma) # Para hallar el área bajo la curva.

data_testing <- function(num_sample, raw, test_results_directory, output_name)
{
   for(k in 1:num_sample)
   {
      
      y <- as.numeric(as.vector(raw[,k*3-2])) # Fuerza (seg)
      z <- as.numeric(as.vector(raw[,k*3-1])) # Distancia (mm)
      x <- as.numeric(as.vector(raw[,k*3]))   # Tiempo (g)
      
      
      y2 <- rollmean(y, rollmean_threshold)
      
      x2 <- x[1:length(y2)]
      
      for(i in 1:length(y2))
      {
         if(y2[i] <= 20)
         {
            y2[i] <- 0
         }
      }
      
      local_max_positions <- find_peaks(y2, m=max_min_threshold)
      time_max <- x2[local_max_positions]
      force_max <- y2[local_max_positions]
      local_min_positions <- find_peaks(-y2, max_min_threshold)
      time_min <- x2[local_min_positions]
      force_min <- y2[local_min_positions]
      
      if((length(local_max_positions) + length(time_max)+ length(force_max)) >= 3 &&
         (length(local_min_positions) + length(time_min)+ length(force_min)) >= 3)
      {
         if(force_max[1] < fitting_shoulder_threshold && time_max[1] < 1) #si existe curva de encaje temprana
         {
            for(i in 1:local_min_positions[1]) #desde la posicion 1 hasta la del minimo de encaje
            {
               pend_encaje <- (y2[local_min_positions[1]]-y2[1])/(x2[local_min_positions[1]]-x2[1])
               
               ordenada <- -pend_encaje*x2[1]+y2[1]
               y2[i] <- pend_encaje*x2[i] + ordenada
            }
         }

      local_max_positions <- find_peaks(y2, max_min_threshold)
      time_max <- x2[local_max_positions]
      force_max <- y2[local_max_positions]
      local_min_positions <- find_peaks(-y2, max_min_threshold)
      time_min <- x2[local_min_positions]
      force_min <- y2[local_min_positions]
      
      fit3 <- lm(y2[1:local_max_positions[1]] ~ poly(x2[1:local_max_positions[1]], 3, raw=TRUE))
      
      pol3 <- function(a) fit3$coefficient[4]*a^3 + fit3$coefficient[3]*a^2 + fit3$coefficient[2]*a + fit3$coefficient[1]
      pol3_der1 <- function(a) fit3$coefficient[4]*3*a^2 + fit3$coefficient[3]*2*a + fit3$coefficient[2]
      pol3_der2 <- function(a) fit3$coefficient[4]*6*a + fit3$coefficient[3]*2
      
      pol3_der2_function <- polynomial(c(fit3$coefficient[3]*2, fit3$coefficient[4]*6))
      constant_pend <- solve(pol3_der2_function)
      
      pol3_der1_function <- polynomial(c(fit3$coefficient[2], fit3$coefficient[3]*2, fit3$coefficient[4]*3))
      pol3_der1_function
      pend_up <- predict(pol3_der1_function, constant_pend)
      
      current_plot <- paste("test_muestra_", k, ".jpeg", sep="")
      plot_directory <- paste(test_results_directory, current_plot, sep="")
      jpeg(file=plot_directory)
      plot(x2, y2, type="l", ylab="Force(g)", xlab="Time(s)", main=paste("Texturometría muestra", k))
      points(x2[local_max_positions], y2[local_max_positions], col="red", cex=1.5, pch=1)
      points(x2[local_min_positions], y2[local_min_positions], col="blue", cex=1.5, pch=1)
      lines(x2[1:local_max_positions[1]], predict(fit3), col="red")
      dev.off()
      } else{print(paste("La semilla", k, "no va: eliminar."))}
   }
}

find_peaks <- function (x, m)
{
   shape <- diff(sign(diff(x, na.pad = FALSE)))
   pks <- sapply(which(shape < 0), FUN = function(i){
      z <- i - m + 1
      z <- ifelse(z > 0, z, 1)
      w <- i + m + 1
      w <- ifelse(w < length(x), w, length(x))
      if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
   })
   pks <- unlist(pks)
   pks
}

texturometry_analysis <- function(num_sample, raw, wrong_samples, results_directory, output_name)
{
   data_sample <- numeric()
   data_break_force <- numeric()
   data_break_step_force <- numeric()
   data_break_energy <- numeric()
   data_break_deformation <- numeric()
   data_slope <- numeric()
   data_total_energy <- numeric()
   data_seed_width <- numeric()
   
   samples_to_analyze <- c(1:num_sample)
   if(length(wrong_samples>=1))
   {
      samples_to_analyze <- c(1:num_sample)[-wrong_samples]
   }
   
   for(k in samples_to_analyze)
   {
      print(paste("Procesando gráfica", k))
      
      y <- as.numeric(as.vector(raw[,k*3-2])) # Fuerza (seg)
      z <- as.numeric(as.vector(raw[,k*3-1])) # Distancia (mm)
      x <- as.numeric(as.vector(raw[,k*3]))   # Tiempo (g)
      
      ## Suavizamos la linea. Al hacerlo la curva se acorta y es inevitable.
      y2 <- rollmean(y, rollmean_threshold)
      
      ## Puesto que y se ha acortado, hay que redefinir x para que se ajuste a y2 en longitud
      ## A partir de aquí hay que trabajar siempre sobre "x2" e "y2".
      x2 <- x[1:length(y2)]
      
      ## Con esta instrucción se convierten en 0 todos los valores muy cercanos a 0. De esta
      ## forma se eliminan los máximos y mínimos locales producidos por muy pequeños montes
      ## en la zona final de la curva.
      for(i in 1:length(y2))
      {
         if(y2[i] <= 20)
         {
            y2[i] <- 0
         }
      }
      
      ## MAXIMOS Y MÍNIMOS LOCALES
      ## Ahora hallamos los máximos locales
      plot(x2,y2, type="l")
      ## Con esto hallamos los máximos y mínimos locales poniendo un umbral.
      local_max_positions <- find_peaks(y2, m=max_min_threshold)
      time_max <- x2[local_max_positions]
      force_max <- y2[local_max_positions]
      local_min_positions <- find_peaks(-y2, max_min_threshold)
      time_min <- x2[local_min_positions]
      force_min <- y2[local_min_positions]
      
      ## ELIMINAR CURVA(HOMBRO) DE ENCAJE
      if(force_max[1] < fitting_shoulder_threshold && time_max[1] < 1) #si existe curva de encaje temprana
      {
         for(i in 1:local_min_positions[1]) #desde la posicion 1 hasta la del minimo de encaje
         {
            pend_encaje <- (y2[local_min_positions[1]]-y2[1])/(x2[local_min_positions[1]]-x2[1])
            
            ordenada <- -pend_encaje*x2[1]+y2[1]
            y2[i] <- pend_encaje*x2[i] + ordenada
         }
      }
      
      ## REDEFINIR MÁXIMOS Y MÍNIMOS LOCALES
      local_max_positions <- find_peaks(y2, max_min_threshold)
      time_max <- x2[local_max_positions]
      force_max <- y2[local_max_positions]
      local_min_positions <- find_peaks(-y2, max_min_threshold)
      time_min <- x2[local_min_positions]
      force_min <- y2[local_min_positions]
      
      ######################
      ## FUERZA DE ROTURA ##
      ######################
      ## Se llamará break_force.
      ## La fuerza que se aplica sobre la semilla en el momento de romperse
      ## corresponde al primer máximo local
      break_mass <- force_max[1]
      
      ## Las unidades de break_force es g, pero lo queremos en newton.
      ## F(N)=F(kg*m/s^2)=m(kg)*a(9.80665m/s^2). Obtenemos la fuerza de rotura en Newton.
      break_force <- break_mass/1000*9.80665
      
      #########################
      ## BAJADA EN LA ROTURA ##
      #########################
      ## Se llamará break_step.
      ## Corresponde con la diferencia entre el primer máximo local
      ## y el primer mínimo local.
      break_step_mass <- force_max[1]-force_min[1]
      break_step_force <- break_step_mass/1000*9.80665
      
      #########################
      ## PENDIENTE DE SUBIDA ##
      #########################
      ## Con lm añadiendo el parametro "poly", se hace el ajuste a una funcion polinomica
      ## del grado que quiera. A mayor grado, mejor ajuste (R-squared). Con un ajuste
      ## a una polinomica de grado 3 da un buen resultado entre el primer punto y el
      ## maximo local.
      fit3 <- lm(y2[1:local_max_positions[1]] ~ poly(x2[1:local_max_positions[1]], 3, raw=TRUE))
      
      ## Cuando y=0 en la segunda derivada (pol3_der2), la pendiente de pol3 es constante.
      ## Para conocer la pendiente de subida más significativa, hay que hallar cuál
      ## es el máximo local de la primera derivada. Esto realmente no funciona, porque el
      ## ajuste de la funcion puede ser mucho más extenso que la curva de subida. Por esta razón
      ## es necesario buscar la pendiente de la funcion en el punto medio de la subida de la curva.
      ## Para resolver las funciones habrá que emplear el paquete polynom.
      ## Escribimos la funcion como un polinomio para poder resolver.
      pol3_function <- polynomial(c(fit3$coefficient[1], fit3$coefficient[2], fit3$coefficient[3], fit3$coefficient[4]))                              
      pol3_der1_function <- polynomial(c(fit3$coefficient[2], fit3$coefficient[3]*2, fit3$coefficient[4]*3))
      med_point <- x2[local_max_positions[1]]/2
      pend_up <- predict(pol3_der1_function, med_point)
      
      ## Esta pendiente está en unidades de g/s, hay que transformar las
      ## unidades a N/mm.
      ## F(N)= m(g)/1000*9.80665(m^2/s)
      ## espacio(mm)=velocidad(mm/s)/tiempo(s)
      slope <- pend_up/1000*9.80665/velocidad
      
      
      #######################
      ## ENERGÍA DE ROTURA ##
      #######################
      ## Para hallar la energía de rotura se pasa de tiempo a espacio multiplicando por
      ## la velocidad. Se pasa de masa en gramos a kg*m/s2. Para ello se pasa de g a kg y se
      ## 1 Julio = 1 Kg*m2/s2. Se obtiene multiplicando Fuerza (N) por metro.
      ## multiplica por la aceleración g. El espacio se deja en mm y se obtienen mJ.
      break_energy <- trapz(x2[1:local_max_positions[1]]*velocidad, y2[1:local_max_positions[1]]*9.80665/1000)
      break_energy ## Está en mJ (N/mm)
      
      
      #################################
      ## ENERGÍA TOTAL DE COMPRESIÓN ##
      #################################
      end_max <- length(local_max_positions)
      #plot(x2[1:local_max_positions[end_max]], y2[1:local_max_positions[end_max]], type="l")
      total_energy <- trapz(x2[1:local_max_positions[end_max]]*velocidad, y2[1:local_max_positions[end_max]]*9.8065/1000)
      total_energy ## Está en mJ (N/mm)
      
      ###########################
      ## ANCHURA DE LA SEMILLA ##
      ###########################
      ## Es la anchura de la semilla. Es el doble (en este caso) de lo que la máquina aplasta.
      ## Recordamos que la variable end_max guarda el último máximo, el punto de finalización.
      ## Sirve para hallar después la deformación de rotura.
      semi_width_time <- x2[local_max_positions[end_max]]
      semi_width_distance <- semi_width_time*velocidad
      seed_width <- semi_width_distance/deformation
      seed_width
      
      ###########################
      ## DEFORMACIÓN DE ROTURA ##
      ###########################
      ## Se trata del porcentaje de deformación al cual la semilla se rompe, es decir, ofrece
      ## su primer máximo local.
      break_time <- x2[local_max_positions[1]]
      break_distance <- break_time*velocidad
      break_deformation <- break_distance/seed_width
      break_deformation
      
      ## HEMOS OBTENIDO ESTOS DATOS que hay que almacenar:
      ## Gráfica
      current_plot <- paste("muestra_", k, ".jpeg", sep="")
      plot_directory <- paste(results_directory, current_plot, sep="")
      jpeg(file=plot_directory)
      plot(x2, y2, type="l", ylab="Force(g)", xlab="Time(s)", main=paste("Texturometría muestra", k))
      points(x2[local_max_positions[1]], y2[local_max_positions[1]], col="red", cex=1.2, pch=16)
      points(x2[local_max_positions[end_max]], y2[local_max_positions[end_max]], col="blue", cex=1.2, pch=16)
      lines(x2[1:local_max_positions[1]], predict(fit3), col="red")
      points(med_point, predict(pol3_function, med_point))
      dev.off()
      
      ## Número de muestra
      data_sample[k] <- k
      ## Fuerza de rotura (N)
      data_break_force[k] <- break_force
      ## Escalón de bajada (N)
      data_break_step_force[k] <- break_step_force
      ## Energía de rotura (N*mm=mJ)
      data_break_energy[k] <- break_energy
      ## Deformación de rotura (%)
      data_break_deformation[k] <- break_deformation
      ## Pendiente (N/mm)
      data_slope[k] <- slope
      ## Energía compresión 50% (N*mm=mJ)
      data_total_energy[k] <- total_energy
      ## Anchura de la semilla (mm)
      data_seed_width[k] <- seed_width
      print(paste("Completada la muestra", k))
   }
   
   output_data <- data.frame(sample= data_sample,
                             break_force= data_break_force,
                             break_step_force= data_break_step_force,
                             break_energy= data_break_energy,
                             break_deformation= data_break_deformation,
                             slope= data_slope,
                             total_energy= data_total_energy,
                             seed_width= data_seed_width,
                             stringsAsFactors=FALSE)
   
   output_data <- output_data[-wrong_samples,]
   
   write.table(output_data, file=paste(results_directory, output_name, sep=""))
   return(output_data)
}

###########################
## PARÁMETROS DE ENTRADA ##
###########################

## RUTA DEL ESPACIO DE TRABAJO (DONDE SE ENCUENTRA EL ARCHIVO .TAB DE ENTRADA)
setwd("~/Académico/Universidad/04_TFG/03_DATOS/input1")

## VELOCIDAD DE ENSAYO (mm/s)
velocidad <- 0.5 #vel1: 0.5 ; vel2: 1 ; vel3: 2 ; vel4: 0.2
## GRADO DE DEFORMACIÓN (%)
deformation <- 0.50
## UMBRAL DE SUAVIZADO MEDIANTE MEDIA MÓVIL
rollmean_threshold <- 20
## UMBRAL DE DETECCIÓN DE MÁXIMOS Y MÍNIMOS LOCALES
max_min_threshold <- 20
## UMBRAL DE DETECCIÓN DE LA CURVA/HOMBRO DE ENCAJE
fitting_shoulder_threshold <- 500

## NOMBRE DE LOS ARCHIVOS .TAB DE ENTRADA
all_data <- c("IA01X", "IA02X", "IA03X", "IM01A", "LM01A", "LM02A", "LM03X", "LM04M",
              "LM05A", "LM06A", "LM07X", "LM08M", "LM09X", "LM10M", "LM11X", "LP01X", 
              "LP02X", "LP03X", "LP05X", "LP06X", "LP07X", "LP08X", "LP09X", "LP10X",
              "LP11X", "LP12X", "LP13X", "LP14X", "LP15X", "LP16X", "LP17X")

file_extension <- ".TAB"
#raw <-  read.csv(file="enterafresca_vel2.TAB", sep="", header=TRUE)

data1 <- read.csv(file=paste("IA01X", file_extension, sep=""), sep="", header=TRUE)

##################################
## PREPROCESAMIENTO DE MUESTRAS ##
##################################

for(i in 1:length(all_data))
{

data_name <- all_data[i]
print(paste("Procesando muestra", data_name))
raw <- read.csv(file=paste(data_name, file_extension, sep=""), sep="", header=TRUE)
print("leído")
## NOMBRE DEL ARCHIVO .TXT DE SALIDA
output_name <- paste("output_",data_name, ".txt", sep="")

## DIRECTORIO DE SALIDA DE LAS GRÁFICAS PARA EL TESTEO PREVIO
test_results_directory <- paste("C:/Users/Carlos-1/Documents/Académico/Universidad/04_TFG/03_DATOS/test_", data_name, "/", sep="")
dir.create(test_results_directory)

## DIRECTORIO DE SALIDA DE LOS RESULTADOS REALES
results_directory <- paste("C:/Users/Carlos-1/Documents/Académico/Universidad/04_TFG/03_DATOS/results_", data_name, "/", sep="")
dir.create(results_directory)

## NÚMERO DE PARÁMETROS PARA CADA MUESTRA (normalmente son 3: tiempo, fuerza, distancia)
num_param <- 3

## NÚMERO DE MUESTRAS O RÉPLICAS DE LA MISMA CONDICIÓN
num_sample <- ncol(raw)/num_param
# num_sample <- 50 podrían introducirse manualmente si resultara en algún error


################################
## PRE-PROCESAMIENTO DE DATOS ##
################################
   
## RENOMBRAR LA MATRIZ
## num_param es el numero de parametros que se tiene para cada semilla (force, distance, time)
## Esto funciona para una dataframe que tenga tantas filas como puntos y sus columnas sean: parametro1.1,
## parametro2.1, parametro3.1, parametro1.2, parametro2.2, etc.
# Hay 3 parámetros para cada semilla (tiempo, distancia y fuerza). En esta prueba analizamos 5
# muestras.
## A continuacion colocamos los nombres adecuados. Para ellos hay que poner "muestra", "parametro",
## "unidad".
print("renombrando datos")
names <- as.vector(as.matrix(raw)[1,])
colnames(raw) <- names

# Ponemos las unidades
for (i in 1:length(names))
{
   if (names[i] == "Force")
   {
      names[i] <- paste(names[i], "(g)", sep="")
   }
   
   if (names[i] == "Distance")
   {
      names[i] <- paste(names[i], "(mm)", sep="")
   }
   
   if (names[i] == "Time")
   {
      names[i] <- paste(names[i], "(sec)", sep="")
   }
}

# Ponemos las muestras
for(i in 1:num_sample)
{
   names[((i-1)*3)+1] <- paste("Sample_", paste(i, names[((i-1)*3)+1]), sep="")
   names[((i-1)*3)+2] <- paste("Sample_", paste(i, names[((i-1)*3)+2]), sep="")
   names[((i-1)*3)+3] <- paste("Sample_", paste(i, names[((i-1)*3)+3]), sep="")
}
names(raw) <- names

## Ahora hay que eliminar las filas 1 y 2 desplazando todo hacia arriba.
raw <- raw[-(1:2),]
rownames(raw) <- (1:nrow(raw))   

## Si el archivo es muy grande teniendo una gran resolución, entonces hay que
## simplificar los datos.
print("redimensionando datos")

sample <- seq(from=2, to=length(raw[,1]), by=2)
if(nrow(raw) > 1000 && nrow(raw)<1500)
{
   raw <- raw[-sample,]
}
if(nrow(raw) >= 1500 && nrow(raw)<3000)
{
   raw <- raw[-sample,]
   raw <- raw[-sample,]
}
if(nrow(raw) >= 3000)
{
   raw <- raw[-sample,]
   raw <- raw[-sample,]
   raw <- raw[-sample,]
}

## RELLENAR HUECOS
## Con este bucle rellenamos los huecos "" con el numero cero, de esta forma todas las columnas
## tendran la misma longitud.
print("rellenando huecos")
raw <- data.frame(lapply(raw, as.character), stringsAsFactors=FALSE)
head(raw)
for (i in 1:nrow(raw))
{
   for(j in 1:ncol(raw))
   {
      if(as.vector(raw[i,j])=="" || is.na(raw[i,j]) || is.null(raw[i,j]))
      {
         raw[i,j] <- 0
      }
   }
}

## HACER 0 LOS VALORES NEGATIVOS
## Con la siguiente función vamos a convertir en 0 los valores negativos.
print("haciendo 0 los valores negativos")
for (i in 1:nrow(raw))
{
   for(j in 1:ncol(raw))
   {
      if(as.vector(raw[i,j]) < "0" || is.na(raw[i,j]))
      {
         raw[i,j] <- 0
      }
   }
}
head(raw)
plot(raw[,3], raw[,1], type="l")
is.list(raw)
raw[1]
raw[2]
mean(as.vector(raw[2,]))


new_raw <- matrix(0, ncol=3, nrow=nrow(raw))


for(i in 1:ncol(raw))
{
   new_raw[i] <- mean(as.numeric(as.matrix(raw)[i,]))
   
}


rowMeans(as.matrix(raw))


###################################
## PURGA DE GRÁFICAS DEFECTUOSAS ##
###################################
print("haciendo data_testing")
data_testing(num_sample, raw, test_results_directory, output_name)

}

## VER LAS GRÁFICAS GENERADAS EN LA NUEVA CARPETA "test_..." CREADA.
## INTRODUCIR AQUÍ (dentro de c(n1, n2, ...)) LOS NÚMEROS DE LAS MUESTRAS QUE SE DESEAN ELIMINAR Y EJECUTAR
## en esta evaluación manual se ha visto que hay muestras que es mejor descartar por su mala calidad.
wrong_samples_list <- list(IA01X=c(2,5,8,11,12,13,17,18,19,24,25,27,28,31,35,38,46,48,49,50,51,53,54,55,56,58,59,60,61,62,63,65,67,68,69,70,73,79,80), ## MEJOR QUITAR LAS IA_
                           IA02X=c(2,4,5,6,7,8,9,10,12,13,14,16,17,20,22,23,24,25,26,27,29,30,31,32,35,36,37,41,42,43,44, 46,47,48,51,52,53,54,55,56,57,58,59,60,61,62,63,64,66,68,69), ##MEJOR ELIMINARLA TODA
                           IA03X=c(57), ## QUITARLO
                           IM01A=c(3,5, 6,7, 8, 9,10,11, 14,13,17,18,25), ## SE PUEDE QUITAR, NO TIENE OTRAS DE IM_.
                           LM01A=c(52),
                           LM02A=c(1,3,4,6,10,12,15,17,18,31,32,39,41,43,49,53,55),
                           LM03X=c(4,20,21,23,26,39,40), # OK
                           LM04M=c(5,11, 12,13,15,17,20,22,24,30,31, 32,36,37,43,44,56,57,60,64,69,73,74),
                           LM05A=c(24), # OK
                           LM06A=c(6,8,12, 15,19,24,27,31,34,35,37,50,54,55), # OK
                           LM07X=c(11, 38, 41), # OK
                           LM08M=c( 4,19,56,66,67), # OK
                           LM09X=c(2,3,4,10, 11,14,18,25,31,41,43,49), # OK,
                           LM10M=c(2,9,16,33,34,41,43,44,77,107,114), # OK
                           LM11X=c(), # SÓLO TIENE 1 MUESTRA: DESCARTAR
                           LP01X=c(2,5,6,8,20,42,52), # OK
                           LP02X=c(), ## ES DEL MISMO DÍA QUE LP01X: DESCARTAR
                           LP03X=c(9,25,40,50), # OK
                           LP05X=c(1,22,26,27,30,33,39,43,46,58),
                           LP06X=c(3,12,15,18,29,32,33,35,44,48,52),
                           LP07X=c(12,16,20,24,31,33,39,43, 46,53,63,64),
                           LP08X=c(47), # UN DESASTRE: DESCARTAR
                           LP09X=c(1,2,4,9,15,19,22,23,26,32,35,36,49,51,52, 56,57,60,63,64, 67,68,69,75,77,78,80,86, 88,89,90,91,92,97,98,100,101), #POCO FIABLE
                           LP10X=c( 3,9,22,30,36,41,48,50,51,54, 57,62),
                           LP11X=c(1,2,6,26,32,46,47,58),
                           LP12X=c(1,5,6,7,9,10,12,17,19,29,30,31,33,34,35, 37,38,43,45,48,51,52,53),
                           LP13X=c(), # PINTA FATAL, NI SELECCIONAMOS
                           LP14X=c(8,19,42),
                           LP15X=c( 1,3,4,5,6,10,25,31,32,33,34,35,36,38,39,40,42,47,53,66, 67),
                           LP16X=c(), # SOLO HAY UNA GRÁFICA: SE DESCARTA
                           LP17X=c(3,4,14,16,17,19,21, 27,31,33,36, 37,39,43,44,54,55,56,57,58,60))

## Pero ahora va a ser esto:
clean_data <- c("LM01A", "LM02A", "LM03X", "LM04M",
              "LM05A", "LM06A", "LM07X", "LM08M", "LM09X", "LM10M", "LP01X", 
              "LP03X", "LP05X", "LP06X", "LP07X", "LP10X",
              "LP11X", "LP12X", "LP14X", "LP15X", "LP17X")
## Puesto que hemos eliminado datos de entrada, también hay que eliminar sus erróneos:
## Quitamos: IA01X, IA02X, IA03X, IM01A, LM11X, LP02X, LP08X, LP09X, LP13X, LP16X.
wrong_samples_list <- list(LM01A=c(52),
                           LM02A=c(1,3,4,6,10,12,15,17,18,31,32,39,41,43,49,53,55),
                           LM03X=c(4,20,21,23,26,39,40), # OK
                           LM04M=c(5,11, 12,13,15,17,20,22,24,30,31, 32,36,37,43,44,56,57,60,64,69,73,74),
                           LM05A=c(24), # OK
                           LM06A=c(6,8,12, 15,19,24,27,31,34,35,37,50,54,55), # OK
                           LM07X=c(11, 38, 41), # OK
                           LM08M=c( 4,19,56,66,67), # OK
                           LM09X=c(2,3,4,10, 11,14,18,25,31,41,43,49), # OK,
                           LM10M=c(2,9,16,33,34,41,43,44,77,107,114), # OK
                           LP01X=c(2,5,6,8,20,42,52), # OK
                           LP03X=c(9,25,40,50), # OK
                           LP05X=c(1,22,26,27,30,33,39,43,46,58),
                           LP06X=c(3,12,15,18,29,32,33,35,44,48,52),
                           LP07X=c(12,16,20,24,31,33,39,43, 46,53,63,64),
                           LP10X=c( 3,9,22,30,36,41,48,50,51,54, 57,62),
                           LP11X=c(1,2,6,26,32,46,47,58),
                           LP12X=c(1,5,6,7,9,10,12,17,19,29,30,31,33,34,35, 37,38,43,45,48,51,52,53),
                           LP14X=c(8,19,42),
                           LP15X=c( 1,3,4,5,6,10,25,31,32,33,34,35,36,38,39,40,42,47,53,66, 67),
                           LP17X=c(3,4,14,16,17,19,21, 27,31,33,36, 37,39,43,44,54,55,56,57,58,60))


############################
## PROCESAMIENTO DE DATOS ##
############################

for(i in 1:length(clean_data))
{
   data_name <- clean_data[i]
   print(paste("Procesando muestra", data_name))
   raw <- read.csv(file=paste(data_name, file_extension, sep=""), sep="", header=TRUE)

   ## NOMBRE DEL ARCHIVO .TXT DE SALIDA
   output_name <- paste("output_",data_name, ".txt", sep="")
   
   ## DIRECTORIO DE SALIDA DE LAS GRÁFICAS PARA EL TESTEO PREVIO
   test_results_directory <- paste("C:/Users/Carlos-1/Documents/Académico/Universidad/04_TFG/03_DATOS/test_", data_name, "/", sep="")
   
   ## DIRECTORIO DE SALIDA DE LOS RESULTADOS REALES
   results_directory <- paste("C:/Users/Carlos-1/Documents/Académico/Universidad/04_TFG/03_DATOS/results_", data_name, "/", sep="")

   ## NÚMERO DE PARÁMETROS PARA CADA MUESTRA (normalmente son 3: tiempo, fuerza, distancia)
   num_param <- 3
   
   ## NÚMERO DE MUESTRAS O RÉPLICAS DE LA MISMA CONDICIÓN
   num_sample <- ncol(raw)/num_param
   # num_sample <- 50 podrían introducirse manualmente si resultara en algún error
   
   ################################
   ## PRE-PROCESAMIENTO DE DATOS ##
   ################################
   
   ## RENOMBRAR LA MATRIZ
   ## num_param es el numero de parametros que se tiene para cada semilla (force, distance, time)
   ## Esto funciona para una dataframe que tenga tantas filas como puntos y sus columnas sean: parametro1.1,
   ## parametro2.1, parametro3.1, parametro1.2, parametro2.2, etc.
   # Hay 3 parámetros para cada semilla (tiempo, distancia y fuerza). En esta prueba analizamos 5
   # muestras.
   ## A continuacion colocamos los nombres adecuados. Para ellos hay que poner "muestra", "parametro",
   ## "unidad".

   names <- as.vector(as.matrix(raw)[1,])
   colnames(raw) <- names
   
   # Ponemos las unidades
   for (i in 1:length(names))
   {
      if (names[i] == "Force")
      {
         names[i] <- paste(names[i], "(g)", sep="")
      }
      
      if (names[i] == "Distance")
      {
         names[i] <- paste(names[i], "(mm)", sep="")
      }
      
      if (names[i] == "Time")
      {
         names[i] <- paste(names[i], "(sec)", sep="")
      }
   }
   
   # Ponemos las muestras
   for(i in 1:num_sample)
   {
      names[((i-1)*3)+1] <- paste("Sample_", paste(i, names[((i-1)*3)+1]), sep="")
      names[((i-1)*3)+2] <- paste("Sample_", paste(i, names[((i-1)*3)+2]), sep="")
      names[((i-1)*3)+3] <- paste("Sample_", paste(i, names[((i-1)*3)+3]), sep="")
   }
   names(raw) <- names
   
   ## Ahora hay que eliminar las filas 1 y 2 desplazando todo hacia arriba.
   raw <- raw[-(1:2),]
   rownames(raw) <- (1:nrow(raw))   
   
   ## Si el archivo es muy grande teniendo una gran resolución, entonces hay que
   ## simplificar los datos.

   sample <- seq(from=2, to=length(raw[,1]), by=2)
   if(nrow(raw) > 1000 && nrow(raw)<1500)
   {
      raw <- raw[-sample,]
   }
   if(nrow(raw) >= 1500 && nrow(raw)<3000)
   {
      raw <- raw[-sample,]
      raw <- raw[-sample,]
   }
   if(nrow(raw) >= 3000)
   {
      raw <- raw[-sample,]
      raw <- raw[-sample,]
      raw <- raw[-sample,]
   }
   
   ## RELLENAR HUECOS
   ## Con este bucle rellenamos los huecos "" con el numero cero, de esta forma todas las columnas
   ## tendran la misma longitud.
   print("rellenando huecos")
   raw <- data.frame(lapply(raw, as.character), stringsAsFactors=FALSE)
   
   for (i in 1:nrow(raw))
   {
      for(j in 1:ncol(raw))
      {
         if(as.vector(raw[i,j])=="" || is.na(raw[i,j]) || is.null(raw[i,j]))
         {
            raw[i,j] <- 0
         }
      }
   }
   
   ## HACER 0 LOS VALORES NEGATIVOS
   ## Con la siguiente función vamos a convertir en 0 los valores negativos.
   print("haciendo 0 los valores negativos")
   for (i in 1:nrow(raw))
   {
      for(j in 1:ncol(raw))
      {
         if(as.vector(raw[i,j]) < "0" || is.na(raw[i,j]))
         {
            raw[i,j] <- 0
         }
      }
   }
   
   ###################################
   ## PURGA DE GRÁFICAS DEFECTUOSAS ##
   ###################################
   texturometry_analysis(num_sample=num_sample, raw=raw, wrong_samples=wrong_samples_list[[data_name]],
                         results_directory=results_directory, output_name=output_name)
   
}


## En este punto ya hemos obtenido los parámetros físicos de cada una de las muestras que
## hemoso analizado. Estos son los parámetros extraídos con los que se trabajará posteriormente:
##    Fuerza de rotura (N): break_force
##    Escalón de bajada (N): break_step_force
##    Energía de rotura (N*mm=mJ): break_energy
##    Deformación de rotura (%): break_deformation
##    Pendiente (N/mm): slope
##    Energía compresión 50% (N*mm=mJ): total_energy
##    Anchura de la semilla (mm): seed_width


############ PARA TESTEAR SÓLO MUESTRAS CONCRETAS ##############
k=1
      y <- as.numeric(as.vector(raw[,k*3-2])) # Fuerza (seg)
      z <- as.numeric(as.vector(raw[,k*3-1])) # Distancia (mm)
      x <- as.numeric(as.vector(raw[,k*3]))   # Tiempo (g)
      plot(x,y,ylab="Fuerza (g)", xlab="Tiempo (s)", type="l")
      y2 <- rollmean(y, rollmean_threshold)
      x2 <- x[1:length(y2)]
      for(i in 1:length(y2))
      {
         if(y2[i] <= 20)
         {
            y2[i] <- 0
         }
      }
      plot(x2,y2, ylab="Fuerza (g)", xlab="Tiempo (s)", type="l")
      
      local_max_positions <- find_peaks(y2, m=max_min_threshold)
      time_max <- x2[local_max_positions]
      force_max <- y2[local_max_positions]
      local_min_positions <- find_peaks(-y2, max_min_threshold)
      time_min <- x2[local_min_positions]
      force_min <- y2[local_min_positions]
      ## ELIMINAR CURVA(HOMBRO) DE ENCAJE
      if(force_max[1] < fitting_shoulder_threshold && time_max[1] < 1) #si existe curva de encaje
      {
         print("Hay curva de encaje!")
         for(i in 1:local_min_positions[1]) #desde la posicion 1 hasta la del minimo de encaje
         {
            pend_encaje <- (y2[local_min_positions[1]]-y2[1])/(x2[local_min_positions[1]]-x2[1])
            
            ordenada <- -pend_encaje*x2[1]+y2[1]
            y2[i] <- pend_encaje*x2[i] + ordenada
         }
      }
      plot(x2, y2, ylab="Fuerza (g)", xlab="Tiempo (s)", type="l")
      ## REDEFINIR MÁXIMOS Y MÍNIMOS LOCALES
      local_max_positions <- find_peaks(y2, max_min_threshold)
      time_max <- x2[local_max_positions]
      force_max <- y2[local_max_positions]
      local_min_positions <- find_peaks(-y2, max_min_threshold)
      time_min <- x2[local_min_positions]
      force_min <- y2[local_min_positions]
      ######################
      ## FUERZA DE ROTURA ##
      ######################
      break_mass <- force_max[1]
      break_force <- break_mass/1000*9.80665
      #########################
      ## BAJADA EN LA ROTURA ##
      #########################
      break_step_mass <- force_max[1]-force_min[1]
      break_step_force <- break_step_mass/1000*9.80665
      #########################
      ## PENDIENTE DE SUBIDA ##
      #########################
      fit3 <- lm(y2[1:local_max_positions[1]] ~ poly(x2[1:local_max_positions[1]], 3, raw=TRUE))
      pol3_function <- polynomial(c(fit3$coefficient[1], fit3$coefficient[2], fit3$coefficient[3], fit3$coefficient[4]))                              
      pol3_der1_function <- polynomial(c(fit3$coefficient[2], fit3$coefficient[3]*2, fit3$coefficient[4]*3))

      #pend_up <- predict(pol3_der1_function, constant_pend)
      med_point <- x2[local_max_positions[1]]/2
      pend_up <- predict(pol3_der1_function, med_point)
      ## pend_up está en g/s. Ahora vamos a pasarlo a N/mm dividiendo entre la velocidad
      ## y multiplicando por 1000 (de g a kg) y por la aceleracion g (pasar a N)
      slope <- pend_up/1000*9.80665/velocidad
      
      #######################
      ## ENERGÍA DE ROTURA ##
      #######################
      ## Para hallar la energía de rotura se pasa de tiempo a espacio multiplicando por
      ## la velocidad. Se pasa de masa en gramos a kg*m/s2. Para ello se pasa de g a kg y se
      ## 1 Julio = 1 Kg*m2/s2. Se obtiene multiplicando Fuerza (N) por metro.
      ## multiplica por la aceleración g. El espacio se deja en mm y se obtienen mJ.
      break_energy <- trapz(x2[1:local_max_positions[1]]*velocidad, y2[1:local_max_positions[1]]*9.80665/1000)
      break_energy ## Está en mJ (N/mm)

      ###########################
      ## ENERGÍA DE COMPRESIÓN ##
      ###########################
      end_max <- length(local_max_positions)
      #plot(x2[1:local_max_positions[end_max]], y2[1:local_max_positions[end_max]], type="l")
      total_area <- trapz(x2[1:local_max_positions[end_max]]*velocidad, y2[1:local_max_positions[end_max]]*9.8065/1000)
      total_energy ## Está en mJ (N/mm)
      
      ###########################
      ## ANCHURA DE LA SEMILLA ##
      ###########################
      semi_width_time <- x2[local_max_positions[end_max]]
      semi_width_distance <- semi_width_time*velocidad
      seed_width <- semi_width_distance/deformation
      seed_width
      
      ###########################
      ## DEFORMACIÓN DE ROTURA ##
      ###########################
      break_time <- x2[local_max_positions[1]]
      break_distance <- break_time*velocidad
      break_deformation <- break_distance/seed_width
      break_deformation
      
      ## HEMOS OBTENIDO ESTOS DATOS que hay que almacenar:
      ## Gráfica
      current_plot <- paste("muestra_", k, ".jpeg", sep="")
      plot_directory <- paste(results_directory, current_plot, sep="")
      jpeg(file=plot_directory)
      plot(x2, y2, type="l", ylab="Force(g)", xlab="Time(s)", main=paste("Texturometría muestra", k))
      points(x2[local_max_positions[1]], y2[local_max_positions[1]], col="red", cex=0.5, pch=16)
      points(x2[local_max_positions[end_max]], y2[local_max_positions[end_max]], col="blue", cex=0.75, pch=16)
      lines(x2[1:local_max_positions[1]], predict(fit3), col="red")
      dev.off()
      
      ## Número de muestra
      data_sample[k] <- k
      ## Fuerza de rotura (N)
      data_break_force[k] <- break_force
      ## Escalón de bajada (N)
      data_break_step_force[k] <- break_step_force
      ## Energía de rotura (N*mm=mJ)
      data_break_energy[k] <- break_energy
      ## Deformación de rotura (%)
      data_break_deformation[k] <- break_deformation
      ## Pendiente (N/mm)
      data_slope[k] <- slope
      ## Energía compresión 50% (N*mm=mJ)
      data_total_energy[k] <- total_energy
      ## Anchura de la semilla (mm)
      data_seed_width[k] <- seed_width
      print(paste("Completada la muestra", k))
   }
   
   output_data <- data.frame(sample= data_sample,
                             break_force= data_break_force,
                             break_step_force= data_break_step_force,
                             break_energy= data_break_energy,
                             break_deformation= data_break_deformation,
                             slope= data_slope,
                             total_energy= data_total_energy,
                             seed_width= data_seed_width,
                             stringsAsFactors=FALSE)
   
   write.table(output_data, file=paste(results_directory, output_name, sep=""))
   return(output_data)
}






