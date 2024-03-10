#TAKEHOME PRIMERA PARTE SIMULACIÓN 
#JUAN ESTEBAN GARZÓN BORDA-SARA MORALES CARRILLO
#Preámbulo
options(scipen = 999)#Para quitar notación científica
#Archivo con las funciones previamente elaboradas
source(file.choose())
#Librerias 
#Para los autocorrelogramas
library(ggplot2)

#ARMA (p, q)----
#PROCESO 1----
# Coeficientes AR
phi <- c(0, -0.2)
phi1 <- phi[2]
phi2 <- phi[1]
p <- length(phi)
#Verificar estacionariedad
rootAR <- Re(polyroot(c(1, -phi1, -phi2)))
rootAR
# Coeficientes MA
theta <- c(0.5)
theta1 <- theta[1]
q <- length(theta)
#Verificar invertibilidad
rootMA <- Re(polyroot(c(1, -theta1)))
rootMA
#Establecer semilla
set.seed(123)  
# Longitud de la serie de tiempo
n <- 7000        
#Ruido blanco de la serie
at <- rnorm(n, mean = 0)
# Iniciar la serie temporal para un modelo ARMA (p, q) el código quedó flexible en los órdenes, se definen arriba y código es capaz de modelo cualquier arma, para el ejercicio se escogió un ARMA (2;1)
arma_series <- numeric(n)  

for (i in (max(p, q) + 1):n) {
   parte_AR <- sum(phi * arma_series[(i - p):(i - 1)])
   parte_MA <- sum(theta * at[(i - q):(i - 1)])
   arma_series[i] <- parte_AR + parte_MA + at[i]
}

#Estimación del arma por comando directo

arima(arma_series, c(p, 0, q), fixed = c(NA, 0, NA, NA))

#Gráfica proceso 1
m <- 6000
arma_series1 <- arma_series[(m+1):n]
p1 <- plot(arma_series1, type = "l", col = "#1874CD", lwd=1.5, main = "Proceso 1", ylab = "Simulación", xlab = "Tiempo", sub = paste("ARMA (", p, ";" ,q, ")")) 

#Autocorrelogramas
autocorrelogramas(arma_series1,30,pacf_(rho_k(arma_series1,30)),acf_(arma_series1,30),q)

#Prueba de  Box-Pierce
bp(arma_series1,2,1,30)

#Comprobación con comando directo
Box.test(arma_series1,30, "Box-Pierce")

#PROCESO 2----
#Defino la constante 
D <- 10
arma_series2 <- arma_series1 + D 

df_proc12 <- data.frame(arma_series1, arma_series2) #Se crea un data-frame para poder poner las gráficas enla misma salida y con los mismos ejes

#GRÁFICAS
#Gráfica del proceso 2
p2 <- plot(arma_series2, type = "l", main = "Proceso 2", col = "#6C2AFA", lwd = 1.5, ylab = "Simulación", xlab = "Tiempo", sub = paste("ARMA (", p, ";", q, ")"))
colores <- hcl.colors(10, "Berlin")
#Gráfica comparando proceso 1 con proceso 2
p12 <- matplot(df_proc12, type="l", lty = 1, lwd = 1.5, col = colores, ylab = "Simulación", xlab = "Tiempo", sub = paste("ARMA (", p, ";" ,q, ")"))+ title(main = "Comparación Proceso 1 y Proceso 2")

#Autocorrelogramas
autocorrelogramas(arma_series2, length(arma_series2)/4, pacf_(rho_k(arma_series2, length(arma_series2)/4)), acf_(arma_series2, length(arma_series2)/4), length(arma_series2)/4)

#Prueba de Box-Pierce
bp(arma_series2,2,1,40)
#Comprobación con comando directo
Box.test(arma_series2,40,"Box-Pierce")


#PROCESO 3----
#Definir parámetros
delta <- 0.009
tiempo <- seq(from = 1,to = n, by = 1)
tendencia <- delta * tiempo
arma_series3 <- arma_series2 + tendencia 
#Gráfica proceso 3
p3 <- plot(arma_series3[6000:7000], type = "l", main = "Proceso 3", col = "#FA9FB5", lwd = 1.5, ylab = "Simulación", xlab = "Tiempo", sub = paste("ARMA (", p, ";", q, ")")) 
#Data frame para graficar las 3 series 
df_proc3 <- data.frame(arma_series[1001:2000], arma_series2, arma_series3[1001:2000])
#Vector colores para las gráficas
colores <- hcl.colors(10, "Berlin")
#Gráfica con el proceso 1, proceso 2, proceso 3
p123 <- matplot(df_proc3, type = "l", col = colores, lty = 1, lw = 1.5, ylab = "Simulación", xlab = "Tiempo") + title(main = "Proceso 1, Proceso 2, Proceso 3" )

#Autocorrelogramas
T3 <- 100
autocorrelogramas(arma_series3, T3/4, pacf_(rho_k(arma_series3, T3/4)), acf_(arma_series3, T3/4), length(arma_series3)/4)

#Prueba de  Box-Pierce
bp(arma_series3,2,1,30)

#Comprobación con comando directo
Box.test(arma_series3,30, "Box-Pierce")


#PROCESO 4 con dos estacionalidades----

comp_estacional1 <- rep(c(6,0,0,0,0,0,0,0,0,0), length.out=n)
comp_estacional2 <- rep(c(1,2,3,4,5,6,7,8,9,10), length.out=n)
#Proceso 4.1
arma_series4_1 <- arma_series1 + comp_estacional1

#Proceso 4.2
arma_series4_2 <- arma_series1 + comp_estacional2

#Gráfica proceso 4.1
p41 <- plot(arma_series4_1[6001:7000], type = "l", main = "Proceso 4.1", col = "#FA8072", lwd = 1.5, xlab = "Tiempo", ylab = "Simulación", sub = paste("ARMA (", p, ";", q, ")"))

#Gráfica proceso 4.2
p42 <- plot(arma_series4_2[6001:7000], type = "l", main = "Proceso 4.2", col = "#43CD80", lwd = 1.5, ylab="Simulación",xlab = "Tiempo")
#Data frame proceso 1 y proceso 4.1
df141 <- data.frame(arma_series1, arma_series4_1[6901:7000])
#Vector colores
colores1 <- hcl.colors(2, "Berlin")

#Gráfica comparación Proceso 1 y Proceso 4.1
p141 <- matplot(df141, type = "l", col = colores1, main = "Proceso 1 VS Proceso 4.1", lty = 1, ylab = "Simulación", xlab="Tiempo")

#Dataframe para el proceso 1 y proceso 4.2
df142 <- data.frame(arma_series1, arma_series4_2[6901:7000])
#Gráfica comparación Proceso 1 y Proceso 4.2
p142 <- matplot(df142, type = "l", col = colores1, main = "Proceso 1 VS Proceso 4.2", lty = 1, xlab = "Tiempo", ylab = "Simulación")

#Prueba de  Box-Pierce
bp(arma_series4_1,2,1,30)
bp(arma_series4_2,2,1,30)
#Comprobación por comando direcro
Box.test(arma_series4_1,30, "Box-Pierce")
Box.test(arma_series4_2,30, "Box-Pierce")
#-

