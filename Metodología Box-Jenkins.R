#Takehome Macroeconometría
####Segundo Punto----
#Preambulo----
#funciones previamemte creadas
source(file.choose())
#librerias a usar
#importar datos
library(readxl)
#Realizar gráficos
library(ggplot2)
#pronótico
library(forecast)

#Importar serie de ventas(3) y proveedor(1)
ventas <- read.csv(file.choose())
ventas <- as.matrix(ventas[,2])

prov <- read.csv(file.choose())
prov <- as.matrix(prov)

#Verificación de supuestos----
plot(ventas, type="l", col="#435194", main="Serie de ventas", ylab="Millones de pesos", xlab="Mensual")
#Parece tener tendencia
plot(prov[,1], type="l", col="#3f7996", main= "Serie de precios", xlab= "Mensual", ylab="Precio por tonelada")

###Ventas----
#I(1)
ventas_dif <- op_diferencia(ventas,1)
T_ventas <- length(ventas_dif)

#Media movil
ventas_ma <- ventas[12:71] - MA(ventas,12)
T_ventas_ <- length(ventas_ma)

##Autocorrelogramas----

#Original
autocorrelogramas(ventas,T_ventas/4, pacf_(rho_k(ventas,T_ventas/4)), acf_(ventas,T_ventas/4),T_ventas/4)
#AR(16) inc - MA(17) in

#Dif (Seleccionado)
autocorrelogramas(ventas_dif,T_ventas/4, pacf_(rho_k(ventas_dif,T_ventas/4)), acf_(ventas_dif,T_ventas/4),T_ventas/4)
#AR(15) inc

#Media movil
autocorrelogramas(ventas_ma,T_ventas_/4, pacf_(rho_k(ventas_ma,T_ventas_/4)), acf_(ventas_ma,T_ventas_/4), T_ventas_/4)
#MA(1,8,9,10)
#AR(1,8)

#Definición de p y q ----
p_v <- 0 
q_v <- 2

#estimación ----
ventas_est <- arima(ventas, order=c(p_v,1,q_v))
#posiblemente un MA 2, al realizar ensayo y error y tener cumplimiento de supuestos

#residuales ventas
res_v <- ventas_est$residuals

#Cumplimiento de supuestos----
#media cero
mcero_v <- media_cero(res_v,p_v)

sgmaa_v <- sigma_a(res_v,p_v,q_v,mcero_v)

#Estadistico de prueba
ep_v <- est_prueba_med_cero(res_v,p_v,mcero_v,sgmaa_v)
#Es menor a 2

#Verificación ruido blanco
rb <- bp(res_v,p_v,q_v,30)

#Normalidad
normalidad(res_v)

#Aberrancia
aberrancia(res_v)

#Admisibilidad
#Solo miramos invertibilidad al representarlo en un modelo MA
coef_ventas <- ventas_est$coef
invertibilidad_ventas <- Re(polyroot(c(1,coef_ventas[1],coef_ventas[2])))
if (abs(invertibilidad_ventas[1]) > 1 & abs(invertibilidad_ventas[2]) > 1) {
  print("Hay invertibilidad")
} else {
  print("No hay invertibilidad")
}

###Proveedores----
#Estabilización de varianza----
#seleccionar el dato del proveedor
data <- read.csv(file.choose(), sep = ";")
data <- as.matrix(data, ncol=1)

# Tamaño de cada grupo
R <- 12

# Calcula el número total de grupos completos (por eso se usa %/% para quitar el residuo) que se pueden sacar de la muestra 
H <- ceiling(nrow(data) %/% R)

# Inicializa vectores para almacenar las medias y las varianzas
medias <- numeric(H)
desviaciones <- numeric(H)

# Itera a través de los grupos para tener cada media y varianza
for (i in 1:H) {
  # Encuentra los índices para el grupo actual
  inicio <- (i - 1) * R + 1
  fin <- min(i * R, nrow(data))
  
  # Extrae los datos del grupo actual
  grupo <- data[inicio:fin]
  
  # Calcula la media y la varianza del grupo actual
  medias[i] <- promedio(grupo)
  desviaciones[i] <- sqrt(varianza(grupo))
}

# Imprime las medias y las varianzas de cada grupo
for (i in 1:H) {
  cat("Grupo", i, ": Media =", medias[i], ", Desviación Estándar =", desviaciones[i], "\n")
}

#Coeficientes de varianza
lambdas <- c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
L <- length(lambdas)
Mlambda <- matrix(data = NA, nrow=1, ncol = L)
SDlambda <- matrix(data = NA, nrow=1, ncol = L)
#Calcular el denominador 
for (l in 1:L) {
  cvpotencia <- matrix(data = NA, nrow= 1, ncol = H)
  for (h in 1:H){
    cvpotencia[h] <- desviaciones[h]/(medias[h]^(1-lambdas[l]))
  }
  Mlambda[l] <- sum(cvpotencia)/H
  SDlambda[l] <- sqrt(sum((cvpotencia-Mlambda[l])^2)/(H-1))
}
rbind(lambdas, Mlambda, SDlambda)
CVlambda <- matrix(data=NA, ncol = L, nrow = 1)
for (l in 1:L){
  CVlambda[l] <- SDlambda[l]/Mlambda[l]
}
rbind(lambdas, CVlambda)
#Escoger la posición en la que se encuentra el valor que minimiza la varianza que en este caso es 0.2260398 equivalente a lambda -2
index <- which.min(CVlambda)
index

#Transformación serie
prov_estabilizado <- ((data^lambdas[index]))
par(mfrow=c(1,2))
plot(data, type="l", col="#40E0D0")
plot(prov_estabilizado, type="l",col="#4682B4")
#Se puede ver la disminución, comparando el rango del eje y

#Otras transformaciones probadas
#I(1)
prov_dif <- op_diferencia(prov,1)
T_prov_dif <- length(prov_dif)

#media movil
prov_ma <- prov[60:71,] - MA(prov,12) #"No tenia,"
T_prov_ma <- length(prov_ma)

#Estacionalidad
#Promedio estacional
prom_estacional <- promedio_estacional(prov,12)

#Prov sin estacionalidad
prov_se <- prov - rep(prom_estacional,length.out=length(prov))
T_prov_se <- length(prov_se)

#estacionalidad con varianza estabilizada
prom_estacional_2 <- promedio_estacional(prov_estabilizado,12)

#Prov sin estacionalidad y con varianza estabilizada
datosest_se <- prov_estabilizado - rep(prom_estacional_2,length.out=length(prov_estabilizado))
T_datosest_se <- length(datosest_se)

##Autocorrelogramas
#original
autocorrelogramas(prov, (nrow(prov)/4), pacf_(rho_k(prov,nrow(prov)/4)), acf_(prov,nrow(prov)/4),nrow(prov)/4)
#AR(13) - MA(17)

#Sin tendencia estacional
autocorrelogramas(prov_se, (T_prov_se/4), pacf_(rho_k(prov_se,T_prov_se/4)), acf_(prov_se,T_prov_se/4),T_prov_se/4)
#MA(11,13)
#AR(11,13)

#I(1)
autocorrelogramas(prov_dif, T_prov_dif/4, pacf_(rho_k(prov_dif,T_prov_dif/4)), acf_(prov_dif,T_prov_dif/4), T_prov_dif/4)
#MA(1,3,10)
#AR(1,2,9,13,15)

#media movil
autocorrelogramas(prov_ma, T_prov_ma/4, pacf_(rho_k(prov_ma,T_prov_ma/4)), acf_(prov_ma,T_prov_ma/4), T_prov_ma/4)
#MA(12)
#AR(5,11,12,15)

#Con la varianza estabilizada
autocorrelogramas(prov_estabilizado, T_prov_dif/4, pacf_(rho_k(prov_estabilizado,T_prov_dif/4)), acf_(prov_estabilizado,T_prov_dif/4),T_prov_dif/4)
#AR(13) inc
#Se estabiliza la varianza de la serie y se eligé este, el cuál al ser un AR(13) de polinomio incompleto lo representaremos como un MA

#Sin estacionalidad y con varianza estabilizada
autocorrelogramas(datosest_se, T_prov_dif/4, pacf_(rho_k(datosest_se,T_prov_dif/4)), acf_(datosest_se,T_prov_dif/4),T_prov_dif/4)

##Definición de p y q ----
p_p <- 0 # Orden AR (polinomio incompleto)
q_p <- 1 # Orden MA (polinomio incompleto)

#Estimación del modelo para proveedores----
prov_est <- arima(prov_estabilizado, order = c(p_p, 0, q_p))

#residuales proveedores
res_p <- prov_est$residuals

#media cero
mcero_p <- media_cero(res_p,p_p)

sgmaa_p <- sigma_a(res_p,p_p,q_p,mcero_p)

#Estadistico de prueba
ep_p <- est_prueba_med_cero(res_p,p_p,mcero_p,sgmaa_p)
#menor a 2

#Verificación ruido blanco
rb_p <- bp(res_p,p_p,q_p,30)

#Normalidad
normalidad(res_p)

#Aberrancia
aberrancia(res_p)
#IMPORTANTE:
#Tanto para MA(1), MA(2), MA(3) se encontró un solo dato (la posición 22) por fuera del intervalo de las tres desviaciones estandar. Por lo que se decidió realizar un MA(1) por parsimonia.

#Admisibilidad solo miramos invertibilidad porque es un modelo MA******
coef_prov <- prov_est$coef
invertibilidad_prov <- Re(polyroot(c(1,coef_prov[1])))
if (abs(invertibilidad_prov[1]) > 1 ) {
  print("Hay invertibilidad")
} else {
  print("No hay invertibilidad")
}

#Pronóstico fuera de muestra----
#MA2----
par(mfrow=c(1,1))
pron <- ventas_dif
arima(pron, order=c(0,0,2))
theta1 <- as.numeric(arima(pron, order=c(0,0,2))$coef[1])
theta2 <- as.numeric(arima(pron, order = c(0,0,2))$coef[2])
at <- as.numeric(arima(pron, order=c(0,0,2))$res)

p <- 50 #numero de pronosticos a realizar
pp <- length(pron)
pronostico <- numeric(p)
pronostico[1:20] <- pron[52:71]

for (i in 20:p) {
  pronostico[i] <- -theta1*at[i-1]+at[i]-theta2*at[i-2]
}
plot(pronostico, type = "l",las =1, main = "Pronóstico fuera de muestra ", col = "blue")
#Antes de 20 es muestra y despúes es pronóstico

#A partir de la linea azul es muestra
serie_completa <- append(pron,pronostico[21:p])
plot(serie_completa, type="l", main ="Pronóstico fuera de muestra", xlab="Tiempo", ylab ="Pronóstico",col="#BA55D3")+abline(v=pp, lty=2:3, col="cornflowerblue")
plot(forecast(as.numeric(pron)))

#MA1----
par(mfrow=c(1,1))
pron1 <- prov_estabilizado
arima(pron, order=c(0,0,1))
theta1 <- as.numeric(arima(pron, order=c(0,0,1))$coef[1])
at <- as.numeric(arima(pron, order=c(0,0,2))$res)

p <- 50 #numero de pronosticos a realizar
pp <- length(pron)
pronostico <- numeric(p)
pronostico[1:20] <- pron[51:pp]

for (i in 21:p) {
  pronostico[i] <- -theta1*at[i-1]+at[i]
}
plot(pronostico, type = "l",las =1, main = "Pronóstico fuera de muestra ")+abline(v=20, col="lightpink")
#los primeros 20 datos son de la muestra, de ahí para adelante es el pronóstico

serie_completa <- append(pron,pronostico[21:p])
plot(serie_completa, type="l", main ="Pronóstico fuera de muestra completa", xlab="Tiempo", ylab ="Pronóstico",col="#BA55D3")+abline(v=pp, lty=2:3, col="cornflowerblue")
plot(forecast(pron))


#Pronóstico dentro de muestra de ventas de un MA(2)
par(mfrow=c(1,1))
pron <- ventas_dif
arima(pron, order=c(0,0,2))
theta1 <- as.numeric(arima(pron, order=c(0,0,2))$coef[1])
theta2 <- as.numeric(arima(pron, order = c(0,0,2))$coef[2])
at <- as.numeric(arima(pron, order=c(0,0,2))$res)

p <- 50 #numero de pronosticos a realizar
pp <- length(pron)
pronostico <- numeric(p)
pronostico[21:50] <- pron[1:30]

for (i in 2:30) {
  pronostico[i] <- -theta1*at[i-1]+at[i]
}
plot(pronostico, type = "l",las =1, main = "Pronóstico dentro de muestra ", col = "#CDAA7D")+abline(v=30, lty=2:3, col = "#1E90FF")

serie_completa <- append(pron,pronostico[1:30])
plot(serie_completa, type="l", main ="Pronóstico dentro de muestra completa", xlab="Tiempo", ylab ="Pronóstico",col="#BA55D3")+abline(v=pp, lty=2:3, col="cornflowerblue")
plot(forecast(pron))


#PUNTO DE MCO
#OLS para ventas----
#Es un MA(2)
#Defino matriz de diseño
yt_1 <- op_diferencia2(ventas_dif, 1)
yt_1 <- yt_1[1:68]
yt_2 <- op_diferencia2(ventas_dif, 2)
col1 <- rbind(rep(1,68))

matriz_diseno <- t(as.matrix(rbind(col1, yt_1, yt_2)))

#Estimación coeficientes

coe <- as.matrix(solve(t(matriz_diseno)%*%matriz_diseno)%*%t(matriz_diseno)%*%ventas_dif[1:68])

inter <- coe[1]
theta1 <- coe[2]
theta2 <- coe[3]

yestOLS <- matriz_diseno %*% coe
plot(yestOLS, col = "#008f30", ylab = "", xlab = "Tiempo" , main = "Pronóstico por OLS")+ abline(a=inter, b=coe[2], col = "lightgreen")



#Pronóstico ingenuo 
#Para ventas
pro_ingv <- promedio(ventas)
plot(C(pro_ingv))


#Pronóstico dentro de muestra PROVEEDORES MA1 ----
arima(pron1, order=c(0,0,1))
theta1 <- as.numeric(arima(pron1, order=c(0,0,1))$coef[1])
at <- as.numeric(arima(pron1, order=c(0,0,2))$res)

p <- 50 #numero de pronosticos a realizar
pp <- length(pron)
pronostico <- numeric(p)
pronostico[21:50] <- pron1[1:30]

for (i in 2:30) {
  pronostico[i] <- -theta1*at[i-1]+at[i]
}
plot(pronostico, type = "l",las =1, main = "Pronóstico dentro de muestra ", col = "#CDAA7D")+abline(v=30, lty=2:3, col = "#1E90FF")

serie_completa <- append(pron1,pronostico[1:30])
plot(serie_completa, type="l", main ="Pronóstico dentro de muestra", xlab="Tiempo", ylab ="Pronóstico",col="#BA55D3")+abline(v=pp, lty=2:3, col="cornflowerblue")
plot(forecast(pron))



#OLS para proveedores----
#Es un MA(1)
#Defino matriz de diseño
yt_1 <- op_diferencia(prov, 1)

col1 <- cbind(1, length(prov)) 

matriz_diseno <- as.matrix(col1, yt_1)

#Estimación coeficientes

coe <- as.matrix(solve(t(matriz_diseno)%*%matriz_diseno)%*%t(matriz_diseno)%*%ventas)

inter <- coe[1]
theta1 <- coe[2]


yestOLS <- coe%*%matriz_diseno
plot(yestOLS)

#Pronóstico ingenuo 
#Para proveedores
#pronostico ingenuo


#COMPARACIÓN CON MEDIDAS DE AJUSTE 
mse(pron[1:68], yestOLS) 
eam(pron[1:68], yestOLS)


#-