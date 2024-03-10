####----
##Modelo Naive
#Graficar
library(ggplot2)
#Manipulación del df
library(dplyr)
#leer el excel
library(readxl)
excel <- read_xlsx(file.choose(), "Hoja1")

#funciones previamente programadas
source(file.choose())

#flip
flip_excel <- flip(excel[,2])

#Variables del modelo
y <- flip_excel
xdiseno <- sum(y)/nrow(y)
e_est <- y - xdiseno

#Gráfica
df <- data.frame(Tiempo = 1:nrow(y), desempleo = y, naive = xdiseno)
colnames(df) <- c("tiempo", "desempleo", "naive")

grafica_naive <- ggplot(df, aes(x = tiempo)) +
  geom_line(aes(y = desempleo), color = "blue", size = 0.5) +
  geom_line(aes(y = naive), color = "#44414a", linetype = "dashed") +
  labs(title = "Desempleo Colombia 2001/01 - 2023/06", x = "Tiempo (Mensual)", y = "Desempleo (%)") +
  scale_x_continuous(breaks= seq(1, nrow(df), by=24)) +
  theme_minimal()

print(grafica_naive)

####----
##Modelo Naive con dos betas
#variables
y_rz <- as.matrix(y)
xdiseno_rz <- as.matrix(cbind(1,1:nrow(y_rz)))
colnames(xdiseno_rz) <- c("intercepto", "beta")
#estimacion
beta_m2 <- solve(t(xdiseno_rz)%*%xdiseno_rz)%*%t(xdiseno_rz)%*%y_rz

#df para graficar
df_gf2 <- data.frame(1:nrow(y_rz),y_rz,xdiseno_rz)
colnames(df_gf2) <- c("tiempo","desempleo","intercepto","beta")

#grafica m2
grafica_m2 <- ggplot(df_gf2, aes(x = tiempo)) +
  geom_line(aes(y = y_rz), color = "blue", size = 0.5) +
  geom_abline(slope = beta_m2[2,1], intercept = beta_m2[1,1], color = "red") +
  labs(title = "Desempleo Colombia 2001/01 - 2023/06", x = "Tiempo (Mensual)", y = "Desempleo (%)") +
  scale_x_continuous(breaks= seq(1, nrow(df_gf2), by=24)) +
  theme_minimal()

print(grafica_m2)

####----
##Modelo Naive con dos betas y dummies por mes
#variables
y_m3 <- as.matrix(y)
#dummies
dummies<- diag(1,12,12)

T_m3 <- nrow(y_m3)

dummies_rep <- ceiling(T_m3/nrow(dummies))

dummies_matriz <- matrix(0, nrow = T_m3, ncol = ncol(dummies))

for (i in 1:dummies_rep) {
  primera_fila <- (i-1) * nrow(dummies) + 1
  ultima_fila <- min(i * nrow(dummies), T_m3)
  dummies_matriz[primera_fila:ultima_fila,] <- dummies[1:(ultima_fila-primera_fila + 1),]
}

xdiseno_m3 <- as.matrix(cbind(1,1:nrow(y_m3),dummies_matriz))
colnames(xdiseno_m3) <- c("intercepto","t","enero","febrero","marzo","abril","mayo","junio","julio","agosto","septiembre","octubre","noviembre","diciembre")

xdiseno_m3 <- xdiseno_m3[,-3]

#estimación
beta_m3 <- solve(t(xdiseno_m3)%*%xdiseno_m3)%*%t(xdiseno_m3)%*%y_m3
y_est_m3 <- xdiseno_m3 %*% beta_m3

#df para graficar
df_gf3 <- data.frame(1:nrow(y_m3),y_m3,xdiseno_m3)
colnames(df_gf3) <- c("tiempo","desempleo","intercepto","t","febrero","marzo","abril","mayo","junio","julio","agosto","septiembre","octubre","noviembre","diciembre")

#gráfica
grafica_m3 <- ggplot(df_gf3, aes(x = tiempo)) +
  geom_line(aes(y = y_m3), color = "blue", size = 0.5) +
  geom_line(aes(y = y_est_m3), color = "red", size = 0.5) +
  labs(title = "Desempleo Colombia 2001/01 - 2023/06", x = "Tiempo (Mensual)", y = "Desempleo (%)") +
  scale_x_continuous(breaks= seq(1, nrow(df_gf3), by=24)) +
  theme_minimal()

print(grafica_m3)



####----
###Ejercicio estacionalidad

##Importar librerias
#Usado para importar el excel
library(readxl)
##Importar datos y funciones
excel <- read_xlsx(file.choose(), "Hoja1")
source(file.choose())

##definir variables
desempleo <- excel[,2]

##Media movil
media_movil <- MA(datos=excel[,2],ventana=12)

##Serie sin tendencia
serie_sin_tendencia <- desempleo[12:nrow(desempleo),] - media_movil

##Promedios estacionales
num_filas <- ceiling(270/12)
num_datos <- num_filas*12

datos_completos <- c(excel[,2],rep(NA, num_datos - length(excel[,2])))

matriz_promedios <- matrix(data=,ncol=12,byrow=TRUE)

vector_promedios <- numeric(ncol(matriz_promedios))

for (i in 1:ncol(matriz_promedios)) {
  promedio <- mean(matriz_promedios[,i],na.rm = TRUE)
  vector_promedios[i] <- promedio
}

print(vector_promedios)

##Serie desestacionalizada

for (i in 1:nrow(excel[,2])) {
  promedios_completos <- c(excel[,2],rep(NA, ))
}

