#DGP1-----
DGP1 <- function(miu, phi) {
  #Semilla
  set.seed(1028)
  #generación del error
  et <- rnorm(5000,0,1)
  #generación de Xt
  xt <- c()
  for (i in 1:length(et)) {
    xt[i] <- ifelse(i == 1, miu + et[i], miu + phi*xt[i-1] + et[i])
  }
  #quitar la última fila
  et <- et[-1]
  xt <- xt[-1]
  #dejar solo las últimas obs
  xt <- xt[4500:5000]
  #graficar
  g_DGP1 <- plot(xt, type="l", col="blue", ylab = "Xt", xlab= "Obs", main= "DGP1 | AR con media no nula")
}

DGP1(0.5,0.4)

#DGP2----
DGP2 <- function(miu, phi) {
  #semilla
  set.seed(1028)
  #genración del error
  et <- rnorm(5000,0,1)
  #generación de la tendencia
  beta_t <- c()
  for (i in 1:length(et)) {
    beta_t[i] <- i
  }
  #generación de Xt
  xt <- c()
  for (i in 1:length(et)) {
    xt[i] <- ifelse(i == 1, miu + beta_t[i] + et[i], miu + beta_t[i] + phi*xt[i-1] + et[i])
  }
  #quitar la última fila
  et <- et[-1]
  beta_t <- beta_t[-1]
  xt <- xt[-1]
  #dejar obs finales
  xt <- xt[4500:5000]
  #grafica
  g_DGP2 <- plot(xt, type="l", col="blue", ylab = "Xt", xlab= "Obs", main= "DGP2 | AR con tendencia (TS)")
}

DGP2(0.5, 0.4)

#DGP3----
DGP3 <- function(t_obs) {
  #semilla
  set.seed(1028)
  #generacion del error
  et <- rnorm(t_obs, 0, 1)
  #generación de Xt
  xt <- c()
  for (i in 1:t_obs) {
    xt[i] <- ifelse(i == 1, et[i], xt[i-1] + et[i])
  }
  #quitar la última fila
  et <- et[-1]
  xt <- xt[-1]
  #obs finales
  xt <- xt[(t_obs - 500):t_obs]
  #grafica
  g_DGP3 <- plot(xt, type="l", col="blue", ylab = "Xt", xlab= "Tiempo", main= "DGP3 | Caminata aleatoria sin deriva (DS)")
}

yt <- DGP3(1000)

#DGP4----
DGP4 <- function(miu) {
  #semilla
  set.seed(1028)
  #generación del error
  et <- rnorm(5000, 0, 1)
  #generación xt
  xt <- c()
  for (i in 1:length(et)) {
    xt[i] <- ifelse(i == 1, miu + et[i], miu + xt[i-1] + et[i])
  }
  #Quitar la última fila
  et <- et[-1]
  xt <- xt[-1]
  #obs finales
  xt <- xt[4500:5000]
  #grafica
  g_DGP4 <- g_DGP4 <- plot(xt, type="l", col="blue", ylab = "Xt", xlab= "Obs", main= "DGP4 | Caminata aleatoria con deriva (DS)") 
}

DGP4(0.5)

#taos------------------------------------------------------------------------------------
######################################### tau ###############################################
tau <- function(observaciones, semilla, obs_finales) { 
  #semilla
  set.seed(semilla)
  #simulación del error
  et <- rnorm(observaciones,0,1)
  #Simulación Yt
  yt <- c()
  for (i in 1:length(et)) {
    yt[i] <- ifelse(i == 1, et[i], yt[i-1] + et[i])
  }
  #diferencia de yt (y)
  dif_yt <- delta_yt(yt)
  #generacion Yt-1 (xdiseño)
  yt_1 <- as.matrix(op_back(yt,1))
  #eliminación de la primera fila
  dif_yt <- dif_yt[-1]
  yt_1 <- yt_1[-1]
  et <- et[-1]
  #dejar solo valores finales
  dif_yt <- dif_yt[(observaciones-obs_finales):(observaciones-1)]
  yt_1 <- yt_1[(observaciones-obs_finales):(observaciones-1)]
  #regresión
  gamma_est <- solve(t(yt_1)%*%yt_1) %*% t(yt_1) %*% dif_yt
  #varcov
  yest <- yt_1 %*% gamma_est
  eest <- dif_yt - yest
  sigma_e <- (t(eest) %*% eest) / (length(eest) - 1)
  vcv <- sigma_e * solve(t(yt_1)%*%yt_1)
  #tau
  tau <- gamma_est/sqrt(vcv)
}

t <- tau(500,10,100)

tau_sim <- function(simulaciones, observaciones, semilla, obs_finales) {
  sim <- c()
  for (s in 1:simulaciones) {
    sim[s] <- tau(observaciones, semilla+s, obs_finales)
  }
  a <- simulaciones-100
  b <- simulaciones
  plot(density(sim[a:b]), col="blue", ylab="Densidad", xlab="Tau", main="Simulaciones Tau")
  return(sim)
}

x <- tau_sim(500,500,10,100)

vc_tau <- function(sim_tau) {
  vc_tau <- c()
  alpha <- as.matrix(seq(0.01, 0.99, by=0.01))
  for (i in 1:length(alpha)) {
    vc_tau[i] <- quantile(sim_tau, alpha[i])  
    vc_tau <- as.matrix(vc_tau)
  }
  mvc_tau <- as.data.frame(vc_tau, alpha)
  return(mvc_tau)
}

vc_tau(x)

#solo para el 1% y 5%
vc_tau_ <- function(sim_tau) {
  alpha <- c(0.01, 0.05, 0.10)  # Define los cuantiles a calcular
  vc_tau <- numeric(length(alpha))
  
  for (i in 1:length(alpha)) {
    vc_tau[i] <- quantile(sim_tau, alpha[i])
  }
  
  mvc_tau <- data.frame(Quantile = alpha, Value = vc_tau)
  return(mvc_tau)
}

vc_tau_(x)
  

###################################### tao_miu ###########################################
tau_miu <- function(observaciones, semilla, miu, obs_finales) {
  #semilla
  set.seed(semilla)
  #simlulacion del error
  et <- rnorm(observaciones,0,1)
  #Simulacion Yt
  yt <- c()
  for (i in 1:length(et)){
    yt[i] <- ifelse(i == 1, miu + et[i], miu + yt[i-1] + et[i])
  }
  #diferencia yt (y)
  y <- delta_yt(yt)
  #generacion xdiseño
  x <- as.matrix(cbind(1,op_back(yt,1)))
  #eliminando la primera fila
  y <- y[-1]
  x <- as.matrix(x[-1,])
  et <- et[-1]
  #eliminar valores iniciales
  y <- y[(observaciones-obs_finales):(observaciones-1)]
  x <- x[(observaciones-obs_finales):(observaciones-1),]
  #regresion
  coef <- solve(t(x)%*%x) %*% t(x) %*% y
  alpha <- coef[1,1]
  gamma <- coef[2,1]
  #varcov
  y_est <- x %*% coef
  e_est <- y - y_est
  sigma_e <- as.numeric((t(e_est) %*% e_est) / (length(e_est)-1))
  vcv <- sigma_e * solve(t(x)%*%x)
  tau_miu <- gamma / sqrt(vcv[2,2])
}

taumiu <- tau_miu(200,10,5,100)

taumiu_sim <- function(simulaciones,observaciones,semilla,miu,obs_finales) {
  sim <- c()
  for (s in 1:simulaciones) {
    sim[s] <- tau_miu(observaciones,semilla+s,miu,obs_finales)
  }
  plot(density(sim[(observaciones-obs_finales):(observaciones-1)]), col="blue", ylab="Densidad", xlab="Tau_miu", main="Simulaciones Tau miu")
  return(sim)
}

tm <- taumiu_sim(5000,500,10,5,100)

vc_taumiu <- function(sim_tau) {
  alpha <- c(0.01, 0.05, 0.10)  # Define los cuantiles a calcular
  vc_tau <- numeric(length(alpha))
  
  for (i in 1:length(alpha)) {
    vc_tau[i] <- quantile(sim_tau, alpha[i])
  }
  
  mvc_tau <- data.frame(Quantile = alpha, Value = vc_tau)
  return(mvc_tau)
}

vc_taumiu(tm)

###################################### tao_tao ###########################################  
tao_tao <- function(semilla, observaciones, obs_finales, miu) {
  #semilla
  set.seed(semilla)
  #simlulacion del error
  et <- rnorm(observaciones,0,1)
  #beta (contador)
  beta <- 1:length(et)
  #Simulacion Yt
  yt <- c()
  for (i in 1:length(et)){
    yt[i] <- ifelse(i == 1, miu + beta[i] + et[i], miu + beta[i] + yt[i-1] + et[i])
  }
  #diferencia yt (y)
  y <- delta_yt(yt)
  #generación xdiseño
  x <- as.matrix(cbind(1,beta, op_back(yt,1)))
  #eliminando la primera fila
  y <- y[-1]
  x <- as.matrix(x[-1,])
  et <- et[-1]
  #eliminar valores iniciales
  y <- y[(observaciones-obs_finales):(observaciones-1)]
  x <- x[(observaciones-obs_finales):(observaciones-1),]
  #regresion
  coef <- solve(t(x)%*%x) %*% t(x) %*% y
  alpha <- coef[1,1]
  beta_est <- coef[2,1]
  gamma <- coef[3,1]
  #varvov
  y_est <- x %*% coef
  e_est <- y - y_est
  sigma_e <- as.numeric((t(e_est) %*% e_est) / (length(e_est)-3))
  vcv <- sigma_e * solve(t(x)%*%x) 
  tao <- gamma / sqrt(vcv[3,3])
}

tt <- tao_tao(10, 500, 100, 10)

taotao_sim <- function(simulaciones,semilla, observaciones, obs_finales, miu) {
  sim <- c()
  for (s in 1:simulaciones){
    sim[s] <- tao_tao(semilla+s, observaciones, obs_finales, miu)
  }
  plot(density(sim[(observaciones-obs_finales):(observaciones-1)]), col="blue", ylab="Densidad", xlab="Tau_tao", main="Simulaciones Tau tao")
  return(sim)
}

ttt <- taotao_sim(5000, 10, 500, 100, 10)

vc_taotao <- function(simulacion) {
  alpha <- c(0.01, 0.05, 0.10)  # Define los cuantiles a calcular
  vc_tau <- numeric(length(alpha))
  
  for (i in 1:length(alpha)) {
    vc_tau[i] <- quantile(simulacion, alpha[i])
  }
  
  mvc_tau <- data.frame(Quantile = alpha, Value = vc_tau)
  return(mvc_tau)
}

tttt <- vc_taotao(ttt)

##3 serie cointegradas
#yt
et <- rnorm(100, 0, 1)
yt <- c()
for (i in 1:100) {
  yt[i] <- ifelse(i == 1, rnorm(1) + et[i], yt[i-1] + et[i])
}
#xt
xt <- 5 + yt + rnorm(100, 0, 1.2)
#zt
zt <- 15 + yt + rnorm(100, 0, 1.5)
plot(yt, type = "l", col="blue", main="3 series cointegradas", xlim=c(1,100), ylim=range(c(yt, xt, zt)))
lines(xt, type = "l", col="#4c5b5e", main="Xt")
lines(zt, type = "l", col="#213e75", main="Zt")

#2 cointegradas y 1 no
#yt
et <- rnorm(100, 0, 1)
yt <- c()
for (i in 1:100) {
  yt[i] <- ifelse(i == 1, rnorm(1) + et[i], yt[i-1] + et[i])
}
#xt
xt <- 5 + yt + rnorm(100, 0, 2)
#zt
zt <- c()
et2 <- rnorm(100,0,2)
for (i in 1:100) {
  zt[i] <- ifelse(i == 1, rnorm(1) + et2[i], zt[i-1] + et2[i])
}
plot(yt, type = "l", col="blue", main="2 series cointegradas y 1 no", xlim=c(1,100), ylim=range(c(yt, xt, zt)))
lines(xt, type = "l", col="black", main="Xt")
lines(zt, type = "l", col="purple", main="Zt")