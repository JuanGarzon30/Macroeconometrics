####----ACF
acf_ <- function(Zt, max_k) {
  T <- length(Zt)
  acf_vector <- c()
  acf_plots <- list()
  
  if (T > max_k) {
    for (k in 1:max_k) {
      numerador <- 0
      for (t in 1:(T - k)) {
        numerador <- numerador + (Zt[t] - mean(Zt)) * (Zt[t + k] - mean(Zt))
      }
      denominador <- sum((Zt - mean(Zt))^2)
      resultado_acf <- numerador / denominador
      acf_vector <- append(acf_vector, resultado_acf)
      
      # grafica
      df_acf <- data.frame(x = 1:k, y = acf_vector)
      g_acf <- ggplot(data = df_acf, aes(x = x, y = y)) +
        geom_segment(aes(xend = x, yend = 0), color = "blue", size=1.5) +
        labs(title = "Correlograma función de autocorrelación",
             subtitle = paste("Con k =", k),
             x = "Distancia (k)",
             y = "Correlación (ρ)") +
        ylim(-1, 1) +
        theme_minimal() +
        scale_x_continuous(breaks = 1:k, minor_breaks = NULL)
      
      acf_plots[[k]] <- g_acf
    }
    print(acf_plots)
    return(acf_vector)
  } else {
    return("El valor asignado a k no puede ser mayor a T")
  }
}

####---->Bandas ACF
bandas <- function(datos, q) {
  T <- length(datos)
  
  for (j in 1:q) {
    suma <- rho_k(datos, j)^2
  }
  
  neg <- -2 * (1/T) * (1 + 2 * suma)
  pos <- 2 * (1/T) * (1 + 2 * suma)
  
  var_estimada <- rbind(neg,pos)
  
  return(var_estimada)
}

####----Phikk
phi_kk <- function(datos,k) {
  #vector phis
  phi_vec <- c()
  
  #matriz inicial
  for (i in 1:k) {
    matriz <- matrix(NA,nrow=k, ncol=k)
    
    #diagonal
    for (i in 1:k) {
      matriz[i, i] <- 1
    }
    
    #llenar la matriz
    rhos <- rho_k(datos,k)
    
    for (z in 1:i) {
      for (j in 1:i) {
        if (z != j) {
          matriz[z, j] <- rhos[abs(z - j)]
          matriz[j, z] <- rhos[abs(z - j)]
        }
      }
    }
    matriz_num <- matriz
    matriz_num[,i] <- rho_k(datos,i)
    # matriz denominador
    matriz_den <- matriz
    #phi kk
    phi_kk <- det(matriz_num)/det(matriz_den)
    phi_vec <- append(phi_vec,phi_kk)
  }
  return(phi_vec[i])  
}

####----PACF
pacf_ <- function(datos, max_lag) {
  rhos_ <- numeric(max_lag) 
  #Hallar los phi_kk hasta el k deseado
  for (k in 1:max_lag) {
    rhos_[k] <- phi_kk(datos, k)
  }
  print(rhos_)
  
  # Grafica
  df_pacf <- data.frame(x = 1:max_lag, y = rhos_)
  g_pacf <- ggplot(data = df_pacf, aes(x = x, y = y)) +
    geom_segment(aes(xend = x, yend = 0), color = "blue", size = 1.5) +
    geom_hline(yintercept = c(2/sqrt(length(datos)), -2/sqrt(length(datos))), linetype = "dashed") +
    labs(title = "Correlograma función de autocorrelación parcial",
         subtitle = paste("Para k =", max_lag),
         x = "Distancia(k)",
         y = "Correlación (Φ)") +
    ylim(-1, 1) +
    theme_minimal() +
    scale_x_continuous(breaks = 1:max_lag, minor_breaks = NULL)
  print(g_pacf)
}

####Flip----
flip <- function(datos) {
  #creacion del nuevo vector invertido
  T <- nrow(datos)
  nuevo_vector <- data.frame(matrix(NA, nrow=T, ncol=ncol(datos)))
  #for que itere a lo largo de las filas del vector inicial  
  for (i in (1:T)) {
    nuevo_vector[i,] <- datos[T-i+1,]
  }
  return(nuevo_vector)
}

####gamma_k----
gamma_k <- function(datos, k) {
  T <- length(datos)
  c <- numeric()
  #covarianza
  for (i in 1:(T-k)) {
    cov <- ((datos[i]-mean(datos))*(datos[i+k]-mean(datos)))
    c <- append(c,cov)
  }
  datos <- sum(c)/T
  return(datos)
}

####rho_k----
rho_k <- function(datos, k) {
  T <- length(datos)
  c_ <- c()
  for (i in 1:k) {
    gma_k <-  gamma_k(datos, i)
    c_ <- append(c_,gma_k)
  }
  rho <- c_/gamma_0(datos)
  return(rho)
}

####gamma_0----
gamma_0 <- function(datos) {
  T <- length(datos)
  media <- sum(datos)/T
  var_ <- sum((datos-media)^2)/T
  return(var_)
}

####----Operador Back
op_back <- function(datos, exponente) {
  #creación del nuevo vector vacio con el nuevo de filas de "datos"
  T <- length(datos)
  nuevo_vector <- numeric(T)
  #iteración para que llene el nuevo vector con el desplazamiento del dato o NA según   corresponda
  for (i in 1:T) {
    nuevo_vector[i] <- ifelse(i - exponente >= 1, datos[i - exponente], NA)
  }
  
  return(nuevo_vector)
}

####----Operador diferencia
op_diferencia <- function(datos, spam) {
#Creación del vector para el máximo de diferencias (T-k)
  T <- nrow(datos)
  diferencias <- numeric(T - spam)
#Convertir los datos en un dataframe
  datos_t <- as.data.frame(datos)
#Dejar un vector vacio para rellenarlo con la resta mediante el for
  nuevo_vector <- c()
#Es necesario que el spam sea menor a la cantidad de observaciones
  if (T > spam) {
    for (i in (spam + 1):T) { 
      resta <- datos_t[i - spam, ] - datos_t[i, ]
      nuevo_vector <- append(nuevo_vector, resta)
    }
    return(nuevo_vector) 
  } else {
    return("El spam de la serie no puede ser mayor al número de observaciones")
  }
}

####----Operador forward
op_forward <- function(datos,exponente) {
#Se crea un vector vacio que tenga el numero de filas de los datos utilizados
  T <- nrow(datos)
  vector_nuevo <- numeric(T)
#Se utiliza un for para que itere por todas las filas desplazando los valores de acuerdo con el exponente ingresado, llenando de NA para valores posteriores a T que fueron desplazados
  for (i in 1:T) {
    if (i + exponente <= T) {
      vector_nuevo[i] <- datos[i+exponente,]  
    }
    else {
      vector_nuevo[i] <- NA
    }
  }
  return(vector_nuevo)
}

####----Variación porcentual
variacion_porcentual <- function(datos, frecuencia,columna) {
#creación del nuevo vector con el numero de filas en las que se puede encontrar la variación
  T <- nrow(datos)
  nuevo_vector <- numeric(T - frecuencia)
#Solo se podrá realizar si el número de observaciones es mayor a la frecuencia
  if (T > frecuencia) {
#mediante el for se está calculando la variación como: (valor final / valor inicial) - 1
    for (i in 1:(T-frecuencia)) {
      operacion <- (datos[i,columna]/ datos[i + frecuencia,columna]) - 1
      nuevo_vector[i + frecuencia] <- operacion
    }
    return(nuevo_vector)
  } else {
    return("La frecuencia en la que se está comparando no debe ser mayor al número de observaciones de la serie de datos")
  }
}

####---- MA
MA <- function(datos,ventana) {
#crear una nueva matriz con el numero de filas de los datos seleccionados
  T <- nrow(datos)
  matriz <- matrix(NA,nrow=(T-ventana+1), ncol=1)
#realizar un for desde 1 hasta T-k+1 siendo k la ventana de tiempo sobre la que se está realizando la media movil
  for (i in 1:(T-ventana+1)) {
    media <- sum(datos[i:(i+ventana-1),1]) / ventana
    matriz[i,1] <- media
  }
  return(matriz)
}
  #c:Es necesario poner datos=excel[,columna a usar]
####----PACF

var_p <- function(datos, q) {
  T <- length(datos)
  
  positivo <- numeric(q)
  negativo <- numeric(q)
  
  pos <- c()
  
  for (i in 1:q) {
    rhos <- rho_k(datos, i)
    positivo <- 2 *(1 / T) * (1 + 2 * sum(rhos^2))
    pos <- append(positivo,pos)
  }
  
  neg <- c()
  
  for (i in 1:q) {
    rhos <- rho_k(datos, i)
    negativo <- -2 *(1 / T) * (1 + 2 * sum(rhos[1:i]^2))
    neg <- append(negativo,neg)
  }
  
  #df para graficar
  df_bandas <- data.frame(pos, neg)
  
  return(df_bandas)  
}

####Autocorrelogramas----
autocorrelogramas <- function(datos,k, pacf_, acf_, q) {
  
  #df almacenando pacf
  pacf <- pacf_(rho_k(datos,k))
  df_pacf <- data.frame(x1 = 1:k, y1 = pacf)
  
  #Gráfica pacf
  g_pacf <- ggplot(data = df_pacf, aes(x1,y1.y)) +
    geom_bar(stat = "identity", fill = "#5999e2", color = "#5999e2") +
    labs(x = "Distancia (k)", y = "Φkk", title = "PACF") +
    theme_minimal() +
    geom_hline(yintercept = c(-2/sqrt(length(datos)), 2/sqrt(length(datos))), linetype="dashed", color= "#999999") +
    scale_x_continuous(breaks = 1:k, minor_breaks = NULL)
  
  print(g_pacf)
  
  #df almacenando acf
  acf <- acf_(datos,k)
  df_acf <- data.frame(x= 1:k, y= acf)
  
  #intervalo de confianza
  df_bandas <- var_p(datos,q)
  
  #IC
  ic_p <- qnorm(0.975)/sqrt(length(datos))
  ic_n <- -(qnorm(0.975)/sqrt(length(datos)))
  
  #Gráfica acf
  g_acf <- ggplot(data = df_acf, aes(x,y)) +
    geom_bar(stat = "identity", fill = "#5999e2", color = "#5999e2") +
    labs(x = "Distancia (k)", y = "ρk", title = "ACF") +
    geom_hline(yintercept = c(ic_p, ic_n), linetype="dashed", color= "#999999") +
    theme_minimal() +
    scale_x_continuous(breaks = 1:k, minor_breaks = NULL)
  
  print(g_acf)
}

####Ljung-Box----
ljung_box <- function(datos,p,k) {
  T <- length(datos)
  suma <- (sum(rho_k(datos,k))^2)/(T-p-k)
  resultado <- (T-p)*(T-p+2)*suma
  return(resultado)
}

#####Box-Pierce
bp <- function(res, p, q, K){
  #largo del residual
  T <- length(res)
  #renombrar el residual
  res_1 <- res
  #var
  var <- sum(res^2)
  #sumatoria
  resultado <- c()
  for(k in 1:K){
    #residual rezagado
    res_2 <- op_back(res_1,k)
    #multiplicación
    m <- na.omit(res_1*res_2)
    #covarianza
    cov <- sum(m)
    #Resultado
    resultado<- append(resultado,cov)
  }
  #Correlación
  corr <- resultado/var
  #Calculo Box-Pierce
  bp_ <- (T-p) * sum(corr^2)
  #Calculo p-value
  p_value <- 1-pchisq(bp_,(K-p-q))
  #Es necesario el cumplimiento del mensaje
  Mensaje <- "Recuerde que K>20"
  
  if(0.05 < p_value){
    resultados <- data.frame(bp_,p_value,"Ruido blanco",Mensaje)
    return(resultados)
  }
  else{
    resultados <- data.frame(bp_,p_value,"No hay ruido blanco",Mensaje)
    return(resultados)
  }
}

#####Normalidad
normalidad <- function(res) {
  T <- length(res)
  #desviación estandar del error
  sd <- sqrt(varianza(res))
  #Intervalos de confianza
  ic_pos <- 2*sd
  ic_neg <- -2*sd
  #¿Normalidad?
  res_dentro <- sum(res >= ic_neg & res <= ic_pos)
  porcentaje_dentro <- (res_dentro/T) * 100
  if (porcentaje_dentro >= 95) {
    print("Hay normalidad")
  } else {
    print("No hay normalidad")
  }
}

#####Media cero----
media_cero <- function(res,p) {
  T <- length(res)
  m_c <- sum(res[(p+1):T]) / (T - p)
  return(m_c)
}

####Varianza residuales----
sigma_a <- function(res,p,q,med_cero) {
  T <- length(res)
  sg_a <- sqrt(sum((res[(p + 1):T] - med_cero)^2) / (T - p - q))
  return(sg_a)
}

#estadistico media cero----
est_prueba_med_cero <- function(res,p,med_cero,sgma_a) {
  T <- length(res)
  e_p <- abs((sqrt(T-p)*med_cero)/(sgma_a))
  print(e_p)
  if (e_p>=2) {
    print("No se rechaza H1")
  }
  else {
    print("Se rechaza H1")
  }
}

#rho_error
rho_k_error <- function(res,k,p,med_cero) {
  T <- length(res)
  for (t in (p+1):(T-k)) {
    numerador <- res[t]*res[t+k]
  }
  for (t in (p+1):T) {
    denominador <- res[t]^2
  }
  rhok_error <- numerador/denominador
  return(rhok_error)
}

#Aberrancia----
aberrancia <- function(res) {
  T <- length(res)
  #desviación estandar del error
  sd <- sqrt(varianza(res))
  #Intervalos de confianza
  ic_pos <- 3*sd
  ic_neg <- -3*sd
  #¿Aberrancia?
  if(any(res>ic_pos | res<ic_neg)) {
    print("Hay algun residual fuera de las 3 desviaciones estandar")
  } else {
    print("No hay algun residual fuera de las 3 desviaciones estandar")
  }
}

#Estabilidad----
estabilidad <- function(data_a, data_b) {
  num <- covarianza(data_a, data_b, n=length(data_a))
  den <- desv_est(data_a, n = length(data_a))*desv_est(data_b, n = length(data_b))
  num/den
}

pacf__ <- function(datos,k) {
  #vector phis
  phi_vec <- c()
  
  #matriz inicial
  for (i in 1:k) {
    matriz <- matrix(NA,nrow=k, ncol=k)
    
    #diagonal
    for (i in 1:k) {
      matriz[i, i] <- 1
    }
    
    #llenar la matriz
    rhos <- rho_k(datos,k)
    
    for (z in 1:i) {
      for (j in 1:i) {
        if (z != j) {
          matriz[z, j] <- rhos[abs(z - j)]
          matriz[j, z] <- rhos[abs(z - j)]
        }
      }
    }
    matriz_num <- matriz
    matriz_num[,i] <- rho_k(datos,i)
    # matriz denominador
    matriz_den <- matriz
    #phi kk
    phi_kk <- det(matriz_num)/det(matriz_den)
    phi_vec <- append(phi_vec,phi_kk)
  }
  return(phi_vec[i])  
}

#Promedio Estacional
promedio_estacional <- function(datos,frecuencia){
  
  max <- (length(datos)/frecuencia)
  
  #Matriz vacia
  m1 <- matrix(NA, nrow = frecuencia, ncol = ceiling(max))
  
  #matriz llena
  m1[1:length(datos)] <- datos
  
  resultado <- c()
  
  for(i in 1:nrow(m1)){
    res <- promedio(na.omit(m1[i,]))
    resultado <- append(resultado,res)
  }
  #matriz
  resultado <- as.matrix(resultado) 
  d <- as.matrix(resultado[1,])
  resultado <- resultado[-1,]
  resultado <- as.matrix(resultado)
  resultado <- matrix(rbind(resultado,d),ncol=1)
  
  promedio <- as.matrix(resultado)
  
  return(promedio)
}

delta_yt <- function(datos) {
  yt_lag <- op_back(datos,1)
  resta <- as.matrix(datos-yt_lag)
  return(resta)
}

#Operador diferencias (DF)
op_diferencias_ <- function(datos, spam) {
  # Creación del vector para el máximo de diferencias (T-k)
  T <- length(datos)
  
  # Es necesario que el spam sea menor a la cantidad de observaciones
  if (T > spam) {
    nuevo_vector <- numeric(T - spam)
    
    for (i in (spam + 1):T) {
      resta <- datos[i] - datos[i - spam]
      nuevo_vector[i - spam] <- resta
    }
    
    return(nuevo_vector)
  } else {
    return("El spam de la serie no puede ser mayor al número de observaciones")
  }
}

#Operador diferencia con NA p veces
op_diferencias <- function(datos, spam) {
  T <- length(datos)
  nuevo_vector <- numeric(T)
  
  if (T > spam) {
    for (i in (spam + 1):T) {
      resta <- datos[i] - datos[i - spam]
      nuevo_vector[i] <- resta
    }
  } else {
    return("El spam de la serie no puede ser mayor al número de observaciones")
  }
  
  nuevo_vector[1:spam] <- NA  # Establecer las primeras posiciones como NA
  
  return(nuevo_vector)
}

#regresion
regresion <- function(x,y) {
  coef <- solve(t(x) %*% x) %*% t(x) %*% y
  return(coef)
}

#Ruido blanco
rb <- function(e_est) {
  rb_p1_E3 <- bp(e_est,1,0,30)
  if (rb_p1_E3[1,3] == "Ruido blanco") {
    rta1_p1_E3 <- TRUE
  } else {
    rta2_p1_E3 <- FALSE
  }
}

#parsimonia
parsimonia <- function(e_est, x) {
  #parametros iniciales
  T <- nrow(x)
  k <- ncol(x)
  #sigma_e
  sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
  vcv <- sgm_e * solve(t(na.omit(x))%*%na.omit(x))
  #IC
  ic <- c(gamma - 2*sqrt(vcv[3,3]), gamma + 2*sqrt(vcv[3,3]))
  #significancia
  if (0 >= ic[1] & 0 <= ic[2]) {
    nosignificante_p1_E3 <- TRUE
  } else {
    nosignificante_p1_E3 <- FALSE
  }
  return(nosignificante_p1_E3)
}

parsimonia_MCE <- function(x, y, p_opt) {
  #parametros iniciales
  t <- nrow(x)
  k <- ncol(x)
  #estimacion
  coef <- regresion(x,y)
  e_est <- y - (x%*%coef)
  #vcv
  sgm_e <- as.numeric((t(e_est)%*%(e_est))/(t-k))
  vcv <- sgm_e * solve(t(x)%*%(x))
  gamma <- coef[1,]
  #IC
  ic <- c(gamma - 2*sqrt(vcv[1,1]), gamma + 2*sqrt(vcv[1,1]))
  #significancia
  if (0 >= ic[1] & 0 <= ic[2]) {
    significante_p1_E3 <- FALSE
  } else {
    significante_p1_E3 <- TRUE
  }
  return(significante_p1_E3)
}

#Columnas ADF
col_aum <- function(yt,p) {
  y_ <- op_back(yt,1)
  m <- matrix(NA, length(y_), (p_opt-1), TRUE)
  for (i in 2:p) {
    m[,(i-1)] <- as.matrix(op_diferencias(y_,(i-1)))
  }
  return(m)
}

################
#Básicas----
promedio <- function(datos) {
  T <- length(datos)
  media <- sum(datos)/T
  return(media)
}

varianza<-function(datos){
  T<-length(datos)
  varianza<-(sum((datos-promedio(datos))^2))/(T-1)
  return(varianza)
}

sd <- function(datos) {
  varianza <- varianza(datos)
  resultado <- sqrt(varianza)
  return(resultado)
}

covarianza <- function(data_a, data_b, n = length(data_a)) {
  sum((data_a - promedio(data_a))*(data_b - promedio(data_b))/(n-1))
}


