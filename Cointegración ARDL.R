#El presente modelo propuesto busca realizar un análisis de las tasas de corto y largo plazo bajo la cointegración por ARDL
#Preambulo
#FUNCIONES#Básicas-------------------
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

#Operadores
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

delta_yt <- function(datos) {
   yt_lag <- op_back(datos,1)
   resta <- as.matrix(datos-yt_lag)
   return(resta)
}

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

regresion_na <- function(x,y) {
   coef <- solve(t(na.omit(x)) %*% na.omit(x)) %*% t(na.omit(x)) %*% na.omit(y)
   return(coef)
}

#Box-Pierce / RB
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
   p_value <- qchisq(0.95, K - p)
   #Es necesario el cumplimiento del mensaje
   Mensaje <- "Recuerde que K>20"
   
   if(bp_ < p_value){
      resultados <- data.frame(bp_,p_value,"Ruido blanco",Mensaje)
      return(resultados)
   }
   else{
      resultados <- data.frame(bp_,p_value,"No hay ruido blanco",Mensaje)
      return(resultados)
   }
}

rb <- function(e_est, p_opt) {
   rb_p1_E3 <- bp(e_est,p_opt,0,120)
   if (rb_p1_E3[1,3] == "Ruido blanco") {
      rta1_p1_E3 <- TRUE
   } else {
      rta2_p1_E3 <- FALSE
   }
}

#Parsimonia
parsimonia <- function(e_est, x) {
   #parametros iniciales
   T <- nrow(x)
   k <- ncol(x)
   #sigma_e
   sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
   vcv <- sgm_e * solve(t(x)%*%x)
   #IC
   ic <- c(gamma - 2*sqrt(vcv[k,k]), gamma + 2*sqrt(vcv[k,k]))
   #significancia
   if (0 >= ic[1] & 0 <= ic[2]) {
      significante_p1_E3 <- FALSE
   } else {
      significante_p1_E3 <- TRUE
   }
   return(significante_p1_E3)
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

## Dickey Fuller Aumentado
#Columnas ADF
col_aum <- function(yt,p) {
   y_ <- op_back(yt,1)
   m <- matrix(NA, length(y_), (p_opt-1), TRUE)
   for (i in 2:p) {
      m[,(i-1)] <- as.matrix(op_diferencias(y_,(i-1)))
   }
   return(m)
}

col_aum_MCE <- function(yt,p_opt) {
   y_ <- op_back(yt,1)
   m <- matrix(NA, length(y_), (p_opt-1), TRUE)
   for (i in 2:p_opt) {
      m[,(i-1)] <- as.matrix(op_diferencias(y_,(i-1)))
   }
   return(m)
}

# Selección p óptimo
#E3
P_opt_E3 <- function(yt) {
   #Pasar datos a Yt y revisar length/nrow
   yt <- as.matrix(yt)
   #Y diferenciado
   y <- as.matrix(delta_yt(yt))
   
   #xdiseño
   contador <- 1:length(yt)
   col_fijas <- as.matrix(cbind(1,contador, op_back(yt,1)))
   
   #eliminar NAs iniciales
   y <- as.matrix(y[-1,])
   x <- as.matrix(col_fijas[-1,])
   
   #regresion
   coef <- regresion(x,y)
   gamma <- coef[3,1]
   
   #Error estimado
   e_est <- y - (x %*% coef)
   
   p_opt <- 1
   
   #Ruido blanco
   RB <- rb(e_est, p_opt)
   
   #Parsimonia
   Parsimonia <- parsimonia(e_est, x)
   
   #¿DF o ADF?
   if ( RB & Parsimonia == TRUE) {
      p_opt <- 1
      rta_p1 <- TRUE
   } else {
      p_opt <- 2
      rta_p1 <- FALSE
   }
   
   #P óptimo en ADF
   if (rta_p1 == FALSE) {
      while(TRUE) {
         #columnas ADF
         m <- col_aum(yt, p_opt)
         #matriz completa
         m_completa <- as.matrix(cbind(col_fijas,m))
         #quitar NAs iniciales
         x <- as.matrix(m_completa[-(1:p_opt),])
         y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
         #regresión
         coef <- regresion(x,y)
         gamma <- coef[nrow(coef), 1]
         #error
         e_est <- y - x %*% coef
         #ruido blanco
         RB <- rb(e_est, p_opt)
         #parsimonia
         Parsimonia <- parsimonia(e_est, x)
         #P óptimo
         if (RB & Parsimonia == TRUE) {
            p_optimo <- p_opt
            rta <- print(paste("El P óptimo bajo el esquema 3 es:", p_optimo))
            break
         } else {
            p_opt <- p_opt + 1
         }
      }
   } else {
      p_optimo <- p_opt
      rta <- print(paste("El P óptimo bajo el esquema 3 es:", p_optimo))
   }
   return(rta)
}

#E2
P_opt_E2 <- function(yt) {
   #Pasar datos a Yt y revisar length/nrow
   yt <- as.matrix(yt)
   #Y diferenciado
   y <- as.matrix(delta_yt(yt))
   
   #xdiseño
   col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
   
   #eliminar NAs iniciales
   y <- as.matrix(y[-1,])
   x <- as.matrix(col_fijas[-1,])
   
   #regresion
   coef <- regresion(x,y)
   gamma <- coef[2,1]
   
   #Error estimado
   e_est <- y - (x %*% coef)
   
   p_opt <- 1
   
   #Ruido blanco
   RB <- rb(e_est, p_opt)
   
   #Parsimonia
   Parsimonia <- parsimonia(e_est, x)
   
   #¿DF o ADF?
   if ( RB & Parsimonia == TRUE) {
      p_opt <- 1
      rta_p1 <- TRUE
   } else {
      p_opt <- 2
      rta_p1 <- FALSE
   }
   
   #P óptimo en ADF
   if (rta_p1 == FALSE) {
      while(TRUE) {
         #columnas ADF
         m <- col_aum(yt, p_opt)
         #matriz completa
         m_completa <- as.matrix(cbind(col_fijas,m))
         #quitar NAs iniciales
         x <- as.matrix(m_completa[-(1:p_opt),])
         y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
         #regresión
         coef <- regresion(x,y)
         gamma <- coef[nrow(coef), 1]
         #error
         e_est <- y - x %*% coef
         #ruido blanco
         RB <- rb(e_est, p_opt)
         #parsimonia
         Parsimonia <- parsimonia(e_est, x)
         #P óptimo
         if (RB & Parsimonia == TRUE) {
            p_optimo <- p_opt
            rta <- print(paste("El P óptimo bajo el esquema 2 es:", p_optimo))
            break
         } else {
            p_opt <- p_opt + 1
         }
      }
   } else {
      p_optimo <- p_opt
      rta <- print(paste("El P óptimo bajo el esquema 2 es:", p_optimo))
   }
   return(rta)
}

#E1
P_opt_E1 <- function(yt) {
   #Pasar datos a Yt y revisar length/nrow
   yt <- as.matrix(yt)
   #Y diferenciado
   y <- as.matrix(delta_yt(yt))
   
   #xdiseño
   col_fijas <- as.matrix(cbind(op_back(yt,1)))
   
   #eliminar NAs iniciales
   y <- as.matrix(y[-1,])
   x <- as.matrix(col_fijas[-1,])
   
   #regresion
   coef <- regresion(x,y)
   gamma <- coef[1,1]
   
   #Error estimado
   e_est <- y - (x %*% coef)
   
   p_opt <- 1
   
   #Ruido blanco
   RB <- rb(e_est, p_opt)
   
   #Parsimonia
   Parsimonia <- parsimonia(e_est, x)
   
   #¿DF o ADF?
   if ( RB & Parsimonia == TRUE) {
      p_opt <- 1
      rta_p1 <- TRUE
   } else {
      p_opt <- 2
      rta_p1 <- FALSE
   }
   
   #P óptimo en ADF
   if (rta_p1 == FALSE) {
      while(TRUE) {
         #columnas ADF
         m <- col_aum(yt, p_opt)
         #matriz completa
         m_completa <- as.matrix(cbind(col_fijas,m))
         #quitar NAs iniciales
         x <- as.matrix(m_completa[-(1:p_opt),])
         y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
         #regresión
         coef <- regresion(x,y)
         gamma <- coef[nrow(coef), 1]
         #error
         e_est <- y - x %*% coef
         #ruido blanco
         RB <- rb(e_est, p_opt)
         #parsimonia
         Parsimonia <- parsimonia(e_est, x)
         #P óptimo
         if (RB & Parsimonia == TRUE) {
            p_optimo <- p_opt
            rta <- print(paste("El P óptimo bajo el esquema 1 es:", p_optimo))
            break
         } else {
            p_opt <- p_opt + 1
         }
      }
   } else {
      p_optimo <- p_opt
      rta <- print(paste("El P óptimo bajo el esquema 1 es:", p_optimo))
   }
   return(rta)
}


# ADF
#Paso 1 - E3
Paso1_E3 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1,contador, op_back(yt,1)))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[3,1]
      
      #vcv
      #parametros iniciales
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      #taotao
      tt <- gamma / sqrt(vcv[k,k])
      tt_5 <- -3.41
      if (tt < tt_5) {
         print("Es un DGP1 | No hay RU")
      } else {
         print("Es un DGP3 | Hay RU")
      }
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1,contador, op_back(yt,1)))
      m <- col_aum(yt, p_optimo)
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #eliminar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_optimo),])
      y <- as.matrix(delta_yt(yt)[-(1:p_optimo),])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[3,1]
      
      #vcv
      #parametros iniciales
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      #taotao
      tt <- gamma / sqrt(vcv[k,k])
      tt_5 <- -3.41
      if (tt < tt_5) {
         print("Es un DGP1 | No hay RU")
      } else {
         print("Es un DGP3 | Hay RU")
      }
   }
}

#Paso 2 - E3
Paso2_E3 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1,contador, op_back(yt,1)))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #parametros
      j <- 2
      T <- nrow(x)
      k <- ncol(x)
      
      #regresion
      coef_nr <- regresion(x,y)
      
      #No restringido
      SRCnr <- sum((y - (x %*% coef_nr))^2)
      
      #Restringido
      x_restringido <- x[,1]
      coef_r <- regresion(x_restringido, y)
      SRCr <- sum((y - (x_restringido %*% coef_r))^2)
      
      #Phi3
      phi3 <- ((SRCr - SRCnr)/(j))/((SRCnr)/(T - k))
      
      #Valor crítico al 95%
      if (phi3 > 6.49) {
         print("DGP2 | No hay RU")
      } else {
         print("DGP4 | Hay RU")
      }
      
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1,contador, op_back(yt,1)))
      m <- col_aum(yt, p_optimo)
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #restringida
      x_restringido <- cbind(x[,1], m)
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_opt),])
      y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
      
      #parametros
      j <- 2
      T <- nrow(x)
      k <- ncol(x)
      
      #regresión
      coef_nr <- regresion(x,y)
      
      #No restringido
      SRCnr <- sum((y - (x %*% coef_nr))^2)
      
      #Restringido
      x_restringido <- as.matrix(x_restringido[-(1:p_opt),])
      coef_r <- regresion(x_restringido, y)
      SRCr <- sum((y - (x_restringido %*% coef_r))^2)
      
      #Phi3
      phi3 <- ((SRCr - SRCnr)/(j))/((SRCnr)/(T - k))
      
      #Valor crítico al 95%
      if (phi3 > 6.49) {
         print("DGP2 | No hay RU")
      } else {
         print("DGP4 | Hay RU")
      }
   }
}

#Paso 3 - E3
Paso3_E3 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1, contador))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #regresion
      coef <- regresion(x,y)
      beta <- coef[2]
      
      #Estadístico t
      t <- beta / sd(contador)
      vc <- abs(qt(1 - (0.05/2), length(y) - 1))
      
      if (t > vc) {
         print("Hay RU")
      } else {
         print("No hay RU")
      }
      
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1, contador))
      m <- col_aum(yt, p_optimo)
      
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_opt),])
      y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
      
      #regresion
      coef <- regresion(x,y)
      beta <- coef[2]
      
      #Estadístico t
      t <- beta / sd(contador)
      vc <- abs(qt(1 - (0.05/2), length(y) - 1))
      
      if (t > vc) {
         print("Hay RU")
      } else {
         print("No hay RU")
      }
   }
}

#Paso 1 - E2
Paso1_E2 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[2,1]
      
      #error
      e_est <- y - (x %*% coef)
      
      #vcv
      #parametros iniciales
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      
      #tau_miu
      tm <- gamma / sqrt(vcv[k,k])
      
      if (tm < -2.90) {
         print("No hay RU")
      } else {
         print("Probar bajo phi1")
      }
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      m <- col_aum(yt, p_optimo)
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_opt),])
      y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[2,1]
      
      #error
      e_est <- y - (x %*% coef)
      
      #vcv
      #parametros iniciales
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      
      #tau_miu
      tm <- gamma / sqrt(vcv[k,k])
      
      if (tm < -2.90) {
         print("No hay RU")
      } else {
         print("Probar bajo phi1")
      }
   }
}

#Paso 2 - E2
Paso2_E2 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #parametros
      j <- 2
      T <- nrow(x)
      k <- ncol(x)
      
      #regresion
      coef_nr <- regresion(x,y)
      
      #No restringido
      SRCnr <- sum((y - (x %*% coef_nr))^2)
      
      #restringido
      SRCr <- sum(y^2)
      
      #phi1
      phi1 <- ((SRCr - SRCnr)/(j))/((SRCnr)/(T-k))
      
      #Valor critico al 95%
      if (phi1>4.71) {
         print("Rechazo | Probar que alpha es igual a 0")
      } else {
         print("No rechazo | Ir al esquema 1")
      }
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      m <- col_aum(yt, p_optimo)
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_opt),])
      y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
      
      #parametros
      j <- 2
      T <- nrow(x)
      k <- ncol(x)
      
      #regresion
      coef_nr <- regresion(x,y)
      
      #No restringido
      SRCnr <- sum((y - (x %*% coef_nr))^2)
      
      #restringido
      x_restringido <- as.matrix(m[-(1:p_opt),])
      coef_r <- regresion(x_restringido, y)
      SRCr <- sum((y - (x_restringido %*% coef_r))^2)
      
      #phi1
      phi1 <- ((SRCr - SRCnr)/(j))/((SRCnr)/(T-k))
      
      #Valor critico al 95%
      if (phi1>4.71) {
         print("Rechazo | Probar que alpha es igual a 0")
      } else {
         print("No rechazo | Ir al esquema 1")
      }
   }
}

#Paso 3 - E3
Paso3_E2 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #regresion
      coef <- regresion(x,y)
      alpha = coef[1,1]
      
      #estadistico
      e_est <- y - (x %*% coef)
      sgm_e <- as.numeric((t(e_est)%*%e_est)/(nrow(x))-ncol(x))
      vcv <- sgm_e * solve(t(x) %*% x)
      std_alpha <- vcv[1,1]
      est <- alpha/std_alpha
      
      #prueba t
      vc <- qt((1-0.05/2), (length(yt) - nrow(coef)))
      
      if (abs(est) > abs(vc)) {
         print("Rechazo | No hay RU")
      } else {
         print("No rechazo | probar RU en el esquema 1")
      }
      
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      m <- col_aum(yt, p_optimo)
      
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_opt),])
      y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
      
      #regresion
      coef <- regresion(x,y)
      alpha = coef[1,1]
      
      #estadistico
      e_est <- y - (x %*% coef)
      sgm_e <- as.numeric((t(e_est)%*%e_est)/(nrow(x))-ncol(x))
      vcv <- sgm_e * solve(t(x) %*% x)
      std_alpha <- vcv[1,1]
      est <- alpha/std_alpha
      
      #prueba t
      vc <- qt((1-0.05/2), (length(yt) - nrow(coef)))
      
      if (abs(est) > abs(vc)) {
         print("Rechazo | No hay RU")
      } else {
         print("No rechazo | probar RU en el esquema 1")
      } 
   }
}

#Paso 1 - E1
Paso1_E1 <- function(yt, p_optimo){
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(op_back(yt,1))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[1,1]
      
      #error
      e_est <- y - (x %*% coef)
      
      #vcv
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      
      #tao
      t <- gamma / sqrt(vcv[k,k])
      
      if (t < -1.95) {
         print("No hay RU")
      } else {
         print("DGP3 | Hay RU")
      }
      
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(op_back(yt,1))
      m <- col_aum(yt, p_optimo)
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_optimo),])
      y <- as.matrix(delta_yt(yt)[-(1:p_optimo),])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[1,1]
      
      #error
      e_est <- y - (x %*% coef)
      
      #vcv
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      
      #tao
      t <- gamma / sqrt(vcv[k,k])
      
      if (t < -1.95) {
         print("No hay RU")
      } else {
         print("DGP3 | Hay RU")
      }
   }
}


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


#Importación de datos----------
library(readxl)
#Para anÁlisis de cointegración por ARDL
library(ARDL)
data <- read.csv(file.choose())
#Escoger la columna donde están los rendimientos de bonos a 10 años
r10 <- as.matrix(data[,2], ncol = 1)
#Escoger la columna donde están los rendimientos de bonos a 3 meses
r3 <- as.matrix(data[,3], ncol = 1)




#CÓDIGO DICKEY FULLEN PARA COMPROBAR EXISTENCIA DE RAÍZ UNITARIA EN LA VARIABLE DEPENDIENTE-----

#########COMENTARIOS ADF:
#Para calibrar cada p óptimo de cada esquema es necesario correrlo linea a linea el código aparece enunciado en cada paso indicando a qué esquema y paso seguir, sin embargo solo es necesario correr los if() correspondientes para obtener el resultado

#Esquema 3
#Probar el P óptimo para el esquema (es necesario correr la funcion linea a linea, en cambio de correr directamente la función)
#r10
P_opt_E3 <- function(yt) {
   #Pasar datos a Yt y revisar length/nrow
   yt <- as.matrix(r10)
   #Y diferenciado
   y <- as.matrix(delta_yt(yt))
   
   #xdiseño
   contador <- 1:length(yt)
   col_fijas <- as.matrix(cbind(1,contador, op_back(yt,1)))
   
   #eliminar NAs iniciales
   y <- as.matrix(y[-1,])
   x <- as.matrix(col_fijas[-1,])
   
   #regresion
   coef <- regresion(x,y)
   gamma <- coef[3,1]
   
   #Error estimado
   e_est <- y - (x %*% coef)
   
   p_opt <- 1
   
   #Ruido blanco
   RB <- rb(e_est, p_opt)
   
   #Parsimonia
   Parsimonia <- parsimonia(e_est, x)
   
   #¿DF o ADF?
   if ( RB & Parsimonia == TRUE) {
      p_opt <- 1
      rta_p1 <- TRUE
   } else {
      p_opt <- 2
      rta_p1 <- FALSE
   }
   
   #P óptimo en ADF
   if (rta_p1 == FALSE) {
      while(TRUE) {
         #columnas ADF
         m <- col_aum(yt, p_opt)
         #matriz completa
         m_completa <- as.matrix(cbind(col_fijas,m))
         #quitar NAs iniciales
         x <- as.matrix(m_completa[-(1:p_opt),])
         y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
         #regresión
         coef <- regresion(x,y)
         gamma <- coef[nrow(coef), 1]
         #error
         e_est <- y - x %*% coef
         #ruido blanco
         RB <- rb(e_est, p_opt)
         #parsimonia
         Parsimonia <- parsimonia(e_est, x)
         #P óptimo
         if (RB & Parsimonia == TRUE) {
            p_optimo <- p_opt
            rta <- print(paste("El P óptimo bajo el esquema 3 es:", p_optimo))
            break
         } else {
            p_opt <- p_opt + 1
         }
      }
   } else {
      p_optimo <- p_opt
      rta <- print(paste("El P óptimo bajo el esquema 3 es:", p_optimo))
   }
   return(rta)
}

#1) El primer paso del E3, es probar gamma = 0 bajo taotao
Paso1_E3 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1,contador, op_back(yt,1)))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[3,1]
      
      #vcv
      #parametros iniciales
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      #taotao
      tt <- gamma / sqrt(vcv[k,k])
      tt_5 <- -3.41
      if (tt < tt_5) {
         print("Es un DGP1 | No hay RU")
      } else {
         print("Es un DGP3 | Hay RU")
      }
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1,contador, op_back(yt,1)))
      m <- col_aum(yt, p_optimo)
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #eliminar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_optimo),])
      y <- as.matrix(delta_yt(yt)[-(1:p_optimo),])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[3,1]
      
      #vcv
      #parametros iniciales
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      #taotao
      tt <- gamma / sqrt(vcv[k,k])
      tt_5 <- -3.41
      if (tt < tt_5) {
         print("Es un DGP1 | No hay RU")
      } else {
         print("Es un DGP3 | Hay RU")
      }
   }
}
#Rta
if (Paso1_E3(yt, p_optimo) == "Es un DGP1 | No hay RU") {
   print("No hay RU, por lo tanto el DF para acá")
} else {
   print("Ir al paso 2, en el que se probará gamma = beta = 0")
}

#2) Usando phi3 se probará gamma = beta = 0 
Paso2_E3 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1,contador, op_back(yt,1)))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #parametros
      j <- 2
      T <- nrow(x)
      k <- ncol(x)
      
      #regresion
      coef_nr <- regresion(x,y)
      
      #No restringido
      SRCnr <- sum((y - (x %*% coef_nr))^2)
      
      #Restringido
      x_restringido <- x[,1]
      coef_r <- regresion(x_restringido, y)
      SRCr <- sum((y - (x_restringido %*% coef_r))^2)
      
      #Phi3
      phi3 <- ((SRCr - SRCnr)/(j))/((SRCnr)/(T - k))
      
      #Valor crítico al 95%
      if (phi3 > 6.49) {
         print("DGP2 | No hay RU")
      } else {
         print("DGP4 | Hay RU")
      }
      
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1,contador, op_back(yt,1)))
      m <- col_aum(yt, p_optimo)
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #restringida
      x_restringido <- cbind(x[,1], m)
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_opt),])
      y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
      
      #parametros
      j <- 2
      T <- nrow(x)
      k <- ncol(x)
      
      #regresión
      coef_nr <- regresion(x,y)
      
      #No restringido
      SRCnr <- sum((y - (x %*% coef_nr))^2)
      
      #Restringido
      x_restringido <- as.matrix(x_restringido[-(1:p_opt),])
      coef_r <- regresion(x_restringido, y)
      SRCr <- sum((y - (x_restringido %*% coef_r))^2)
      
      #Phi3
      phi3 <- ((SRCr - SRCnr)/(j))/((SRCnr)/(T - k))
      
      #Valor crítico al 95%
      if (phi3 > 6.49) {
         print("DGP2 | No hay RU")
      } else {
         print("DGP4 | Hay RU")
      }
   }
}
#Rta
if (Paso2_E3(yt, p_optimo) == "DGP2 | No hay RU") {
   print("Ir al tercer paso, para probar beta = 0")
} else {
   print("Probar RU en el esquema 2")
}

#3)Si en el paso 2 rechazó, i.e. "DGP2 | No hay RU", se procederá a probar beta = 0 
Paso3_E3 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1, contador))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #regresion
      coef <- regresion(x,y)
      beta <- coef[2]
      
      #Estadístico t
      t <- beta / sd(contador)
      vc <- abs(qt(1 - (0.05/2), length(y) - 1))
      
      if (t > vc) {
         print("Hay RU")
      } else {
         print("No hay RU")
      }
      
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1, contador))
      m <- col_aum(yt, p_optimo)
      
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_opt),])
      y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
      
      #regresion
      coef <- regresion(x,y)
      beta <- coef[2]
      
      #Estadístico t
      t <- beta / sd(contador)
      vc <- abs(qt(1 - (0.05/2), length(y) - 1))
      
      if (t > vc) {
         print("Hay RU")
      } else {
         print("No hay RU")
      }
   }
}
#Rta
if (Paso3_E3(yt, p_optimo) == "No hay RU") {
   print("Puedo concluir que no hay RU bajo el esquema 3")
} else {
   print("Probar RU en el esquema 2")
}

#Pasamos al esquema 2, por lo que es necesario estimar nuevamente el p óptimo. Al igual que en el esquema anterior, es necesario correrlo linea a linea
P_opt_E2 <- function(yt) {
   #Pasar datos a Yt y revisar length/nrow
   yt <- as.matrix(yt)
   #Y diferenciado
   y <- as.matrix(delta_yt(yt))
   
   #xdiseño
   col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
   
   #eliminar NAs iniciales
   y <- as.matrix(y[-1,])
   x <- as.matrix(col_fijas[-1,])
   
   #regresion
   coef <- regresion(x,y)
   gamma <- coef[2,1]
   
   #Error estimado
   e_est <- y - (x %*% coef)
   
   p_opt <- 1
   
   #Ruido blanco
   RB <- rb(e_est, p_opt)
   
   #Parsimonia
   Parsimonia <- parsimonia(e_est, x)
   
   #¿DF o ADF?
   if ( RB & Parsimonia == TRUE) {
      p_opt <- 1
      rta_p1 <- TRUE
   } else {
      p_opt <- 2
      rta_p1 <- FALSE
   }
   
   #P óptimo en ADF
   if (rta_p1 == FALSE) {
      while(TRUE) {
         #columnas ADF
         m <- col_aum(yt, p_opt)
         #matriz completa
         m_completa <- as.matrix(cbind(col_fijas,m))
         #quitar NAs iniciales
         x <- as.matrix(m_completa[-(1:p_opt),])
         y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
         #regresión
         coef <- regresion(x,y)
         gamma <- coef[nrow(coef), 1]
         #error
         e_est <- y - x %*% coef
         #ruido blanco
         RB <- rb(e_est, p_opt)
         #parsimonia
         Parsimonia <- parsimonia(e_est, x)
         #P óptimo
         if (RB & Parsimonia == TRUE) {
            p_optimo <- p_opt
            rta <- print(paste("El P óptimo bajo el esquema 2 es:", p_optimo))
            break
         } else {
            p_opt <- p_opt + 1
         }
      }
   } else {
      p_optimo <- p_opt
      rta <- print(paste("El P óptimo bajo el esquema 2 es:", p_optimo))
   }
   return(rta)
}

#4)El primer paso del esquema 2 es probar gamma = 0 usando tao miu
Paso1_E2 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[2,1]
      
      #error
      e_est <- y - (x %*% coef)
      
      #vcv
      #parametros iniciales
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      
      #tau_miu
      tm <- gamma / sqrt(vcv[k,k])
      
      if (tm < -2.90) {
         print("No hay RU")
      } else {
         print("Probar bajo phi1")
      }
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      m <- col_aum(yt, p_optimo)
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_opt),])
      y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[2,1]
      
      #error
      e_est <- y - (x %*% coef)
      
      #vcv
      #parametros iniciales
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      
      #tau_miu
      tm <- gamma / sqrt(vcv[k,k])
      
      if (tm < -2.90) {
         print("No hay RU")
      } else {
         print("Probar bajo phi1")
      }
   }
}
#Rta
if (Paso2_E3(yt, p_optimo) == "No hay RU") {
   print("Puedo concluir que bajo el esquema 2 no hay RU")
} else {
   print("Ir al paso 2 del esquema 2 a probar gamma = alpha = 0, bajo phi1")
}

#5)Se prueba que gamma = alpha = 0, bajo phi1
Paso2_E2 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #parametros
      j <- 2
      T <- nrow(x)
      k <- ncol(x)
      
      #regresion
      coef_nr <- regresion(x,y)
      
      #No restringido
      SRCnr <- sum((y - (x %*% coef_nr))^2)
      
      #restringido
      SRCr <- sum(y^2)
      
      #phi1
      phi1 <- ((SRCr - SRCnr)/(j))/((SRCnr)/(T-k))
      
      #Valor critico al 95%
      if (phi1>4.71) {
         print("Rechazo | Probar que alpha es igual a 0")
      } else {
         print("No rechazo | Ir al esquema 1")
      }
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      m <- col_aum(yt, p_optimo)
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_opt),])
      y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
      
      #parametros
      j <- 2
      T <- nrow(x)
      k <- ncol(x)
      
      #regresion
      coef_nr <- regresion(x,y)
      
      #No restringido
      SRCnr <- sum((y - (x %*% coef_nr))^2)
      
      #restringido
      x_restringido <- as.matrix(m[-(1:p_opt),])
      coef_r <- regresion(x_restringido, y)
      SRCr <- sum((y - (x_restringido %*% coef_r))^2)
      
      #phi1
      phi1 <- ((SRCr - SRCnr)/(j))/((SRCnr)/(T-k))
      
      #Valor critico al 95%
      if (phi1>4.71) {
         print("Rechazo | Probar que alpha es igual a 0")
      } else {
         print("No rechazo | Ir al esquema 1")
      }
   }
}
#Rta
if (Paso2_E2(yt, p_optimo) == "Rechazo | Probar que alpha es igual a 0"){
   print("Probar mediante una t que alpha = 0")
} else {
   print("Estimar nuevamente el p óptimo e ir al esquema 1 a probar RU")
}

#6)Si se rechaza que el Paso2 - E2 se debe probar bajo una t que alpha = 0
#Revisar
Paso3_E2 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #regresion
      coef <- regresion(x,y)
      alpha = coef[1,1]
      
      #estadistico
      e_est <- y - (x %*% coef)
      sgm_e <- as.numeric((t(e_est)%*%e_est)/(nrow(x))-ncol(x))
      vcv <- sgm_e * solve(t(x) %*% x)
      std_alpha <- vcv[1,1]
      est <- alpha/std_alpha
      
      #prueba t
      vc <- qt((1-0.05/2), (length(yt) - nrow(coef)))
      
      if (abs(est) > abs(vc)) {
         print("Rechazo | No hay RU")
      } else {
         print("No rechazo | probar RU en el esquema 1")
      }
      
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      m <- col_aum(yt, p_optimo)
      
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_opt),])
      y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
      
      #regresion
      coef <- regresion(x,y)
      alpha = coef[1,1]
      
      #estadistico
      e_est <- y - (x %*% coef)
      sgm_e <- as.numeric((t(e_est)%*%e_est)/(nrow(x))-ncol(x))
      vcv <- sgm_e * solve(t(x) %*% x)
      std_alpha <- vcv[1,1]
      est <- alpha/std_alpha
      
      #prueba t
      vc <- qt((1-0.05/2), (length(yt) - nrow(coef)))
      
      if (abs(est) > abs(vc)) {
         print("Rechazo | No hay RU")
      } else {
         print("No rechazo | probar RU en el esquema 1")
      } 
   }
}
#Rta
if (Paso3_E2(yt, p_optimo) == "Rechazo | No hay RU") {
   print("No hay RU bajo el esquema 2")
} else {
   print("Estimar nuevamente el p óptimo e ir al esquema 1 a probar RU")
}

#Pasamos al último esquema (1), por lo que debemos estimar nuevamente el p óptimo
P_opt_E1 <- function(yt) {
   #Pasar datos a Yt y revisar length/nrow
   yt <- as.matrix(yt)
   #Y diferenciado
   y <- as.matrix(delta_yt(yt))
   
   #xdiseño
   col_fijas <- as.matrix(cbind(op_back(yt,1)))
   
   #eliminar NAs iniciales
   y <- as.matrix(y[-1,])
   x <- as.matrix(col_fijas[-1,])
   
   #regresion
   coef <- regresion(x,y)
   gamma <- coef[1,1]
   
   #Error estimado
   e_est <- y - (x %*% coef)
   
   p_opt <- 1
   
   #Ruido blanco
   RB <- rb(e_est, p_opt)
   
   #Parsimonia
   Parsimonia <- parsimonia(e_est, x)
   
   #¿DF o ADF?
   if ( RB & Parsimonia == TRUE) {
      p_opt <- 1
      rta_p1 <- TRUE
   } else {
      p_opt <- 2
      rta_p1 <- FALSE
   }
   
   #P óptimo en ADF
   if (rta_p1 == FALSE) {
      while(TRUE) {
         #columnas ADF
         m <- col_aum(yt, p_opt)
         #matriz completa
         m_completa <- as.matrix(cbind(col_fijas,m))
         #quitar NAs iniciales
         x <- as.matrix(m_completa[-(1:p_opt),])
         y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
         #regresión
         coef <- regresion(x,y)
         gamma <- coef[nrow(coef), 1]
         #error
         e_est <- y - x %*% coef
         #ruido blanco
         RB <- rb(e_est, p_opt)
         #parsimonia
         Parsimonia <- parsimonia(e_est, x)
         #P óptimo
         if (RB & Parsimonia == TRUE) {
            p_optimo <- p_opt
            rta <- print(paste("El P óptimo bajo el esquema 1 es:", p_optimo))
            break
         } else {
            p_opt <- p_opt + 1
         }
      }
   } else {
      p_optimo <- p_opt
      rta <- print(paste("El P óptimo bajo el esquema 1 es:", p_optimo))
   }
   return(rta)
}

#7)Bajo tao, probaremos que gamma = 0 en el esquema 1
Paso1_E1 <- function(yt, p_optimo){
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(op_back(yt,1))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[1,1]
      
      #error
      e_est <- y - (x %*% coef)
      
      #vcv
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      
      #tao
      t <- gamma / sqrt(vcv[k,k])
      
      if (t < -1.95) {
         print("No hay RU")
      } else {
         print("DGP3 | Hay RU")
      }
      
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(op_back(yt,1))
      m <- col_aum(yt, p_optimo)
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_optimo),])
      y <- as.matrix(delta_yt(yt)[-(1:p_optimo),])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[1,1]
      
      #error
      e_est <- y - (x %*% coef)
      
      #vcv
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      
      #tao
      t <- gamma / sqrt(vcv[k,k])
      
      if (t < -1.95) {
         print("No hay RU")
      } else {
         print("DGP3 | Hay RU")
      }
   }
}
#Rta
if (Paso1_E1(yt, p_optimo) == "No hay RU") {
   print("Bajo el esquema 1 no hay RU")
} else {
   print("Bajo el esquema 1 hay RU")
}

#Para r10 hay raiz unitaria en el esquema 1
#Ahora se probará para r3



#MODELO-----
#Se realizan transformaciones a la variable de bonos a 3 meses para crear la X que se adapta al modelo
#Sacar la media móvil para las expectativas (en el documento se explicará)
expr3 <- MA(r3, 40)
#Teniendo en cuenta el modelo de la expectativas puras, se realiza la siguiente modificación 
X <- 0.5*(r3[40:247] + expr3)


#COMPROBACIÓN DICKEY FULLEN PARA LA VARIABLE EXPLICATIVA----

#Esquema 3

#########COMENTARIOS ADF:
#Para calibrar cada p óptimo de cada esquema es necesario correrlo linea a linea el código aparece enunciado en cada paso indicando a qué esquema y paso seguir, sin embargo solo es necesario correr los if() correspondientes para obtener el resultado


#Probar el P óptimo para el esquema (es necesario correr la funcion linea a linea, en cambio de correr directamente la función)
P_opt_E3 <- function(x) {
   #Pasar datos a Yt y revisar length/nrow
   yt <- as.matrix(X)
   #Y diferenciado
   y <- as.matrix(delta_yt(yt))
   
   #xdiseño
   contador <- 1:length(yt)
   col_fijas <- as.matrix(cbind(1,contador, op_back(yt,1)))
   
   #eliminar NAs iniciales
   y <- as.matrix(y[-1,])
   x <- as.matrix(col_fijas[-1,])
   
   #regresion
   coef <- regresion(x,y)
   gamma <- coef[3,1]
   
   #Error estimado
   e_est <- y - (x %*% coef)
   
   p_opt <- 1
   
   #Ruido blanco
   RB <- rb(e_est, p_opt)
   
   #Parsimonia
   Parsimonia <- parsimonia(e_est, x)
   
   #¿DF o ADF?
   if ( RB & Parsimonia == TRUE) {
      p_opt <- 1
      rta_p1 <- TRUE
   } else {
      p_opt <- 2
      rta_p1 <- FALSE
   }
   
   #P óptimo en ADF
   if (rta_p1 == FALSE) {
      while(TRUE) {
         #columnas ADF
         m <- col_aum(yt, p_opt)
         #matriz completa
         m_completa <- as.matrix(cbind(col_fijas,m))
         #quitar NAs iniciales
         x <- as.matrix(m_completa[-(1:p_opt),])
         y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
         #regresión
         coef <- regresion(x,y)
         gamma <- coef[nrow(coef), 1]
         #error
         e_est <- y - x %*% coef
         #ruido blanco
         RB <- rb(e_est, p_opt)
         #parsimonia
         Parsimonia <- parsimonia(e_est, x)
         #P óptimo
         if (RB & Parsimonia == TRUE) {
            p_optimo <- p_opt
            rta <- print(paste("El P óptimo bajo el esquema 3 es:", p_optimo))
            break
         } else {
            p_opt <- p_opt + 1
         }
      }
   } else {
      p_optimo <- p_opt
      rta <- print(paste("El P óptimo bajo el esquema 3 es:", p_optimo))
   }
   return(rta)
}

#1) El primer paso del E3, es probar gamma = 0 bajo taotao
Paso1_E3 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1,contador, op_back(yt,1)))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[3,1]
      
      #vcv
      #parametros iniciales
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      #taotao
      tt <- gamma / sqrt(vcv[k,k])
      tt_5 <- -3.41
      if (tt < tt_5) {
         print("Es un DGP1 | No hay RU")
      } else {
         print("Es un DGP3 | Hay RU")
      }
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1,contador, op_back(yt,1)))
      m <- col_aum(yt, p_optimo)
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #eliminar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_optimo),])
      y <- as.matrix(delta_yt(yt)[-(1:p_optimo),])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[3,1]
      
      #vcv
      #parametros iniciales
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      #taotao
      tt <- gamma / sqrt(vcv[k,k])
      tt_5 <- -3.41
      if (tt < tt_5) {
         print("Es un DGP1 | No hay RU")
      } else {
         print("Es un DGP3 | Hay RU")
      }
   }
}
#Rta
if (Paso1_E3(yt, p_optimo) == "Es un DGP1 | No hay RU") {
   print("No hay RU, por lo tanto el DF para acá")
} else {
   print("Ir al paso 2, en el que se probará gamma = beta = 0")
}

#2) Usando phi3 se probará gamma = beta = 0 
Paso2_E3 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1,contador, op_back(yt,1)))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #parametros
      j <- 2
      T <- nrow(x)
      k <- ncol(x)
      
      #regresion
      coef_nr <- regresion(x,y)
      
      #No restringido
      SRCnr <- sum((y - (x %*% coef_nr))^2)
      
      #Restringido
      x_restringido <- x[,1]
      coef_r <- regresion(x_restringido, y)
      SRCr <- sum((y - (x_restringido %*% coef_r))^2)
      
      #Phi3
      phi3 <- ((SRCr - SRCnr)/(j))/((SRCnr)/(T - k))
      
      #Valor crítico al 95%
      if (phi3 > 6.49) {
         print("DGP2 | No hay RU")
      } else {
         print("DGP4 | Hay RU")
      }
      
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1,contador, op_back(yt,1)))
      m <- col_aum(yt, p_optimo)
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #restringida
      x_restringido <- cbind(x[,1], m)
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_opt),])
      y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
      
      #parametros
      j <- 2
      T <- nrow(x)
      k <- ncol(x)
      
      #regresión
      coef_nr <- regresion(x,y)
      
      #No restringido
      SRCnr <- sum((y - (x %*% coef_nr))^2)
      
      #Restringido
      x_restringido <- as.matrix(x_restringido[-(1:p_opt),])
      coef_r <- regresion(x_restringido, y)
      SRCr <- sum((y - (x_restringido %*% coef_r))^2)
      
      #Phi3
      phi3 <- ((SRCr - SRCnr)/(j))/((SRCnr)/(T - k))
      
      #Valor crítico al 95%
      if (phi3 > 6.49) {
         print("DGP2 | No hay RU")
      } else {
         print("DGP4 | Hay RU")
      }
   }
}
#Rta
if (Paso2_E3(yt, p_optimo) == "DGP2 | No hay RU") {
   print("Ir al tercer paso, para probar beta = 0")
} else {
   print("Probar RU en el esquema 2")
}

#3)Si en el paso 2 rechazó, i.e. "DGP2 | No hay RU", se procederá a probar beta = 0 
Paso3_E3 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1, contador))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #regresion
      coef <- regresion(x,y)
      beta <- coef[2]
      
      #Estadístico t
      t <- beta / sd(contador)
      vc <- abs(qt(1 - (0.05/2), length(y) - 1))
      
      if (t > vc) {
         print("Hay RU")
      } else {
         print("No hay RU")
      }
      
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      contador <- 1:length(yt)
      col_fijas <- as.matrix(cbind(1, contador))
      m <- col_aum(yt, p_optimo)
      
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_opt),])
      y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
      
      #regresion
      coef <- regresion(x,y)
      beta <- coef[2]
      
      #Estadístico t
      t <- beta / sd(contador)
      vc <- abs(qt(1 - (0.05/2), length(y) - 1))
      
      if (t > vc) {
         print("Hay RU")
      } else {
         print("No hay RU")
      }
   }
}
#Rta
if (Paso3_E3(yt, p_optimo) == "No hay RU") {
   print("Puedo concluir que no hay RU bajo el esquema 3")
} else {
   print("Probar RU en el esquema 2")
}

#Pasamos al esquema 2, por lo que es necesario estimar nuevamente el p óptimo. Al igual que en el esquema anterior, es necesario correrlo linea a linea
P_opt_E2 <- function(yt) {
   #Pasar datos a Yt y revisar length/nrow
   yt <- as.matrix(X)
   #Y diferenciado
   y <- as.matrix(delta_yt(yt))
   
   #xdiseño
   col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
   
   #eliminar NAs iniciales
   y <- as.matrix(y[-1,])
   x <- as.matrix(col_fijas[-1,])
   
   #regresion
   coef <- regresion(x,y)
   gamma <- coef[2,1]
   
   #Error estimado
   e_est <- y - (x %*% coef)
   
   p_opt <- 1
   
   #Ruido blanco
   RB <- rb(e_est, p_opt)
   
   #Parsimonia
   Parsimonia <- parsimonia(e_est, x)
   
   #¿DF o ADF?
   if ( RB & Parsimonia == TRUE) {
      p_opt <- 1
      rta_p1 <- TRUE
   } else {
      p_opt <- 2
      rta_p1 <- FALSE
   }
   
   #P óptimo en ADF
   if (rta_p1 == FALSE) {
      while(TRUE) {
         #columnas ADF
         m <- col_aum(yt, p_opt)
         #matriz completa
         m_completa <- as.matrix(cbind(col_fijas,m))
         #quitar NAs iniciales
         x <- as.matrix(m_completa[-(1:p_opt),])
         y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
         #regresión
         coef <- regresion(x,y)
         gamma <- coef[nrow(coef), 1]
         #error
         e_est <- y - x %*% coef
         #ruido blanco
         RB <- rb(e_est, p_opt)
         #parsimonia
         Parsimonia <- parsimonia(e_est, x)
         #P óptimo
         if (RB & Parsimonia == TRUE) {
            p_optimo <- p_opt
            rta <- print(paste("El P óptimo bajo el esquema 2 es:", p_optimo))
            break
         } else {
            p_opt <- p_opt + 1
         }
      }
   } else {
      p_optimo <- p_opt
      rta <- print(paste("El P óptimo bajo el esquema 2 es:", p_optimo))
   }
   return(rta)
}

#4)El primer paso del esquema 2 es probar gamma = 0 usando tao miu
Paso1_E2 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[2,1]
      
      #error
      e_est <- y - (x %*% coef)
      
      #vcv
      #parametros iniciales
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      
      #tau_miu
      tm <- gamma / sqrt(vcv[k,k])
      
      if (tm < -2.90) {
         print("No hay RU")
      } else {
         print("Probar bajo phi1")
      }
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      m <- col_aum(yt, p_optimo)
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_opt),])
      y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[2,1]
      
      #error
      e_est <- y - (x %*% coef)
      
      #vcv
      #parametros iniciales
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      
      #tau_miu
      tm <- gamma / sqrt(vcv[k,k])
      
      if (tm < -2.90) {
         print("No hay RU")
      } else {
         print("Probar bajo phi1")
      }
   }
}
#Rta
if (Paso2_E3(yt, p_optimo) == "No hay RU") {
   print("Puedo concluir que bajo el esquema 2 no hay RU")
} else {
   print("Ir al paso 2 del esquema 2 a probar gamma = alpha = 0, bajo phi1")
}

#5)Se prueba que gamma = alpha = 0, bajo phi1
Paso2_E2 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #parametros
      j <- 2
      T <- nrow(x)
      k <- ncol(x)
      
      #regresion
      coef_nr <- regresion(x,y)
      
      #No restringido
      SRCnr <- sum((y - (x %*% coef_nr))^2)
      
      #restringido
      SRCr <- sum(y^2)
      
      #phi1
      phi1 <- ((SRCr - SRCnr)/(j))/((SRCnr)/(T-k))
      
      #Valor critico al 95%
      if (phi1>4.71) {
         print("Rechazo | Probar que alpha es igual a 0")
      } else {
         print("No rechazo | Ir al esquema 1")
      }
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      m <- col_aum(yt, p_optimo)
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_opt),])
      y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
      
      #parametros
      j <- 2
      T <- nrow(x)
      k <- ncol(x)
      
      #regresion
      coef_nr <- regresion(x,y)
      
      #No restringido
      SRCnr <- sum((y - (x %*% coef_nr))^2)
      
      #restringido
      x_restringido <- as.matrix(m[-(1:p_opt),])
      coef_r <- regresion(x_restringido, y)
      SRCr <- sum((y - (x_restringido %*% coef_r))^2)
      
      #phi1
      phi1 <- ((SRCr - SRCnr)/(j))/((SRCnr)/(T-k))
      
      #Valor critico al 95%
      if (phi1>4.71) {
         print("Rechazo | Probar que alpha es igual a 0")
      } else {
         print("No rechazo | Ir al esquema 1")
      }
   }
}
#Rta
if (Paso2_E2(yt, p_optimo) == "Rechazo | Probar que alpha es igual a 0"){
   print("Probar mediante una t que alpha = 0")
} else {
   print("Estimar nuevamente el p óptimo e ir al esquema 1 a probar RU")
}

#6)Si se rechaza que el Paso2 - E2 se debe probar bajo una t que alpha = 0
#Revisar
Paso3_E2 <- function(yt, p_optimo) {
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #regresion
      coef <- regresion(x,y)
      alpha = coef[1,1]
      
      #estadistico
      e_est <- y - (x %*% coef)
      sgm_e <- as.numeric((t(e_est)%*%e_est)/(nrow(x))-ncol(x))
      vcv <- sgm_e * solve(t(x) %*% x)
      std_alpha <- vcv[1,1]
      est <- alpha/std_alpha
      
      #prueba t
      vc <- qt((1-0.05/2), (length(yt) - nrow(coef)))
      
      if (abs(est) > abs(vc)) {
         print("Rechazo | No hay RU")
      } else {
         print("No rechazo | probar RU en el esquema 1")
      }
      
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(cbind(1, op_back(yt,1)))
      m <- col_aum(yt, p_optimo)
      
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_opt),])
      y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
      
      #regresion
      coef <- regresion(x,y)
      alpha = coef[1,1]
      
      #estadistico
      e_est <- y - (x %*% coef)
      sgm_e <- as.numeric((t(e_est)%*%e_est)/(nrow(x))-ncol(x))
      vcv <- sgm_e * solve(t(x) %*% x)
      std_alpha <- vcv[1,1]
      est <- alpha/std_alpha
      
      #prueba t
      vc <- qt((1-0.05/2), (length(yt) - nrow(coef)))
      
      if (abs(est) > abs(vc)) {
         print("Rechazo | No hay RU")
      } else {
         print("No rechazo | probar RU en el esquema 1")
      } 
   }
}
#Rta
if (Paso3_E2(yt, p_optimo) == "Rechazo | No hay RU") {
   print("No hay RU bajo el esquema 2")
} else {
   print("Estimar nuevamente el p óptimo e ir al esquema 1 a probar RU")
}

#Pasamos al último esquema (1), por lo que debemos estimar nuevamente el p óptimo
P_opt_E1 <- function(yt) {
   #Pasar datos a Yt y revisar length/nrow
   yt <- as.matrix(X)
   #Y diferenciado
   y <- as.matrix(delta_yt(yt))
   
   #xdiseño
   col_fijas <- as.matrix(cbind(op_back(yt,1)))
   
   #eliminar NAs iniciales
   y <- as.matrix(y[-1,])
   x <- as.matrix(col_fijas[-1,])
   
   #regresion
   coef <- regresion(x,y)
   gamma <- coef[1,1]
   
   #Error estimado
   e_est <- y - (x %*% coef)
   
   p_opt <- 1
   
   #Ruido blanco
   RB <- rb(e_est, p_opt)
   
   #Parsimonia
   Parsimonia <- parsimonia(e_est, x)
   
   #¿DF o ADF?
   if ( RB & Parsimonia == TRUE) {
      p_opt <- 1
      rta_p1 <- TRUE
   } else {
      p_opt <- 2
      rta_p1 <- FALSE
   }
   
   #P óptimo en ADF
   if (rta_p1 == FALSE) {
      while(TRUE) {
         #columnas ADF
         m <- col_aum(yt, p_opt)
         #matriz completa
         m_completa <- as.matrix(cbind(col_fijas,m))
         #quitar NAs iniciales
         x <- as.matrix(m_completa[-(1:p_opt),])
         y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
         #regresión
         coef <- regresion(x,y)
         gamma <- coef[nrow(coef), 1]
         #error
         e_est <- y - x %*% coef
         #ruido blanco
         RB <- rb(e_est, p_opt)
         #parsimonia
         Parsimonia <- parsimonia(e_est, x)
         #P óptimo
         if (RB & Parsimonia == TRUE) {
            p_optimo <- p_opt
            rta <- print(paste("El P óptimo bajo el esquema 1 es:", p_optimo))
            break
         } else {
            p_opt <- p_opt + 1
         }
      }
   } else {
      p_optimo <- p_opt
      rta <- print(paste("El P óptimo bajo el esquema 1 es:", p_optimo))
   }
   return(rta)
}

#7)Bajo tao, probaremos que gamma = 0 en el esquema 1
Paso1_E1 <- function(yt, p_optimo){
   if (p_optimo == 1) {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(op_back(yt,1))
      
      #eliminar NAs iniciales
      y <- as.matrix(y[-1,])
      x <- as.matrix(col_fijas[-1,])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[1,1]
      
      #error
      e_est <- y - (x %*% coef)
      
      #vcv
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      
      #tao
      t <- gamma / sqrt(vcv[k,k])
      
      if (t < -1.95) {
         print("No hay RU")
      } else {
         print("DGP3 | Hay RU")
      }
      
   } else {
      #Pasar datos a Yt y revisar length/nrow
      yt <- as.matrix(yt)
      #Y diferenciado
      y <- as.matrix(delta_yt(yt))
      
      #xdiseño
      col_fijas <- as.matrix(op_back(yt,1))
      m <- col_aum(yt, p_optimo)
      #matriz completa
      m_completa <- as.matrix(cbind(col_fijas,m))
      
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_optimo),])
      y <- as.matrix(delta_yt(yt)[-(1:p_optimo),])
      
      #regresion
      coef <- regresion(x,y)
      gamma <- coef[1,1]
      
      #error
      e_est <- y - (x %*% coef)
      
      #vcv
      T <- nrow(x)
      k <- ncol(x)
      #sigma_e
      sgm_e <- as.numeric((t(e_est) %*% e_est)/(T-k))
      vcv <- sgm_e * solve(t(x)%*%x)
      
      #tao
      t <- gamma / sqrt(vcv[k,k])
      
      if (t < -1.95) {
         print("No hay RU")
      } else {
         print("DGP3 | Hay RU")
      }
   }
}
#Rta
if (Paso1_E1(yt, p_optimo) == "No hay RU") {
   print("Bajo el esquema 1 no hay RU")
} else {
   print("Bajo el esquema 1 hay RU")
}

#Para r3 hay raiz unitaria en todos los esquemas



#Creación de matriz con solo las variables usadas-----
datos <- cbind(r10[40:247],X)
colnames(datos) <- c("r10", "X")

plot(r10[40:247], type = "l", lwd = 2, col = "#4682B4", ylim = c(0, 15), ylab = "Rendimientos", main = "Rendimientos de bonos a 10 años y 3 meses", xlab = "Trimestre")
lines(X, col = "#698B69", lwd = 2)
legend("topright", legend = c("Bonos 10 años", "Bonos a 3 meses"), col = c("#4682B4", "#698B69"), lty = c(1, 1))
#COINTEGRACIÓN POR ARDL------
#Cómo mi variable dependiente es I(1) puedo hacer ARDL
#creo el modelo, poniendo mi dependiente y mis explicativas
model <- auto_ardl(r10 ~ X, data = datos, max_order = 4)

#Se escoge el modelo con menor AIC (criterio de información) por los rezagos   
ardl_p <- model$best_model
#Se revisan los rezagos por cada variable 
ardl_p$order
#Se muestran los resultados
summary(ardl_p)


#Modelo de correción de error 
#Usamos el restringido para encontrar la velocidad de ajuste
RMCE <- recm(ardl_p, case = 3)
#Velocidad de ajuste 
vel_aju <- tail(RMCE$coefficients, 1)
summary(RMCE)
#se miran los coeficientes del MCE en el largo plazo

#Se revisa la existencia de cointegración
#Si hay I(0) se usa f statistic pero si son I(1) se usa t statistic, en este caso los t statistic
#Se tiene que especificar el esquema (case = 3)

bounds_t_test(MCE, case = 3, alpha = 0.01)
#SE PRUEBA QUE HAY COINTEGRACIÓN Y QUE EL MODELO ESTÁ CORRECTAMENTE ESPECIFICADO

#ANÁLISIS ENGLE-GRANGER-----
y <- r10[40:247]
datos <- cbind(y,X)
colnames(datos) <- c("r10","X")


P_opt <- function(x,y) {
   #Probar p=1 primero
   #xdiseño
   x <- cbind(1,X)
   #OLS
   coef <- regresion(x,y)
   #Error 
   #estimado
   e_est <- y - (x%*%coef)
   #diferenciado
   e_dif <- as.matrix(op_diferencias(e_est,1))
   #rezagado
   e_lag <- as.matrix(op_back(e_est,1))
   
   #quitar valores iniciales
   e_dif <- as.matrix(e_dif[-1,1])
   e_lag <- as.matrix(e_lag[-1,1])
   
   #gamma
   gamma <- regresion(e_lag,e_dif)
   
   #estimados
   e_dif_est <- e_lag %*% gamma
   neta_est <- e_dif - e_dif_est
   
   #probar que p=1 sea el óptimo, de lo contrario iterar hasta encontralo
   #RB
   RB <- rb(neta_est, 1)
   #parsimonia
   Parsimonia <- parsimonia_MCE(x,y,1)
   
   if ( RB & Parsimonia == TRUE) {
      p_opt <- 1
      rta_p1 <- TRUE
   } else {
      p_opt <- 2
      rta_p1 <- FALSE
   }
   rta <- (paste("El P óptimo para el modelo de corrección de error es:", p_opt))
   return(rta)
   #p optimo ampliado
   if(rta_p1 == FALSE) {
      while(TRUE) {
         #columnas aumentadas
         m <- col_aum(y, p_opt)
         #matriz completa
         m_completa <- cbind(x,m)
         #quitar NAs iniciales
         x <- as.matrix(m_completa[-(1:p_opt),])
         y <- as.matrix(y[-(1:p_opt),])
         #OLS
         coef <- regresion(x,y)
         #Error 
         #estimado
         e_est <- y - (x%*%coef)
         #diferenciado
         e_dif <- as.matrix(op_diferencias(e_est,1))
         #rezagado
         e_lag <- as.matrix(op_back(e_est,1))
         #quitar valores iniciales
         e_dif <- as.matrix(e_dif[-1,1])
         e_lag <- as.matrix(e_lag[-1,1])
         
         #gamma
         gamma <- regresion(e_lag,e_dif)
         
         #estimados
         e_dif_est <- e_lag %*% gamma
         neta_est <- e_dif - e_dif_est
         
         #RB
         RB <- rb(neta_est, p_opt)
         #parsimonia
         Parsimonia <- parsimonia_MCE(x,y,p_opt)
         if ( RB & Parsimonia == TRUE) {
            p_optimo <- p_opt
            rta_p1 <- TRUE
         } else {
            p_opt <- p_opt + 1
            rta_p1 <- FALSE
         }
         rta <- (paste("El P óptimo para el modelo de corrección de error es:", p_optimo))
         return(rta)
      }
   }
} 

#Se obtiene p_optimo igual a 1
Resultado_p_optimo <- P_opt(X,y)
p_optimo <- 1

CI <- function(x,y, p_optimo) { 
   if(p_optimo == 1) {
      #Agregar B0
      x <- cbind(1,x)
      #estimacion
      coef <- regresion(x,y)
      #error
      #estimado
      e_est <- as.matrix(y - (x%*%coef))
      #diferenciado
      e_dif <- as.matrix(op_diferencias(e_est,1))
      e_dif <- as.matrix(e_dif[-1,1])
      #lagg
      e_lag <- as.matrix(op_back(e_est,1))
      e_lag <- as.matrix(e_lag[-1,1])
      #coeficientes - regresion auxiliar
      coef_eaux <-regresion(e_lag,e_dif)
      gamma_aux <- coef_eaux[1]
      #error reg aux
      neta_aux <- e_dif - (e_lag %*% coef_eaux)
      #VCV
      sgme_aux <- as.numeric((t(neta_aux)%*%neta_aux)/(length(y) - ncol(e_lag)))
      vcv_aux <- sgme_aux * solve(t(e_lag)%*%e_lag)
      #estadistico
      t <- gamma_aux/sqrt(vcv_aux[1,1])
      #ph
      if (t<4.76) {
         respuesta_CI <- "Cointegradas | El termino de perturbación es I(0)"
         print(respuesta_CI)
      } else {
         respuesta_CI <- "No cointegradas | I(0)"
         print(respuesta_CI)
      }
   } else { 
      #Agregar B0
      x <- cbind(1,x)
      #columnas aumentadas
      m <- col_aum_MCE(y, p_optimo) ###Hay que ajustarlo manual###
      m_completa <- cbind(x,m)
      #quitar NAs iniciales
      x <- as.matrix(m_completa[-(1:p_optimo),])
      y <- as.matrix(y[-(1:p_optimo),])
      #estimacion
      coef <- regresion(x,y)
      #error
      #estimado
      e_est <- as.matrix(y - (x%*%coef))
      #diferenciado
      e_dif <- as.matrix(op_diferencias(e_est,1))
      e_dif <- as.matrix(e_dif[-1,1])
      #lagg
      e_lag <- as.matrix(op_back(e_est,1))
      e_lag <- as.matrix(e_lag[-1,1])
      #coeficientes - regresion auxiliar
      coef_eaux <-regresion(e_lag,e_dif)
      gamma_aux <- coef_eaux[1]
      #error reg aux
      neta_aux <- e_dif - (e_lag %*% coef_eaux)
      #VCV
      sgme_aux <- as.numeric((t(neta_aux)%*%neta_aux)/(length(y) - ncol(e_lag)))
      vcv_aux <- sgme_aux * solve(t(e_lag)%*%e_lag)
      #estadistico
      t <- gamma_aux/sqrt(vcv_aux[1,1])
      #ph
      if (t<4.76) {
         respuesta_CI <- "Cointegradas | I(0)"
         print(respuesta_CI)
      } else {
         respuesta_CI <- "No cointegradas | I(0)"
         print(respuesta_CI)
      }
   }
} 

rta_CI <- CI(X,y,p_optimo)
print(rta_CI)
#Hay cointegración entre las variables seleccionadas en la estimación del modelo MCE

MCE <- function(y,X,P_optimo) {
   #OLS
   coef <- regresion(X,y)
   #Estimar desequilibrio i.e. error estimado
   desequilibrio <- y - (X%*%coef)
   #desequilibrio lag
   desequilibrio_lag <- as.matrix(op_back(desequilibrio,1), ncol=1)
   desequilibrio_lag <- as.matrix(desequilibrio_lag[-1,1])
   desequilibrio_lag <- as.matrix(desequilibrio_lag[-1,1])
   #teorema representarción de Granger
   y_dif <- as.matrix(op_diferencias(y,1))
   y_dif <- as.matrix(y_dif[-1,1])
   y_dif <- as.matrix(y_dif[-1,1])
   #xdiseño
   X_1e <- as.matrix(op_back(X,1))
   X_1e <- as.matrix(X_1e[-1,1])
   X_1e <- as.matrix(X_1e[-1,1])
   ycp <- as.matrix(op_diferencias(op_back(y, 1),1))
   ycp <- as.matrix(ycp[-1,1])
   ycp <- as.matrix(ycp[-1,1])
   x_lag <- cbind(X_1e,ycp,desequilibrio_lag)
   #OLS 2da etapa
   coef_mce <- regresion(x_lag,y_dif)
   #error estimado 2da etapa
   neta_est <- y_dif - (x_lag %*% coef_mce)
   #la garzoneta
   neta_dif <- as.matrix(op_diferencias(neta_est,1))
   neta_dif <- as.matrix(neta_dif[-1,1])
   #neta_lag
   neta_lag <- as.matrix(op_back(neta_est,1))
   neta_lag <- as.matrix(neta_lag[-1,1])
   #OLS neta 3ra etapa
   beta_3e <- regresion(neta_lag,neta_dif)
   #beta ajustado
   coef_adj <- coef - beta_3e
   rta <- print(paste("El coeficiente ajustado de la tercera etapa es:", coef_adj))
}

est_trietapica <- MCE(y,X,p_optimo)

