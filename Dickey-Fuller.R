#### Dickey Fuller Aumentado ----
# Seleccion P óptimo E3
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

# Paso 1 - E3 / ¿gamma = 0?
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
    
    #error
    e_est <- y - (x %*% coef)
    
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
    
    #quitar NAs iniciales
    x <- as.matrix(m_completa[-(1:p_opt),])
    y <- as.matrix(delta_yt(yt)[-(1:p_opt),])
    
    #regresion
    coef <- regresion(x,y)
    gamma <- coef[3,1]
    
    #vcv
    #parametros iniciales
    T <- nrow(x)
    k <- ncol(x)
    
    #error
    e_est <- y - (x %*% coef)
    
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

# Si rechazo -> fin
# Si no rechazo -> Reestimar nuevamente el P óptimo
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

# Paso 2 - E3 / ¿gamma y beta = 0?
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

# Si rechazo -> Paso 3 - E3 ¿Beta = 0?
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
  # Si rechazo -> fin
# Si no rechazo Paso 2 o 3 -> Reestimar nuevamente el P óptimo

P_opt_E2 <- function(yt) {
  #Pasar datos a Yt y revisar length/nrow
  yt <- as.matrix(yt)
  #Y diferenciado
  y <- as.matrix(delta_yt(yt))
  
  #xdiseño
  contador <- 1:length(yt)
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

# Paso 1 - E2 / ¿gamma = 0?
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

# Paso 2 - E2 / ¿alpha y gamma = 0?
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