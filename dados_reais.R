#DADOS REAIS

library(maxLik)
library(countDM)
library(LambertW)


dbellt.mu <- function (x, mu, nu, log = FALSE) 
{
  if (any(mu <= 0)) 
    stop(paste("mu must be greater than 0"))
  if (any(nu <= 0)) 
    stop(paste("nu must be greater than 0"))
  TPx <- c()
  for (i in 1:length(x)){
    TPx[i] <- TP(x[i], theta = nu)
  }
  lf <- x * log_W(mu/nu) + nu * (-exp(W(mu/nu)) + 1) + log(TPx) - 
    lgamma(x+1)
  if (log == TRUE) {
    return(lf)
  }
  else {
    return(exp(lf))
  }
}

############# Função de distribuição acumulada (mu,nu) ###################

pbellt.mu <- function (q, mu, nu, lower.tail = TRUE, log.p = FALSE) 
{
  if (any(mu <= 0)) 
    stop(paste("mu must be greater than 0"))
  if (any(nu <= 0)) 
    stop(paste("nu must be greater than 0"))
  if (any(q < 0)) 
    stop(paste("q must be 0 or greater than 0"))
  prob <- c()
  for (i in 1:length(q)) prob[i] <- sum(dbellt.mu(0:q[i], mu = mu, 
                                                  nu = nu, log = FALSE))
  if (lower.tail == FALSE) {
    if (log.p == FALSE) {
      prob <- 1 - prob
    }
    else prob <- log(1 - prob)
  }
  else if (log.p == TRUE) 
    prob <- log(prob)
  prob
}

################## Função quantilica (mu,nu) ###########################

qbellt.mu <- function (p, mu, nu, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0)) 
    stop(paste("mu must be greater than 0"))
  if (any(nu <= 0)) 
    stop(paste("nu must be greater than 0"))
  if (any(p < 0) | any(p > 1)) 
    warning("pi must be in 0, 1")
  n <- length(p)
  k <- rep(0, n)
  for (i in 1:n) {
    q <- 0
    while (q <= p[i]) {
      q <- pbellt.mu(k[i], mu, nu)
      k[i] <- k[i] + 1
    }
  }
  if (log.p) {
    k <- log(k - 1)
  }
  else k <- k - 1
  k
}

########## Geração de números (mu,nu) ##########################

rbellt.mu <- function (n, mu, nu) {
  if (any(mu <= 0)) 
    stop(paste("mu must greated than 0"))
  if (any(nu <= 0)) 
    stop(paste("nu must greated than 0"))
  if (any(n <= 0)) 
    stop(paste("n must be a positive integer"))
  n <- ceiling(n)
  p <- runif(n)
  z <- qbellt.mu(p, mu = mu, nu = nu, lower.tail = TRUE, 
                 log.p = FALSE)
  as.integer(z)
}

dados <- sih_cnv_mirn192114187_19_254_71
dados$tempo <- 1:nrow(dados)
head(dados)
str(dados)


mes <- iconv(as.character(dados$mes),
             from = "",
             to = "ASCII//TRANSLIT")

ordem_meses <- c("Janeiro", "Fevereiro", "Marco",
                 "Abril", "Maio", "Junho",
                 "Julho", "Agosto", "Setembro",
                 "Outubro", "Novembro", "Dezembro")

dados$mes_ordem <- match(mes, ordem_meses)


dados <- dados[order(dados$Ano, dados$mes_ordem), ]



dados$periodo_chuvoso <- ifelse(
  mes[order(dados$Ano, dados$mes_ordem)] %in%
    c("Marco", "Abril", "Maio", "Junho"),
  1, 0
)


table(dados$periodo_chuvoso)


y <- dados$internacoes
x <- dados$periodo_chuvoso

n <- length(y)

summary(y)
mean(y)
var(y)
mean(y == 0)


plot(dados$tempo, dados$internacoes,
     type = "o",
     pch = 16,
     xlab = "Tempo (meses)",
     ylab = "Número de internações",
     main = "Série temporal de internações por dengue")

######################## FUNCAO DE VEROSSIMILHANCA ########################

loglik <- function(param, y, x) {
  
  alpha <- param[1]
  beta  <- param[2]
  fi    <- param[3]
  nu    <- param[4]
  
  n <- length(y)
  
  #Garantir nu positivo
  if (nu <= 0 || !is.finite(nu)) {
    return(-1e100)
  }
  
  #Preditor linear
  eta <- numeric(n)
  
  #Primeira observacao
  eta[1] <- alpha + beta*x[1]
  
  #Demais observacoes
  for (t in 2:n) {
    
    if (y[t-1] == 0) {
      yest <- 0.1
    } else {
      yest <- y[t-1]
    }
    
    eta[t] <- alpha +
      beta*x[t] +
      fi*(log(yest) - beta*x[t-1])
  }
  
  #Media condicional
  mu <- exp(eta)
  
  if (any(!is.finite(mu)) || any(mu <= 0)) {
    return(-1e100)
  }
  
  #Parametro auxiliar
  a <- W(mu/nu)
  
  #Quantidades de Touchard
  TPy <- numeric(n)
  
  for (t in 1:n) {
    TPy[t] <- TP(y[t], theta = nu)
  }
  
  if (any(!is.finite(TPy)) || any(TPy <= 0)) {
    return(-1e100)
  }
  
  #Log-verossimilhanca
  ll <- nu*(1 - exp(a)) +
    y*log(a) +
    log(TPy) -
    lgamma(y + 1)
  
  #A primeira observacao nao entra na verossimilhanca
  return(sum(ll[2:n]))
}

######################## ESTIMACAO DOS PARAMETROS ########################

valores_iniciais <- c(0.5, 1, 0.4, 0.7)

#Estimacao por maxima verossimilhanca
saida_real <- maxLik(
  loglik,
  y = y,
  x = x,
  start = valores_iniciais
)

#Resultado da estimacao
summary(saida_real)

#Verificar convergencia
saida_real$code

# Estimativas
estimativas <- saida_real$estimate

alpha_hat <- estimativas[1]
beta_hat  <- estimativas[2]
fi_hat    <- estimativas[3]
nu_hat    <- estimativas[4]

# Calcular as médias condicionais estimadas
eta_hat <- numeric(n)

eta_hat[1] <- alpha_hat + beta_hat*x[1]

for (t in 2:n) {
  
  yest <- if (y[t-1] == 0) {
    0.1
  } else {
    y[t-1]
  }
  
  eta_hat[t] <- alpha_hat +
    beta_hat*x[t] +
    fi_hat*(log(yest) - beta_hat*x[t-1])
}

mu_hat <- exp(eta_hat)

# Erros padrão
erros_padrao <- sqrt(diag(vcov(saida_real)))

# Limites do intervalo aproximado
limite_inferior <- estimativas - 2 * erros_padrao
limite_superior <- estimativas + 2 * erros_padrao

# Tabela de resultados
tabela_parametros <- data.frame(
  Parametro = c("alpha", "beta", "fi", "nu"),
  Estimativa = estimativas,
  Erro_Padrao = erros_padrao,
  Limite_Inferior = limite_inferior,
  Limite_Superior = limite_superior
)

# Arredondar resultados
tabela_parametros[, 2:5] <- round(tabela_parametros[, 2:5], 4)

print(tabela_parametros)


######################## INTERPRETACAO DA COVARIAVEL ########################

#Efeito multiplicativo do periodo chuvoso
exp(beta_hat)

#> A média condicional estimada de internações no período chuvoso é aproximadamente
#> 2,54 vezes a média correspondente ao período não chuvoso, 
#> mantendo constante a contribuição dinâmica anterior

#AIC
AIC_modelo <- -2*as.numeric(logLik(saida_real)) +
  2*length(estimativas)

BIC_modelo <- -2*as.numeric(logLik(saida_real)) +
  log(n)*length(estimativas)


  AIC_modelo
  BIC_modelo


######################## RESIDUOS QUANTILICOS ALEATORIZADOS ########################

set.seed(123)

res_quant <- rep(NA, n)

for (t in 2:n) {
  
  # Limite inferior da probabilidade acumulada
  if (y[t] == 0) {
    
    inferior <- 0
    
  } else {
    
    inferior <- sum(
      dbellt.mu(
        0:(y[t]-1),
        mu = mu_hat[t],
        nu = nu_hat
      )
    )
  }
  
  # Limite superior da probabilidade acumulada
  superior <- sum(
    dbellt.mu(
      0:y[t],
      mu = mu_hat[t],
      nu = nu_hat
    )
  )
  
  # Proteção contra erros numéricos
  inferior <- max(0, min(1, inferior))
  superior <- max(inferior, min(1, superior))
  
  # Sorteio uniforme
  u <- runif(1, inferior, superior)
  
  # Evitar valores iguais a zero ou um
  u <- max(.Machine$double.eps,
           min(1 - .Machine$double.eps, u))
  
  # Resíduo quantílico
  res_quant[t] <- qnorm(u)
}
  
  length(res_quant)
  
  sum(!is.na(res_quant))
  
  head(res_quant)
  
  tail(res_quant)
  
  which(!is.na(res_quant))

  ######################## 10. DIAGNOSTICO DOS RESIDUOS QUANTILICOS ########################
  
  
  # Serie temporal
  g1 = plot(res_quant,
       type = "b",
       pch = 16,
       xlab = "Tempo",
       ylab = "Resíduos Quantílicos",
       main = "Resíduos Quantílicos Aleatorizados")
  
  abline(h = 0, col = "red", lty = 2)
  
  # Histograma
 g2= hist(res_quant,
       probability = TRUE,
       main = "Histograma dos Resíduos Quantílicos",
       xlab = "Resíduos Quantílicos")
  
  curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
  
  # Grafico Q-Q
  qqnorm(res_quant,
         main = "Gráfico Q-Q dos Resíduos Quantílicos")
  
  qqline(res_quant, col = "red", lwd = 2)
  
  par(mfrow = c(1, 2))
  # ACF
  acf(res_quant[-1],
      main = "ACF dos Resíduos Quantílicos")
  
  # PACF
  pacf(res_quant[-1],
       main = "PACF dos Resíduos Quantílicos")
  
  
  ######################## 11. TESTE DE LJUNG-BOX ########################
  
  Box.test(
    res_quant[-1],
    lag = 12,
    type = "Ljung-Box"
  )
  
  
  ######################## 12. RESULTADOS FINAIS ########################
  
  cat("\nESTIMATIVAS DOS PARAMETROS\n")
  print(tabela_parametros)
  
  cat("\nEFEITO MULTIPLICATIVO DA COVARIAVEL\n")
  print(exp(beta_hat))
  
  cat("\nLOG-VEROSSIMILHANCA\n")
  print(logLik(saida_real))
  
  cat("\nAIC\n")
  print(AIC_modelo)
  
  cat("\nTESTE DE LJUNG-BOX\n")
  print(Box.test(res_quant[-1],
                 lag = 12,
                 type = "Ljung-Box"))
  
  
  
  #garmafit c/ poisson e binominal negativa-------------
 library(gamlss)
  
  modelo_nb <- gamlss.util::garmaFit(
    formula = internacoes ~ periodo_chuvoso,
    order = c(1, 0),
    data = dados,
    family = gamlss.dist::NBI(),
    alpha = 0.1
  )
  
  summary(modelo_nb)
  
  AIC(modelo_nb)
  
  
  
  
  trace(
    "garmaFit",
    where = asNamespace("gamlss.util"),
    tracer = quote(BItrue <- FALSE),
    print = FALSE
  )

  modelo_pois <- gamlss.util::garmaFit(
    formula = internacoes ~ periodo_chuvoso,
    order = c(1, 0),
    data = dados,
    family = gamlss.dist::PO(),
    alpha = 0.1
  )  

  summary(modelo_pois)
  
  AIC(modelo_pois)
  
  modelo_pois$G.deviance
  
  untrace(
    "garmaFit",
    where = asNamespace("gamlss.util")
  )
  