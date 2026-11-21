library(maxLik)
library(countDM)
library(LambertW)
library(knitr)
library(DescTools)

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

#################likelihood bell touchard###########################
loglik <- function(param , y , x) { #sequencia e quantidade de parametros
  beta <- param[2] 
  alpha <- param[1]
  fi <- param[3]
  nu <- param[4]
  n <- length(y)  
  aux = numeric(0)
  for(i in 1:n){
    aux[i] = TP(y[i] , theta = nu) #TP= quantidade baseada nos numeros de bell
  }
  
  
  gmi = numeric(0) #iniciando vetor vazio
  
  gmi[1] = alpha + x[1]*beta #primeiro elemento de gmi
  
  for(t in 2:n){
    if(y[t-1] == 0){
      yest = 0.1  #evitando problemas com log(0)
    }else{
      yest = y[t-1]
    }
    gmi[t] = alpha + x[t] * beta + fi * (log(yest) - x[t-1]* beta) 
  }
  mit = exp(gmi)
  
  
  ll <- nu * (1 - exp(W(mit/nu))) + y*log(W(mit/nu)) + log(aux) 
  vero = sum(ll[2:n])
  return(vero)
} #fim do loglik 


##################### SIMULAÇÃO ##############################
#Número de simulações de monte carlo

M = 10000

#vetores pra armazenar o q foi estimado
estim_alfa = numeric(M)
estim_beta = numeric(M)
estim_fi = numeric(M)
estim_nu = numeric(M)

sd_alfa = numeric(M)
sd_beta = numeric(M)
sd_fi = numeric(M)
sd_nu = numeric(M)

#Gerando amostras bell-Touchard (Modelo BETO-GARMA)

nu = 1.5
alfa = 0.1
beta = 0.2
fi = 0.1
n = 1000




Ymat = matrix(999 , M , n + 96)


j = 0
m =1 


teve_nas = numeric(0)




# Barra de progresso
pb <- txtProgressBar(min = 0, max = M, style = 3)


while(m <= M) {
  
  cat("\n=================================================\n")
  cat("Amostra válida:", m, "de", M,
      "| Tentativa:", j + 1,
      "|", format(Sys.time(), "%H:%M:%S"), "\n")
  
  tempo.total <- Sys.time()
  tempo.geracao <- Sys.time()
  mu = numeric(0) #vetor vazio pra usar no loop
  y = numeric(0)
  x = rep(0,n + 96)
  x[seq(1, n + 96 , 12)] = 1 #saltos de tamanho 12 do x.
  
  
  
  mu[1] = 3 
  y[1] = rbellt.mu(1,mu[1],nu)
  
  
  
  
  for(i in 2:(n + 96)){
    mu[i] = exp(alfa + x[i]*beta + fi*(log(max(y[i-1] , 0.1)) - x[i-1]*beta))
    y[i] = rbellt.mu(1, mu[i] , nu)
  } #fim do for
  
  
  
  y <- y[-(1:96)]
  x <- x[-(1:96)]
  
  tempo.geracao.fim <- Sys.time()
  
  cat("Tempo para gerar a amostra:",
      round(as.numeric(difftime(tempo.geracao.fim,
                                tempo.geracao,
                                units="secs")),2),
      "segundos\n")
  
  
  tempo.maxlik <- Sys.time()
  
  saida <- maxLik(loglik, y=y, x=x , start = c(1,2,0.4,0.1))
  tempo.maxlik.fim <- Sys.time()
  
  cat("Tempo do maxLik:",
      round(as.numeric(difftime(tempo.maxlik.fim,
                                tempo.maxlik,
                                units="secs")),2),
      "segundos\n")
  
  j = j + 1 
  
  teve_nas[j] = miscTools::stdEr(saida)[4]
  
  
  
  if(!is.na(miscTools::stdEr(saida)[4])){
    if(!is.nan(miscTools::stdEr(saida)[4])){
      if(!is.infinite(miscTools::stdEr(saida)[4])){
        if(!is.na(miscTools::stdEr(saida)[3])){
          if(!is.nan(miscTools::stdEr(saida)[3])){
            if(!is.infinite(miscTools::stdEr(saida)[3])){
              if(!is.na(miscTools::stdEr(saida)[2])){
                if(!is.nan(miscTools::stdEr(saida)[2])){
                  if(!is.infinite(miscTools::stdEr(saida)[2])){
                    if(!is.na(miscTools::stdEr(saida)[1])){
                      if(!is.nan(miscTools::stdEr(saida)[1])){
                        if(!is.infinite(miscTools::stdEr(saida)[1])){
                          
                          
                          
                          
                          estim_alfa[m] = saida$estimate[1]
                          estim_beta[m] = saida$estimate[2]
                          estim_fi[m] = saida$estimate[3]
                          estim_nu[m] = saida$estimate[4]
                          
                          sd_alfa[m] = miscTools::stdEr(saida)[1]
                          sd_beta[m] = miscTools::stdEr(saida)[2]
                          sd_fi[m] = miscTools::stdEr(saida)[3]
                          sd_nu[m] = miscTools::stdEr(saida)[4]
                          
                          m = m + 1
                          setTxtProgressBar(pb, m - 1)
                        }}}}}}}}}}}}
  
  
  
  
  miscTools::stdEr(saida)
  cat("Tempo total da tentativa:",
      round(as.numeric(difftime(Sys.time(),
                                tempo.total,
                                units="secs")),2),
      "segundos\n")
}#fim do while

close(pb)
########fim monte /carlo 

######################### Histogramas ###########################
hist(estim_nu , prob = TRUE)
hist(Winsorize(estim_nu, val= quantile(estim_nu ,probs = c(0.05, 0.95))))
hist(estim_alfa , prob = TRUE)
hist(estim_fi , prob = TRUE)
hist(estim_beta , prob = TRUE)

write.table(estim_nu,file="estimativas_nu")
write.table(estim_alfa,file="estimativas_alfa")
write.table(estim_beta,file="estimativas_beta")
write.table(estim_fi,file="estimativas_fi")

mean(Winsorize(estim_nu, val= quantile(estim_nu ,probs = c(0.05, 0.95))))
mean(Winsorize(estim_alfa, val = quantile(estim_alfa , probs = c(0.05, 0.95))))
mean(Winsorize(estim_beta,val = quantile( estim_beta , probs = c(0.05, 0.95))))
mean(Winsorize(estim_fi, val= quantile(estim_fi, probs = c(0.05, 0.95))))

estim_nu_winsor = Winsorize(estim_nu, val= quantile(estim_nu ,probs = c(0.05, 0.95)))
estim_alfa_winsor = Winsorize(estim_alfa, val = quantile(estim_alfa , probs = c(0.05, 0.95)))
estim_beta_winsor = Winsorize(estim_beta,val = quantile( estim_beta , probs = c(0.05, 0.95)))
estim_fi_winsor   = Winsorize(estim_fi, val= quantile(estim_fi, probs = c(0.05, 0.95)))


# Média
EVM_alfa = mean(Winsorize(estim_alfa, val = quantile(estim_alfa , probs = c(0.05, 0.95))))
EVM_beta = mean(Winsorize(estim_beta,val = quantile( estim_beta , probs = c(0.05, 0.95))))
EVM_fi = mean(Winsorize(estim_fi, val= quantile(estim_fi, probs = c(0.05, 0.95))))
EVM_nu = mean(Winsorize(estim_nu, val= quantile(estim_nu ,probs = c(0.05, 0.95))))



# Cálculo do EQM
EQM_alfa = mean((estim_alfa_winsor - alfa)^2)
EQM_beta = mean((estim_beta_winsor - beta)^2)
EQM_fi = mean((estim_fi_winsor - fi)^2)
EQM_nu = mean((estim_nu_winsor - nu)^2)

#DESVIO PADRAO
SD_alfa = mean(sd_alfa)
SD_beta = mean(sd_beta)
SD_fi = mean(sd_fi)
SD_nu = mean(sd_nu)

#vies
vies_alfa = mean(estim_alfa - alfa)
vies_beta = mean(estim_beta - beta)
vies_fi = mean(estim_fi - fi)
vies_nu = mean(estim_nu - nu)


round(SD_alfa, 4)
round(SD_beta, 4)
round(SD_fi, 4)
round(SD_nu ,4)

round(vies_alfa ,4)
round(vies_beta ,4)
round(vies_fi ,4)
round(vies_nu ,4)

round(EQM_alfa, 4)
round(EQM_beta,4)
round(EQM_fi, 4)
round(EQM_nu,4)

round(EVM_alfa, 4)
round(EVM_beta,4)
round(EVM_fi, 4)
round(EVM_nu,4)




tabela_resultados <- data.frame(
  Parâmetro = c("alfa", "beta", "fi", "nu"),
  Media = round(c(EVM_alfa, EVM_beta, EVM_fi, EVM_nu), 4),
  EQM = round(c(EQM_alfa, EQM_beta, EQM_fi, EQM_nu), 4),
  Desvio_Padrao = round(c(SD_alfa, SD_beta, SD_fi, SD_nu), 4),
  Vies = round(c(vies_alfa, vies_beta, vies_fi, vies_nu), 4)
)


print(tabela_resultados)




