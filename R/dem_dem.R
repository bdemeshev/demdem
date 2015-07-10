


## ----, message=FALSE, warning=FALSE--------------------------------------
library("ggplot2")
library("knitr")
library("MCMCpack")
library("pander")
library("dplyr")
library("MHadaptive")
library("mvtnorm")

library("microbenchmark") # test speed of various approaches


## ----, "define MH functions"---------------------------------------------
li_rho <- function(rho) {
  if ((rho > -1) & (rho<1)) {
    e <- y_star - rho*Wy_star - Z %*% phi # y?
    eVe <- sum(e^2/v) # faster than eVe <- t(e) %*% invV %*% e    
    ans <- log(det(diag(n)-rho*W)) - 0.5*eVe/s2
  } else ans <- -Inf
  return(ans)
}

li_q <- function(q) {
  ans <- ifelse(q>0, 0.5*n*q*log(0.5*q) - n*lgamma(0.5*q) - kappa*q + (a_q-1)*log(q), -Inf)
  return(ans)
}


# загружаем данные
## ------------------------------------------------------------------------
h <- read.csv("../data/regional_data.csv")
W <- read.csv("../data/Wb.csv", header=FALSE)
W <- as.matrix(W)



# Переименуем для удобства:
## ------------------------------------------------------------------------
h <- dplyr::rename(h, y_star = Y, Wy_star = WbY, y_0 = ln.gdppercappp., region = X) %>% 
  dplyr::select(-number)
glimpse(h)

# Удаляем Калининградскую Область! 
# Она ни с кем не граничит, поэтому при пограничной $W$ 
# нарушается свойство $W\vec{1}=\vec{1}$ и определитель $det(I_{n\times n }-\rho W)$ оказывается отрицательным. 
# А он фигурирует в плотностях. 
## ------------------------------------------------------------------------
n_kalin <- which(h$region=="Kaliningrad region")
h <- filter(h, !region=="Kaliningrad region")
W <- W[-n_kalin,-n_kalin]

# Априорные распределения:
#   
#   1. $\rho \sim U[-1;1]$
#   
#   2. $\phi \sim $ diffuse
# 
# 3. $\sigma^2_{\varepsilon} \sim $ standard diffuse ???
# 
# 4. $q \sim \Gamma(a_q, b_q)$
#   
#   5. $v_i^{-1} | q \sim iid \chi^2(q)$, $v_i$ --- diagonal of $V$
#   
#   6. $Var(\varepsilon) = \sigma^2_{\varepsilon} V$ ?
# 
# Упрощенная модель 
# \[
#   y^* = \rho Wy^* + \alpha i + \beta y_0 + X\gamma + \varepsilon
#   \]
# 
# Полная из статьи
# \[
#   y^* = \rho Wy^* + \alpha i + \beta y_0 +\theta Wy_0 + X\gamma + WX\xi +  \varepsilon
#   \]
# 
# 
# 
# Упрощения: $\theta=0$, $\xi=0$
## ------------------------------------------------------------------------
X <- as.matrix(h[,5:17])
WX <- W %*% X
C <- WX # japan notation
Wy_0 <- W %*% h$y_0
y_0 <- h$y_0
y_star <- h$y_star
Wy_star <- h$Wy_star

n <- nrow(h)
m <- ncol(X)
l <- ncol(C)

Z <- cbind( rep(1, n), y_0, Wy_0, X, C) 


# in our case l=m, as C=WX


## ----, eval=FALSE--------------------------------------------------------
## phi <- c(alpha, beta, theta, gamma, xi)
## pars <- c(rho, phi, s2, v, q) # last change: q


# \[
#   (\rho, \alpha, \beta, \theta, \gamma, \xi, \sigma^2, v, q)
#   \]
# 
# 
# $\xi, \gamma \in R^{`r m`}$, $v\in R^{`r n`}$
#   
#   
#   Параметры априорных распределений и кое-какие предрасчеты:

## ----, results='asis'----------------------------------------------------
a_a <- 0.001 # page 63 bottom  or a_sigma top of the same page:)
b_a <- 0.001
a_q <- 1
b_q <- 0.5

r <- rep(0, m+l+3)
S <- 10^12 * diag(m+l+3)

# precalculate
invS <- solve(S)
invSr <- invS %*% r

# Инициализируем параметры случайно по априорному распределению
# 
# Если $q_{init}=0$, то все `rchisq` будут равны нулю. Так нам не надо!


## ------------------------------------------------------------------------
pars_init <- rep(0, 6+m+l+n)


q <- rgamma(a_q,b_q)

rho <- runif(1)
phi <- rmvnorm(1, r, S)
s2 <- rinvgamma(1, a_a/2, b_a/2)
v <- 1/rchisq(n, df=q)

pars_init[1] <- rho
pars_init[2:(m+l+4)] <- phi
pars_init[5+m+l] <- s2
pars_init[(6+m+l): (5+m+l+n)] <- v
pars_init[6+m+l+n] <- q

# model_0 <- lm(h$y_star ~ h$Wy_star + h$y_0 + Wy_0 + X + WX)
# coefs <- coef(model_0)
# pander(model_0)

# y^* = \rho Wy^* + \alpha i + \beta y_0 +\theta Wy_0 + X\gamma + WX\xi +  \varepsilon
# (\rho, \alpha, \beta, \theta, \gamma, \xi, \sigma^2, v)
# rho_init <- coefs[2]
# alpha_init <- coefs[1]
# v_init <- rep(1, n)
# btgxi <- coefs[-2:-1]
#s2_init <- deviance(model_0)/df.residual(model_0)
#pars_init <- c(rho_init, alpha_init, btgxi, s2_init, v_init)

# Именуем вектор параметров:

## ------------------------------------------------------------------------
names(pars_init)[1:4] <- c("rho", "alpha", "beta", "theta")
names(pars_init)[5:(4+m)] <- colnames(X)
names(pars_init)[(5+m):(4+m+l)] <- paste0(rep("w_",l),colnames(X))
names(pars_init)[5+m+l] <- "s2"
names(pars_init)[(6+m+l):(5+m+l+n)] <- paste0(rep("v",n),1:n)
names(pars_init)[6+m+l+n] <-"q"

# MCMC. Уже сохраненные результаты хранятся в `/estimation/`+pars_file.

## ------------------------------------------------------------------------
n_sim <- 20000 # jap: 15000
# n_burnin <- 5000 # jap: 5000 # not used, we remove them later
n_mh_iters <- 20000 # jap?
n_chains <- 9

set.seed(13)  # wish your good luck, MCMC
for (chain_no in 1:n_chains) { # we create 9 chains :)

pars <- matrix(0, nrow=n_sim, ncol=length(pars_init))
colnames(pars) <- names(pars_init)

pars_file <- paste0("pars_chain_",chain_no,".Rds")
pars_full_path <- paste0("./estimation/",pars_file)

# if files already exist we add simultations to them!
if (pars_file %in% list.files("./estimation/")) {
  old_pars <- readRDS(pars_full_path)
  n_sim_done <- sum(old_pars[,"q"]>0)
  pars[1:n_sim_done,] <- old_pars[1:n_sim_done,]
  j_start <- n_sim_done + 1
} else {  
  pars[1, ] <- pars_init
  j_start <- 2
}





# go-go-go (one chain)

## ------------------------------------------------------------------------
time_start <- proc.time()["elapsed"]

# (\rho, \alpha, \beta, \theta, \gamma, \xi, \sigma^2, v)
for (j in j_start:n_sim) {
  time_now <- proc.time()["elapsed"]
  time_per_iter <- round((time_now-time_start)/(j-j_start-1), 2)
  message("Doing mcmc step ",j,"/",n_sim,". Time per iter=",time_per_iter)
 
  rho <- pars[j-1, 1]
  s2 <- pars[j-1, 5+m+l]
  v <- pars[j-1, (6+m+l):(5+l+m+n)]
  q <- pars[j-1,6+m+l+n]
  # phi <- pars[j-1, 2:(m+l+4) ] # one may drop --- generated in step a
   
  
  
  # step a
  invV <- diag(1/v) # solve(diag(v)) is much much slower
  term_a <- t(Z) %*% invV/s2
  
  y_tilde <- y_star-rho*Wy_star # y? page 63 eq 25 in Seya 2012
  S_star <- solve(term_a %*% Z + invS)
  r_star <- S_star %*% (term_a %*% y_tilde + invSr)
  
  phi <- as.vector(rmvnorm(1, r_star, S_star))
  # as.vector is needed because rmvnorm returns matrix of size 1x(m+l+3)
  # and we use phi as vector later in multiplication
  
  # step b
  e <- as.vector(y_star - rho*Wy_star - Z %*% phi) # y? page 63 eq 25 in Seya 2012
  eVe <- sum(e^2/v) # 7 times faster than eVe <- t(e) %*% invV %*% e
  s2 <- rinvgamma(1, (n+a_a)/2, (eVe+b_a)/2) # check parameter specification!
  
  
  # step c
  rw_c <- rchisq(n, df=q+1) 
  v <- (q+e^2/s2)/rw_c # error in Seya 2012 (27), see Kakamura 2008 (3) 
  
  # step d
  mcmc_out <- Metro_Hastings(li_rho, pars = rho, iterations = n_mh_iters, quiet = TRUE)
  rho <- tail(mcmc_out$trace,1)

  # step e
  kappa <- 0.5 * (sum(log(v))+sum(1/v)) + b_q
  mcmc_out <- Metro_Hastings(li_q, pars = q, iterations = n_mh_iters, quiet = TRUE)
  q <- tail(mcmc_out$trace,1)
  
  
  # save
  rho -> pars[j, 1]
  s2 -> pars[j, 5+m+l]
  v -> pars[j, (6+m+l):(5+m+l+n)]
  q -> pars[j, 6+m+l+n]
  phi -> pars[j, 2:(m+l+4) ]
  
  if (j %% 10 == 0) {
    saveRDS(pars, pars_full_path)
    message("Chain progress saved in '", pars_file,"'")
  }
}
saveRDS(pars, pars_full_path)

} # end chain no cycle


## ----, "MH rho test", eval=FALSE-----------------------------------------
## mcmc_out <- Metro_Hastings(li_rho, pars = 0, iterations = n_mh_iters)
## rho <- tail(mcmc_out$trace,1)
## rho


## ----, "MH q test", eval=FALSE-------------------------------------------
## 
## kappa <- 0.5 * (sum(log(v))+sum(1/v)) + b_q
## mcmc_out <- Metro_Hastings(li_q, pars = 1, iterations = n_mh_iters)
## q <- tail(mcmc_out$trace,1)
## q


