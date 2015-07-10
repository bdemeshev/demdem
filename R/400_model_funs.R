library("readr")
library("MHadaptive")
library("MCMCpack")
library("mvtnorm")
library("dplyr")


create_model_list <- function() {
  df <- expand.grid(type="mcmc", chain_no=1:5, status="not estimated", 
                    sample=c("west","east","all"), W=c("border","road"))
  df <- df %>% mutate_each("as.factor",type, status, sample, W) %>% mutate(seed=chain_no+42)
  df <- df %>% mutate(file=paste0(type,"_",W,"_",sample,"_",chain_no,".Rds"), id=row_number())
  df <- reshape2::melt(df, id.vars="id") %>% arrange(id)
  
  return(df)
}


estimate_models <- function(mlist, parallel = parallel, no_reestimation=TRUE) {
  if (no_reestimation) mlist_todo <- dplyr::filter(mlist, !status=="estimated")
  
  model_ids <- unique(mlist_todo$id)
  
  
  if (parallel=="windows") {
    library("doSNOW")
    cl <- makeCluster(10, outfile="") # number of CPU cores
    registerDoSNOW(cl)
    
    foreach(i=model_ids, .packages=c("mvtnorm","MHadaptive","MCMCpack")) %dopar% {
      model_info <- mlist_todo %>% dplyr::filter(id==i)
      status <- estimate_model(model_info)
      mlist$status[mlist$id==i] <- status
    }
    stopCluster(cl)
  }
  
  if (parallel=="off") {
    foreach(i=model_ids, .packages=c("mvtnorm","MHadaptive","MCMCpack")) %do% {
      model_info <- mlist_todo %>% dplyr::filter(id==i)
      status <- estimate_model(model_info)
      mlist$status[mlist$id==i] <- status
    }
  }
  return(mlist) # statuses are updated
}



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



estimate_model <- function(model_info) {
  minfo <- reshape2::dcast(model_info, id~variable)
  model_full_path <- paste0("../estimation/models/",minfo$file)
  set.seed(minfo$seed)  # wish your good luck, MCMC
  
  
  
  if (minfo$W=="border") {
    W <- read_csv("../data/W_border.csv")
    W <- as.matrix(W)
    h <- read_csv("../data/h_nk.csv")
  }

  if (minfo$W=="road") {
    W <- read_csv("../data/W_road.csv")
    W <- as.matrix(W)
    h <- read_csv("../data/h.csv")
  }
  
  if (minfo$sample=="west") {
    use_rows <- h$location=="west"
    W <- W[use_rows,use_rows] # get requested part of W
    h <- h[use_rows,] # get requested part of h
  }
  
  if (minfo$sample=="east") {
    use_rows <- h$location=="east"
    W <- W[use_rows,use_rows] # get requested part of W
    h <- h[use_rows,] # get requested part of h
  }
  
  
  
  X <- as.matrix(h[,5:17])
  WX <- W %*% X
  C <- WX # japan notation
  y_0 <- h$y_0
  Wy_0 <- W %*% y_0
  y_star <- h$y_star
  Wy_star <- W %*% y_star
  
  n <- nrow(h)
  m <- ncol(X)
  l <- ncol(C)
  
  Z <- cbind( rep(1, n), y_0, Wy_0, X, C) 
  
  
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
  n_sim <- 15000 # jap: 15000
  # n_burnin <- 5000 # jap: 5000 # not used, we remove them later
  n_mh_iters <- 10000 # jap?
  n_chains <- 9
  

    pars <- matrix(0, nrow=n_sim, ncol=length(pars_init))
    colnames(pars) <- names(pars_init)
    
  
    
    # if files already exist we add simulations to them!
    if (minfo$file %in% list.files("../estimation/models/")) {
      old_pars <- readRDS(model_full_path)
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
    if (j_start <= n_sim) for (j in j_start:n_sim) {
       
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
      
      
      time_now <- proc.time()["elapsed"]
      time_per_iter <- round((time_now-time_start)/(j-j_start+1), 2)
      message("Mcmc step ",j,"/",n_sim," done. Mean time per iter=",time_per_iter)
      
      
      if (j %% 10 == 0) {
        saveRDS(pars, model_full_path)
        message("Chain progress saved in '", minfo$file,"'")
      }
    }
    saveRDS(pars, model_full_path)
    

  status <- "estimated"
  
  return(status)
}



