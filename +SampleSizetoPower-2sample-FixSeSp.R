SampleSizetoPower2sampleFixSeSp = function(
                              true_p1,
                              true_p2,
                              Se1,
                              Sp1,
                              Se2,
                              Sp2,
                              alpha = 0.05,
                              pwr = 0.8,
                              min_n = 20,
                              max_n = 5000,
                              drprate = 0,
                              R = 20000,
                              met = "Wi"){
  
  # Function to calculate the required sample size for a two-sample binomial test when Se1, Se2, Sp1 and Sp2 are considered known, fixed.
  #  
  # Args:
  #   true_p1: Real parameter value in group 1.
  #   true_p2: Real parameter value in group 2.
  #   Se1: Sensitivity of the diagnostic test used for observing true_p1 in group 1.
  #   Sp1: Specificity of the diagnostic test used for observing true_p1 in group 1.
  #   Se2: Sensitivity of the diagnostic test used for observing true_p2 in group 2.  
  #   Sp2: Specificity of the diagnostic test used for observing true_p2 in group 2.
  #   alpha: Type I error probability. Specify as a number between 0 and 1 rather than a percent value.
  #   pwr: Prescribed power of the test. Specify as a number between 0 and 1 rather than a percent value. 
  #   min_n: Minimum sample size (smaller sample sizes will not be examined).
  #   max_n: Maximum sample size (if power is below pwr for all n less than or equal to max_n, function returns 0).
  #   drp_rate: Highest drop-out rate for which power must not fall below pwr. Specify as a number between 0 and 1 rather than a percent value. 
  #   R: Number of simulations to be used.
  #   met: Method - Available methods: "Wald", "Wilson", "Agresti-Coull", "Clopper-Pearson", "Blaker"
  #
  # Returns:
  #   The minimal sample size n satisfying that for all sample sizes from drp_rate*n to n power is at least pwr.
  
  # Nested functions:
  
  CI_Reiczigel2010_Zou_RD_FixSeSp <- function(
    r1_nprev,
    r1_kprev,
    Se1_const,
    Sp1_const,
    r2_nprev,
    r2_kprev,
    Se2_const,
    Sp2_const,
    conflevel = 0.95,
    met= "Wi",
    alt = "two.sided")
  {
    
    # Function to calculate the Confidence Interval (CI) for the Risk Difference (RD) between two groups with the method of Hársfalvi and Singer (2023).
    #  
    # Args:
    #   r1_nprev: The number of individuals that are tested in group 1.
    #   r1_kprev: Number of individuals with positive test results in group 1.
    #   Se1_const: Sensitivity of the diagnostic test used for testing in group 1.
    #   Sp1_const: Specificity of the diagnostic test used for testing in group 1.
    #   r2_nprev: The number of individuals that are tested in group 1.
    #   r2_kprev: Number of individuals with positive test results in group 1.
    #   Se2_const: Sensitivity of the diagnostic test used for testing in group 2.
    #   Sp2_const: Specificity of the diagnostic test used for testing in group 2.
    #   conf.level: Confidence level.
    #   met: Method - Available methods: "Wald", "Wilson", "Agresti-Coull", "Clopper-Pearson", "Blaker".
    #   alt: alternative - "two.sided", "less", "greater".
    #
    # Returns:
    #   A vector containing the LCL and UCL of the CI for the RD.
    
    # Nested function:
    
    BinomialConfInt <- function(n,
                                y,
                                conf.level = 0.95,
                                alt = "two.sided",
                                met = "Clopper-Pearson") {
      
      # Function to calculate a binomial confidence interval with several different methods
      #  
      # Args:
      #   n: The number of individuals that are tested.
      #   y: Number of individuals with positive test results.
      #   conf.level: Confidence level.
      #   alt: alternative - "two.sided", "less", "greater".
      #   met: Method - Available methods: "Wald", "Wilson", "Agresti-Coull", "Clopper-Pearson", "Blaker".
      #
      # Returns:
      #   A vector containing the LCL and UCL of the binomial CI.
      
      p_obs <- y / n
      alpha <- 1 - conf.level
      z1s <- qnorm(1 - alpha)
      z2s <- qnorm(1 - alpha / 2)
      cilower <- 0
      ciupper <- 1
      
      # Methods
      if (met == "Wald") {
        if (alt == "two.sided") {
          CI <- c(max(p_obs - z2s * sqrt(p_obs * (1 - p_obs)/(n)), 0), min(p_obs + z2s * sqrt(p_obs * (1 - p_obs)/(n)), 1))
        } else {
          if (alt == "less") {
            CI <- c(0, min(p_obs + z1s * sqrt(p_obs * (1 - p_obs)/(n)), 1))
          } else {
            if (alt == "greater") {
              CI = c(max(p_obs - z1s * sqrt(p_obs * (1 - p_obs)/(n)), 0), 1)
            } else {
              stop("Argument alt must be specified!")
            }
          }
        }
      }
      if (met == "Wilson") {
        if (alt == "two.sided") {
          CI <- c((2 * n * p_obs + z2s^2 - z2s * sqrt(z2s^2 + 4 * n * p_obs * (1 - p_obs)))/(2 * (n + z2s^2)),
                  (2 * n * p_obs + z2s^2 + z2s * sqrt(z2s^2 + 4 * n * p_obs * (1-p_obs)))/(2 * (n + z2s^2)))
        } else {
          if (alt == "less") {
            CI <- c(0, (2 * n * p_obs + z1s^2 + z1s * sqrt(z1s^2 + 4 * n * p_obs * (1 - p_obs)))/(2 * (n + z1s^2)))
          } else {
            if (alt == "greater") {
              CI <- c((2 * n * p_obs + z1s^2 - z1s * sqrt(z1s^2 + 4 * n * p_obs * (1-p_obs)))/(2 * (n + z1s^2)), 1)
            } else {
              stop("Argument alt must be specified!")
            }
          }
        }
      }
      if (met == "Agresti-Coull") {
        phat1s = (p_obs * n + (z1s^2) / 2) / (n + z1s^2)
        phat2s = (p_obs * n + (z2s^2) / 2) / (n + z2s^2)
        nhat1s = n + z1s^2
        nhat2s = n + z2s^2
        if (alt == "two.sided") {
          CI <- c(max(phat2s - z2s * sqrt(phat2s * (1 - phat2s)/(nhat2s)), 0),
                  min(phat2s + z2s * sqrt(phat2s * (1 - phat2s)/(nhat2s)), 1))
        } else {
          if (alt == "less") {
            CI <- c(0, min(phat1s + z1s * sqrt(phat1s * (1 - phat1s)/(nhat1s)), 1))
          } else {
            if (alt == "greater") {
              CI <- c(max(0, phat1s - z1s * sqrt(phat1s * (1 - phat1s)/(nhat1s))), 1)
            } else {
              stop("Argument alt must be specified!")
            }
          }
        }
      }
      if (met == "Clopper-Pearson") {
        if (alt == "two.sided") {
          if (y != 0) {
            cilower <- qbeta((1 - conf.level)/2, y, n - y + 1)
          }
          if (y != n) {
            ciupper <- qbeta(1 - (1 - conf.level)/2, y + 1, n - y)
          }
        }
        if (alt == "less") {
          if (y != n) {
            ciupper <- qbeta(1 - (1 - conf.level), y + 1, n - y)
          }
        }
        if (alt == "greater") {
          if (y != 0) {
            cilower <- qbeta((1 - conf.level), y, n - y + 1)
          }
        }
        CI <- c(cilower, ciupper)
      }
      if (met == "Blaker") {
        step.value = 1e-04
        blakeraccept <- function(y, n, p) {
          p1 = 1 - pbinom(y - 1, n, p)
          p2 = pbinom(y, n, p)
          a1 = p1 + pbinom(qbinom(p1, n, p) - 1, n, p)
          a2 = p2 + 1 - pbinom(qbinom(1 - p2, n, p), n, p)
          return(min(a1, a2))
        }
        if (alt == "two.sided") {
          if (y != 0) {
            cilower <- qbeta((1 - conf.level)/2, y, n - (y) + 1)
            {
              while (blakeraccept(y, n, cilower + step.value) < (1 - 
                                                                 conf.level)) cilower = cilower + step.value
            }
          }
          if (y != n) {
            ciupper <- qbeta(1 - (1 - conf.level)/2, (y) + 1, n - (y))
            {
              while (blakeraccept(y, n, ciupper - step.value) < (1 - 
                                                                 conf.level)) ciupper = ciupper - step.value
            }
          }
        }
        if (alt == "less") {
          if (y != n) {
            ciupper <- qbeta(1 - (1 - conf.level), (y) + 1, n - (y))
            {
              while (blakeraccept(y, n, ciupper - step.value) < (1 - 
                                                                 conf.level)) ciupper = ciupper - step.value
            }
          }
        }
        if (alt == "greater") {
          if (y != 0) {
            cilower <- qbeta((1 - conf.level), y, n - (y) + 1)
            {
              while (blakeraccept(y, n, cilower + step.value) < (1 - 
                                                                 conf.level)) cilower = cilower + step.value
            }
          }
        }
        CI <- c(cilower, ciupper)
      }
      return(CI)
    }
    
    # Errors
    if (r1_kprev > r1_nprev | r2_kprev > r2_nprev) {
      stop("k must be <= to n")
    }
    if (Se1_const > 1 | Se1_const < 0) {
      stop("Please select the value of p_alt to be between 0 and 1")
    }
    if (Se2_const > 1 | Se2_const < 0) {
      stop("Please select the value of p_alt to be between 0 and 1")
    }
    if (Sp1_const > 1 | Sp1_const < 0) {
      stop("Please select the value of p_alt to be between 0 and 1")
    }
    if (Sp2_const > 1 | Sp2_const < 0) {
      stop("Please select the value of p_alt to be between 0 and 1")
    }
    
    alt <- match.arg(alt, choices = c("two.sided", "less", "greater"))
    met <- match.arg(met, choices = c("Wald", "Wilson", "Agresti-Coull", "Clopper-Pearson", "Blaker"))
    
    obs.prev1 <- r1_kprev / r1_nprev
    obs.prev2 <- r2_kprev / r2_nprev
    
    r1_hat <- (obs.prev1+Sp1_const-1)/(Se1_const+Sp1_const-1)
    r1_hat <- min(1, max(0, r1_hat))
    r2_hat <- (obs.prev2+Sp2_const-1)/(Se2_const+Sp2_const-1)
    r2_hat <- min(1, max(0, r2_hat))
    
    l1_u1 <- ((BinomialConfInt(n = r1_nprev, y = r1_kprev, conf.level = conflevel, alt = "two.sided", met = met) + Sp1_const - 1) / (Se1_const + Sp1_const - 1))
    l2_u2 <- ((BinomialConfInt(n = r2_nprev, y = r2_kprev, conf.level = conflevel, alt = "two.sided", met = met) + Sp2_const - 1) / (Se2_const + Sp2_const - 1))
    
    CI_lower <- r1_hat - r2_hat - sqrt((r1_hat-l1_u1[1])^2+(l2_u2[2]-r2_hat)^2) # Lower Limit for the CI of the RD
    CI_upper <- r1_hat - r2_hat + sqrt((l1_u1[2]-r1_hat)^2+(r2_hat-l2_u2[1])^2) # Upper Limit for the CI of the RD
    
    ci = c(CI_lower, CI_upper)
    
    # Return
    return(ci)
  }
  
  # Only works for two-sided tests and n1=n2
  # true_p2 must be higher than true_p1
  
  # Errors
  if (Se1 > 1 | Se1 < 0) {
    stop("Please select the value of p_alt to be between 0 and 1")
  }
  if (Se2 > 1 | Se2 < 0) {
    stop("Please select the value of p_alt to be between 0 and 1")
  }
  if (Sp1 > 1 | Sp1 < 0) {
    stop("Please select the value of p_alt to be between 0 and 1")
  }
  if (Sp2 > 1 | Sp2 < 0) {
    stop("Please select the value of p_alt to be between 0 and 1")
  }
  
  met <- match.arg(met, choices = c("Wald", "Wilson", "Agresti-Coull", "Clopper-Pearson", "Blaker"))
  
  conf=1-alpha
  
  obsp1 = Se1*true_p1 + (1-Sp1)*(1-true_p1) 
  obsp2 = Se2*true_p2 + (1-Sp2)*(1-true_p2) 
  obsp1 = min(1,max(0,obsp1))
  obsp2 = min(1,max(0,obsp2))
  
  nn1=min_n; nn2=max_n
  
  repeat{
    nn=floor((nn1+nn2)/2)
    elut=rep(NA,R)
    k1=rbinom(R,nn,obsp1)
    k2=rbinom(R,nn,obsp2)
    for (r in 1:R){
      #cat(nn,k1[r],k2[r],"\n")
      CI=CI_Reiczigel2010_Zou_RD_FixSeSp(nn,k1[r],Se1,Sp1,nn,k2[r],Se2,Sp2,conf=conf,met=met)
      if(CI[1]<=0 & 0<=CI[2]) elut[r]=0 else elut[r]=1
    }
    ero=mean(elut)
    #cat(nn,ero,"\n")
    if(ero >= pwr) nn2=nn else nn1=nn
    if(nn2-nn1<3)break
  }
  
  #cat("---\n")
  n=nn1
  firstn=0	
  for (n in nn:max_n){
    elut=rep(NA,R)
    k1=rbinom(R,n,obsp1)
    k2=rbinom(R,n,obsp2)
    for (r in 1:R){
      CI=CI_Reiczigel2010_Zou_RD_FixSeSp(n,k1[r],Se1,Sp1,n,k2[r],Se2,Sp2,conf=conf,met="Wil")
      if(CI[1]<=0 & 0<=CI[2]) elut[r]=0 else elut[r]=1
    }
    ero=mean(elut)
    #cat(n,ero,firstn,"\n")
    if(ero < pwr) firstn=0 
    if(ero >= pwr & firstn==0) firstn=n 
    if(ero >= pwr & firstn>0 & firstn <= n*(1-drprate)) return(n)
  }
  if (-firstn == 0) {res <- "> Maximum Sample Size"} else {
    res <- -firstn
  }
  return(res)
}

SampleSizetoPower2sampleFixSeSp(
  0.4,
  0.25,
  0.998,
  0.996,
  0.998,
  0.996,
  0.05,
  0.8,
  10,
  4000,
  0.2,
  R = 20000,
  met = "Wi")

p1=c(.01,.02,.03,.05,.10,.20,.30)
p2=c(.10,.12,.13,.15,.20,.35,.50)
for (i in 1:7)
  for (Se1 in c(1,.99,.98,.95, 0.9)) for (Sp1 in c(1,.99,.98,.95, 0.9)) 
    for (Se2 in c(1,.99,.98,.95, 0.9)) for (Sp2 in c(1,.99,.98,.95, 0.9))
      cat("p1:",p1[i],"p2:",p2[i],"Se1:",Se1,"Sp1:",Sp1,"Se2:",Se2,"Sp2:",Sp2,"Sample Size: ",
          sampleSizeReiczigelRD(p1[i], p2[i],  Se1,  Sp1,  Se2,  Sp2,  .05, .8,   10, 4000,    .15, 20000, "Wil"),"\n", 
          file = "interimresults_2sample_FixSeSp-SS.txt", append = TRUE)