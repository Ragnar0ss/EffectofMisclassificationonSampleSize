SampleSizetoPower1sampleEstSeSp <- function(p_null, 
                                               alpha = 0.05, 
                                               pwr = 0.8, 
                                               p_alt, 
                                               min_n = 20, 
                                               max_n, 
                                               drp_rate = 0,
                                               n_se,
                                               Se_exp,
                                               n_sp,
                                               Sp_exp, 
                                               R = 20000) {
  
  # Function to calculate the required sample size for a binomial test when Se and Sp are estimated from validation samples.
  #  
  # Args:
  #   p_null: Parameter value under the null hypothesis.
  #   alpha: Type I error probability.
  #   pwr: Prescribed power of the test. Specify as a number between 0 and 1 rather than a percent value. 
  #   p_alt: Hypothetized true proportion for which the prescribed power "pwr" should be reached (must conform to the argument "alt")
  #   min_n: Minimum sample size (smaller sample sizes will not be examined).
  #   max_n: Maximum sample size (if power is below pwr for all n less than or equal to max_n, function returns 0).
  #   drp_rate: Highest drop-out rate for which power must not fall below pwr. Specify as a number between 0 and 1 rather than a percent value. 
  #   n_se: Sample size for a sensitivity validation study.
  #   Se_exp: Hypothesized sensitivity of the diagnostic test.
  #   n_sp: Sample size for a specificity validation study.
  #   Sp_exp: Hypothesized specificity of the diagnostic test.
  #   R: Number of simulations to be used.
  #
  # Returns:
  #   The minimal sample size n satisfying that for all sample sizes from drp_rate*n to n power is at least pwr.
  
  # Nested function:
  
  BinomialPowerEstSeSp_LangReiczigel <- function(n,
                                                 p_null,
                                                 p_alt,
                                                 conf.level = 0.95,
                                                 Se_exp,
                                                 n_se,
                                                 Sp_exp,
                                                 n_sp,
                                                 R = 20000) {
    
    # Function to calculate the simulated power of a binomial test when Se and Sp are not fixed but estimated from validation samples.
    #  
    # Args:
    #   n: The number of individuals that are tested.
    #   p_null: Parameter value under the null hypothesis.
    #   p_alt: Parameter value under the alternative hypothesis.
    #   conf.level: Confidence level.
    #   Se_exp: Hypothesized sensitivity of the diagnostic test.
    #   n_se: Sample size for a sensitivity validation study.
    #   Sp_exp: Hypothesized specificity of the diagnostic test.
    #   n_sp: Sample size for a specificity validation study.
    #   R: Number of simulations.
    #
    # Returns:
    #   The power of the specified test.
    
    CI_Binom = function(
    nprev,       # Sample size for prevalence
    kprev,       # Frequency of positive diagnoses in sample of size nprev
    nsens,       # Sample size for sensitivity
    ksens,       # Frequency of positive diagnoses in sample of size nsens
    nspec,       # Sample size for specificity
    kspec,       # Frequency of negative diagnoses in sample of size nspec
    conflevel=.95) # Confidence level
    {
      # Observed relative frequencies
      obs.prev = kprev/nprev
      obs.sens = ksens/nsens
      obs.spec = kspec/nspec
      
      # Rogan-Gladen point estimate of true prevalence
      est.prev = (obs.prev+obs.spec-1)/(obs.sens+obs.spec-1)
      est.prev = min(1,max(0,est.prev))
      
      # Adjustments
      zcrit = qnorm((1+conflevel)/2)
      plus  = 2
      
      nprev. = nprev+zcrit^2
      kprev. = kprev+zcrit^2/2
      
      nsens. = nsens+plus
      ksens. = ksens+plus/2
      
      nspec. = nspec+plus
      kspec. = kspec+plus/2
      
      obs.prev. = kprev./nprev.
      obs.sens. = ksens./nsens.
      obs.spec. = kspec./nspec.
      
      est.prev. = (obs.prev.+obs.spec.-1)/(obs.sens.+obs.spec.-1)
      
      # Youden index
      Youden. = obs.sens.+obs.spec.-1
      
      # Standard error of est.prev.
      se.est.prev. = sqrt(
        obs.prev.*(1-obs.prev.)/nprev. +
          obs.sens.*(1-obs.sens.)/nsens. * est.prev.^2 +
          obs.spec.*(1-obs.spec.)/nspec. * (1-est.prev.)^2
      )/abs(Youden.)
      
      # Shift parameter
      dprev = 2*zcrit^2*
        (est.prev.*obs.sens.*(1-obs.sens.)/nsens. - (1-est.prev.)*obs.spec.*(1-obs.spec.)/nspec.)
      
      # Adjusted confidence limits
      LCL = est.prev.+dprev - zcrit*se.est.prev.
      UCL = est.prev.+dprev + zcrit*se.est.prev.
      LCL = min(1,max(0,LCL))
      UCL = min(1,max(0,UCL))
      
      ci <- c(LCL, UCL)
      return(ci)
      # cat("\nSensitivity: ", round(obs.sens, 4), ", adjusted: ", round(obs.sens., 4),
      #     "\nSpecificity: ", round(obs.spec, 4), ", adjusted: ", round(obs.spec., 4),
      #     "\nObserved prevalence: ", round(obs.prev, 8), ", adjusted: ", round(obs.prev., 8),
      #     "\nRogan-Gladen true prevalence: ", round(est.prev,4), "\n",
      #     "Adjusted ", round(100*conflevel), "% CI: ", round(LCL,4), " - ", round(UCL,4), "\n",
      #     sep="")
    }
    
    # Warnings
    if (p_alt > 1 | p_alt < 0) {
      stop("Please select the value of p_alt to be between 0 and 1")
    }
    if (p_null > 1 | p_null < 0) {
      stop("Please select the value of p_null to be between 0 and 1")
    }
    if (Se_exp > 1 | Se_exp < 0) {
      stop("Please select the value of Se_exp to be between 0 and 1")
    }
    if (Sp_exp > 1 | Sp_exp < 0) {
      stop("Please select the value of Sp_exp to be between 0 and 1")
    }
    
    p_obs <- (Se_exp * p_alt) + ((1 - Sp_exp) * (1 - p_alt))
    p_obs = min(1,max(0,p_obs))
    
    elut=rep(NA,R)
    
    for (i in 1:R) {
      k_p <- rbinom(1, n, p_obs)
      k_se <- rbinom(1, n_se, Se_exp)
      k_sp <- rbinom(1, n_sp, Sp_exp)
      
      CI <- CI_Binom(nprev = n, kprev = k_p, nsens = n_se, ksens = k_se,
                     nspec = n_sp, kspec = k_sp, conflevel = conf.level)
      
      if(CI[1]<=p_null & p_null<=CI[2]) elut[i]=0 else elut[i]=1
    }
    ero=mean(elut)
    
    return(ero)
  }
  
  # Initial parameter values:
  
  conf <- 1 - alpha
  p_obs <- (Se_exp * p_alt) + ((1 - Sp_exp) * (1 - p_alt))
  p_obs <- min(1, max(p_obs,0))

  # Errors and Warnings
  if (p_alt > 1 | p_alt < 0) {
    stop("Please select the value of p_alt to be between 0 and 1")
  }
  if (p_null > 1 | p_null < 0) {
    stop("Please select the value of p_null to be between 0 and 1")
  }
  if (pwr > 1 | pwr < 0) {
    stop("Please select the value of pwr to be between 0 and 1")
  }
  if (drp_rate > 1 | drp_rate < 0) {
    stop("Please select the value of drp_rate to be between 0 and 1")
  }
  if (p_null == p_alt) {
    stop("Power is NA as p_null = p_alt!")
  }
  if (Se_exp > 1 | Se_exp < 0) {
    stop("Please select the value of Se_exp to be between 0 and 1")
  }
  if (Sp_exp > 1 | Sp_exp < 0) {
    stop("Please select the value of Sp_exp to be between 0 and 1")
  }
  
  first.good.n <- 0
  act_alpha <- 0
  act_pwr <- 1
  val <- min_n
  
  repeat{
    ntopwr <- BinomialPowerEstSeSp_LangReiczigel(n = val, p_null = p_null, p_alt = p_alt, conf.level = conf,
                                                 Se_exp = Se_exp, n_se = n_se, Sp_exp = Sp_exp, n_sp = n_sp, R = R)
    if(ntopwr < pwr){
      val <- val+1
      if(val > max_n) {
        return("> Maximum Sample Size")
      }
      first.good.n <- 0
      act_pwr <- 1
      next
    }	
    if(first.good.n == 0) first.good.n <- val else if(first.good.n <= val*(1-drp_rate)) break  
    val <- val + 1
    if(val > max_n) {
      return("> Maximum Sample Size")
    }
  }
  
  ## Return options
  return(val)
}

SampleSizetoPower1sampleEstSeSp(p_null = 0.02,
                                 alpha = 0.05, 
                                 pwr = 0.8, 
                                 p_alt = 0.002, 
                                 min_n = 20, 
                                 max_n = 100000, 
                                 drp_rate = 0, 
                                 Se_exp = 0.9,
                                 n_se = 1530,
                                 Sp_exp = 0.9659, 
                                 n_sp = 3230,
                                 R = 200)

