SampleSizetoPower2sampleEstSeSp = function(
    true_p1,
    true_p2,
    Se1_n,
    Se1_exp,
    Sp1_n,
    Sp1_exp,
    Se2_n,
    Se2_exp,
    Sp2_n,
    Sp2_exp,
    alpha = 0.05,
    pwr = 0.8,
    min_n = 20,
    max_n = 5000,
    drprate = 0,
    R = 20000){
  
  # Function to calculate the required sample size for a two-sample binomial test, when Se1, Se2, Sp1 and Sp2 are estimated.
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
  
  CI_Reiczigel2010_Zou_RD_EstSeSp <- function(
    r1_nprev,
    r1_kprev,
    Se1_nprev,
    Se1_kprev,
    Sp1_nprev,
    Sp1_kprev,
    r2_nprev,
    r2_kprev,
    Se2_nprev,
    Se2_kprev,
    Sp2_nprev,
    Sp2_kprev,
    conflevel = 0.95)
  {
    # Function to calculate the Confidence Interval for for the Risk Difference between two groups by combining the methods of Lang and Reiczigel (2014) and Zou and Donner (2008).
    #  
    # Args:
    #   r1_nprev: The number of individuals that are tested in group 1.
    #   r1_kprev: Number of individuals with positive test results in group 1.
    #   Se1_nprev: The size of the sample used for the estimation of the Sensitivity of the test used for testing in group 1.
    #   Se1_kprev: The positive cases of the validation study for the estimation of the Sensitivity of the test used for testing in group 1.
    #   Sp1_nprev: The size of the sample used for the estimation of the Specificity of the test used for testing in group 1.
    #   Sp1_kprev: The positive cases of the validation study for the estimation of the Specificity for the test used for testing in group 1.
    #   r2_nprev: The number of individuals that are tested in group 1.
    #   r2_kprev: Number of individuals with positive test results in group 1.
    #   Se2_nprev: The size of the sample used for the estimation of the Sensitivity of the test used for testing in group 2.
    #   Se2_kprev: The positive cases of the validation study for the estimation of the Sensitivity for the test used for testing in group 2.
    #   Sp2_nprev: The size of the sample used for the estimation of the Specificity of the test used for testing in group 1.
    #   Sp2_kprev: The positive cases of the validation study for the estimation of the Specificity of the test used for testing in group 2.
    #   conf.level: Confidence level.
    #   alt: alternative - "two.sided", "less", "greater".
    #
    # Returns:
    #   A vector containing the LCL and UCL of the binomial CI.
    
    # Nested function:
    
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
    
    # Errors
    if (r1_kprev > r1_nprev | r2_kprev > r2_nprev) {
      stop("k must be <= to n")
    }
    if (Se1_kprev > Se1_nprev | Se2_kprev > Se2_nprev) {
      stop("k must be <= to n")
    }
    if (Sp1_kprev > Sp1_nprev | Sp2_kprev > Sp2_nprev) {
      stop("k must be <= to n")
    }
    
    obs.prev1 <- r1_kprev / r1_nprev
    obs.prev2 <- r2_kprev / r2_nprev
    
    r1_hat <- (obs.prev1+(Sp1_kprev/Sp1_nprev)-1)/((Se1_kprev/Se1_nprev)+(Sp1_kprev/Sp1_nprev)-1)
    r1_hat <- min(1, max(0, r1_hat))
    r2_hat <- (obs.prev2+(Sp2_kprev/Sp2_nprev)-1)/((Se2_kprev/Se2_nprev)+(Sp2_kprev/Sp2_nprev)-1)
    r2_hat <- min(1, max(0, r2_hat))
    
    l1_u1 <- CI_Binom(kprev = r1_kprev, nprev = r1_nprev, 
                      ksens = Se1_kprev, nsens = Se1_nprev,
                      kspec = Sp1_kprev, nspec = Sp1_nprev,
                      conflevel = conflevel)
    l2_u2 <- CI_Binom(kprev = r2_kprev, nprev = r2_nprev, 
                      ksens = Se2_kprev, nsens = Se2_nprev,
                      kspec = Sp2_kprev, nspec = Sp2_nprev,
                      conflevel = conflevel)
    
    CI_lower <- r1_hat - r2_hat - sqrt((r1_hat-l1_u1[1])^2+(l2_u2[2]-r2_hat)^2) # Lower Limit for the CI of the RD
    CI_upper <- r1_hat - r2_hat + sqrt((l1_u1[2]-r1_hat)^2+(r2_hat-l2_u2[1])^2) # Upper Limit for the CI of the RD
    
    ci = c(CI_lower, CI_upper)
    
    # Returns
    
    return(ci)
  }
  
  # Only works for two-sided tests and n1=n2
  # true_p2 must be higher than true_p1
  
  conf=1-alpha
  
  obsp1 = (Se1_exp)*true_p1 + (1-(Sp1_exp))*(1-true_p1) 
  obsp2 = (Se2_exp)*true_p2 + (1-(Sp2_exp))*(1-true_p2)
  obsp1 = min(1,max(0,obsp1))
  obsp2 = min(1,max(0,obsp2))
  
  nn1=min_n; nn2=max_n
  
  repeat{
    nn=floor((nn1+nn2)/2)
    elut=rep(NA,R)
    k1=rbinom(R,nn,obsp1)
    k2=rbinom(R,nn,obsp2)
    se1_k=rbinom(R,Se1_n,Se1_exp)
    se2_k=rbinom(R,Se2_n,Se2_exp)
    sp1_k=rbinom(R,Sp1_n,Sp1_exp)
    sp2_k=rbinom(R,Sp2_n,Sp2_exp)
    
    for (r in 1:R){
      #cat(nn,k1[r],k2[r],"\n")
      CI = CI_Reiczigel2010_Zou_RD_EstSeSp(
        r1_nprev = nn,
        r1_kprev = k1[r],
        Se1_nprev = Se1_n,
        Se1_kprev = se1_k[r],
        Sp1_nprev = Sp1_n,
        Sp1_kprev = sp1_k[r],
        r2_nprev = nn,
        r2_kprev = k2[r],
        Se2_nprev = Se2_n,
        Se2_kprev = se2_k[r],
        Sp2_nprev = Sp2_n,
        Sp2_kprev = sp2_k[r],
        conflevel = conf
      )
      #cat(CI, "\n")
      if(CI[1]<=0 & 0<=CI[2]) elut[r]=0 else elut[r]=1
      #cat(elut[r])
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
      CI = CI_Reiczigel2010_Zou_RD_EstSeSp(
        r1_nprev = n,
        r1_kprev = k1[r],
        Se1_nprev = Se1_n,
        Se1_kprev = se1_k[r],
        Sp1_nprev = Sp1_n,
        Sp1_kprev = sp1_k[r],
        r2_nprev = n,
        r2_kprev = k2[r],
        Se2_nprev = Se2_n,
        Se2_kprev = se2_k[r],
        Sp2_nprev = Sp2_n,
        Sp2_kprev = sp2_k[r],
        conflevel = conf
      )
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

SampleSizetoPower2sampleEstSeSp(true_p1 = 0.4,
                                true_p2 = 0.25,
                                Se1_n = 425,
                                Se1_exp = 0.998,
                                Sp1_n = 1612,
                                Sp1_exp = 0.996,
                                Se2_n = 425,
                                Se2_exp = 0.998,
                                Sp2_n = 1612,
                                Sp2_exp = 0.996,
                                alpha = 0.05,
                                pwr = 0.8,
                                min_n = 20,
                                max_n = 5000,
                                drprate = 0.20,
                                R = 20000)


