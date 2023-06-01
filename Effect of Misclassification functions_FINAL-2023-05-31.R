SampleSizetoPowerSeSp <- function(p_null, 
                                  alt = "two.sided", 
                                  met = "Clopper-Pearson", 
                                  alpha = 0.05, 
                                  max_alpha, 
                                  pwr = 0.8, 
                                  p_alt, 
                                  min_n = 5, 
                                  max_n, 
                                  drp_rate = 0, 
                                  Se = 1, 
                                  Sp = 1) {
  
  # Function to calculate the required sample size for a binomial test.
  #  
  # Args:
  #   p_null: Parameter value under the null hypothesis.
  #   alt: Alternative - "two.sided", "less", "greater".
  #   met: Method - Available methods: "Wald", "Wilson", "Agresti-Coull", "Clopper-Pearson", "Blaker"
  #   alpha: Type I error probability.
  #   max_alpha: Maximum allowed Type I error rate for asymptotic tests (sample sizes where actual alpha exceeds max_alpha are not accepted)
  #   pwr: Prescribed power of the test. Specify as a number between 0 and 1 rather than a percent value. 
  #   p_alt: Hypothetized true proportion for which the prescribed power "pwr" should be reached (must conform to the argument "alt")
  #   min_n: Minimum sample size (smaller sample sizes will not be examined).
  #   max_n: Maximum sample size (if power is below pwr for all n less than or equal to max_n, function returns 0).
  #   drp_rate: Highest drop-out rate for which power must not fall below pwr. Specify as a number between 0 and 1 rather than a percent value. 
  #   Se: Sensitivity of the diagnostic test.
  #   Sp: Specificity of the diagnostic test.
  #
  # Returns:
  #   The minimal sample size n satisfying that for all sample sizes from drp_rate*n to n power is at least pwr, and actual alpha is at most max_alpha.
  
  # Nested functions:
  
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
    #   met: Method - Available methods: "Wald", "Wilson", "Agresti-Coull", "Clopper-Pearson", "Blaker"
    #
    # Returns:
    #   A vector containing the LCL and UCL of the binomial CI.
    
    p_obs <- y / n
    alpha <- 1 - conf.level
    z1s <- qnorm(1 - alpha)
    z2s <- qnorm(1 - alpha / 2)
    cilower <- 0
    ciupper <- 1
    
    alt <- match.arg(alt, choices = c("two.sided", "less", "greater"))
    met <- match.arg(met, choices = c("Wald", "Wilson", "Agresti-Coull", "Clopper-Pearson", "Blaker"))
    
    # Warnings
    if (alt != "two.sided" & alt != "less" & alt != "greater") {
      stop("Please select alt from two.sided, less or greater.")
    }
    if (met != "Wald" & met != "Wilson" & met != "Agresti-Coull" & met != "Clopper-Pearson" & met != "Blaker") {
      stop("Please select met from the following: Wald, Wilson, Agresti-Coull, Clopper-Pearson and Blaker")
    }
    
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
  
  
  BinomialPowerSeSp <- function(n,
                                p_null,
                                p_alt,
                                conf.level = 0.95,
                                alt = "two.sided",
                                met = "Clopper-Pearson",
                                Se = 1,
                                Sp = 1) {
    
    # Function to calculate the power of a binomial test
    #  
    # Args:
    #   n: The number of individuals that are tested.
    #   y: Number of individuals with positive test results.
    #   conf.level: Confidence level.
    #   alt: alternative - "two.sided", "less", "greater".
    #   met: Method - Available methods: "Wald", "Wilson", "Agresti-Coull", "Clopper-Pearson", "Blaker"
    #   Se: Sensitivity of the diagnostic test.
    #   Sp: Specificity of the diagnostic test.
    #
    # Returns:
    #   The power of the specified test.
    
    p_obs <- (Se * p_alt) + ((1 - Sp) * (1 - p_alt))
    yvec <- 0:n
    
    alt <- match.arg(alt, choices = c("two.sided", "less", "greater"))
    met <- match.arg(met, choices = c("Wald", "Wilson", "Agresti-Coull", "Clopper-Pearson", "Blaker"))
    
    # Warnings
    if (alt != "two.sided" & alt != "less" & alt != "greater") {
      stop("Please select alt from two.sided, less or greater.")
    }
    if (met != "Wald" & met != "Wilson" & met != "Agresti-Coull" & met != "Clopper-Pearson" & met != "Blaker") {
      stop("Please select met from the following: Wald, Wilson, Agresti-Coull, Clopper-Pearson and Blaker")
    }
    if (p_alt > 1 | p_alt < 0) {
      stop("Please select the value of p_alt to be between 0 and 1")
    }
    if (p_null > 1 | p_null < 0) {
      stop("Please select the value of p_null to be between 0 and 1")
    }
    
    # Alternatives
    if (alt == "less") {
      powvec <- numeric(length = length(yvec))
      probvec <- numeric(length = length(yvec))
      for (i in 1:length(yvec)) {
        CI <- ((BinomialConfInt(y = yvec[i], n = n, conf.level = conf.level, alt = alt, met = met) + Sp - 1) / (Se + Sp - 1))
        powvec[i] <- CI[1] > p_null | CI[2] < p_null
        probvec[i] <- dbinom(x = yvec[i], size = n, prob = p_obs)
      }
      power = sum(powvec * probvec)
    }
    if (alt == "greater") {
      powvec <- numeric(length = length(yvec))
      probvec <- numeric(length = length(yvec))
      for (i in 1:length(yvec)) {
        CI <- ((BinomialConfInt(y = yvec[i], n = n, conf.level = conf.level, alt = alt, met = met) + Sp - 1) / (Se + Sp - 1))
        powvec[i] <- CI[1] > p_null | CI[2] < p_null
        probvec[i] <- dbinom(x = yvec[i], size = n, prob = p_obs)
      }  
      power = sum(powvec * probvec)
    }
    if (alt == "two.sided") {
      powvec <- numeric(length = length(yvec))
      probvec <- numeric(length = length(yvec))
      for (i in 1:length(yvec)) {
        CI <- ((BinomialConfInt(y = yvec[i], n = n, conf.level = conf.level, alt = alt, met = met) + Sp - 1) / (Se + Sp - 1))
        powvec[i] <- CI[1] > p_null | CI[2] < p_null
        probvec[i] <- dbinom(x = yvec[i], size = n, prob = p_obs)
      }
      power = sum(powvec * probvec)
    }
    
    #out <- list(power = power)
    return(power)
  }
  
  # Initial parameter values:
  
  conf <- 1 - alpha 
  ntopwr <- 0
  act_alpha <- 0
  meta <- 0
  val <- min_n
  alpha_vec <- 0
  act_power <- 0
  pwr_vec <- 0
  res <- 0
  
  alt <- match.arg(alt, choices = c("two.sided", "less", "greater"))
  met <- match.arg(met, choices = c("Wald", "Wilson", "Agresti-Coull", "Clopper-Pearson", "Blaker"))
  
  # Errors and Warnings
  if (alt != "two.sided" & alt != "less" & alt != "greater") {
    stop("Please select alt from two.sided, less or greater.")
  }
  if (met != "Wald" & met != "Wilson" & met != "Agresti-Coull" & met != "Clopper-Pearson" & met != "Blaker") {
    stop("Please select met from the following: Wald, Wilson, Agresti-Coull, Clopper-Pearson and Blaker")
  }
  if (p_null < p_alt & alt == "less") {
    stop("Power is NA as p_null < p_alt!")
  }
  if (p_null > p_alt & alt == "greater") {
    stop("Power is NA as p_null > p_alt!")
  }
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
  
  if (alt == "two.sided") {
      if (p_null == p_alt) {
      } else {
		first.good.n <- 0
		act_alpha <- 0
		act_pwr <- 1
		
        repeat{
          ntopwr <- BinomialPowerSeSp(n = val, p_null = p_null, p_alt = p_alt, 
                                      conf.level = conf, alt = alt, met = met, Se = Se, Sp = Sp)
          if(ntopwr < pwr){
			val <- val+1
			first.good.n <- 0
			act_alpha <- 0
			act_pwr <- 1
			next
		  }	
		  
#		  ntoalpha <- BinomialPowerSeSp(n = val, p_null = p_null, p_alt = p_null, 
#                                                         conf.level = conf, alt = alt, met = met, Se = Se, Sp = Sp)
#         if(ntoalpha > max_alpha){
#			val <- val+1
#			first.good.n <- 0
#			act_alpha <- 0
#			act_pwr <- 1
#			next
#		  }	
#		  
#		  if(ntopwr < act_pwr) act_pwr <- ntopwr
#		  if(ntoalpha > act_alpha) act_alpha <- ntoalpha
		  if(first.good.n == 0) first.good.n <- val else if(first.good.n <= val*(1-drp_rate)) break  
		  val <- val + 1
        }
      }
    }
  if (alt == "less") {
      if (p_null == p_alt | p_null < p_alt) {
      } else {
		first.good.n <- 0
		act_alpha <- 0
		act_pwr <- 1
		
        repeat{
          ntopwr <- BinomialPowerSeSp(n = val, p_null = p_null, p_alt = p_alt, 
                                      conf.level = conf, alt = alt, met = met, Se = Se, Sp = Sp)
          if(ntopwr < pwr){
			val <- val+1
			first.good.n <- 0
			act_alpha <- 0
			act_pwr <- 1
			next
		  }	
		  
#		  ntoalpha <- BinomialPowerSeSp(n = val, p_null = p_null, p_alt = p_null, 
#                                                         conf.level = conf, alt = alt, met = met, Se = Se, Sp = Sp)
#         if(ntoalpha > max_alpha){
#			val <- val+1
#			first.good.n <- 0
#			act_alpha <- 0
#			act_pwr <- 1
#			next
#		  }	
#		  
#		  if(ntopwr < act_pwr) act_pwr <- ntopwr
#		  if(ntoalpha > act_alpha) act_alpha <- ntoalpha
		  if(first.good.n == 0) first.good.n <- val else if(first.good.n <= val*(1-drp_rate)) break  
		  val <- val + 1
        }
      }
    }
  if (alt == "greater") {
      if (p_null == p_alt | p_null > p_alt) {
      } else {
		first.good.n <- 0
		act_alpha <- 0
		act_pwr <- 1
		
        repeat{
          ntopwr <- BinomialPowerSeSp(n = val, p_null = p_null, p_alt = p_alt, 
                                      conf.level = conf, alt = alt, met = met, Se = Se, Sp = Sp)
          if(ntopwr < pwr){
			val <- val+1
			first.good.n <- 0
			act_alpha <- 0
			act_pwr <- 1
			next
		  }	
		  
#		  ntoalpha <- BinomialPowerSeSp(n = val, p_null = p_null, p_alt = p_null, 
#                                                        conf.level = conf, alt = alt, met = met, Se = Se, Sp = Sp)
#         if(ntoalpha > max_alpha){
#		    val <- val+1
#			first.good.n <- 0
#			act_alpha <- 0
#			act_pwr <- 1
#			next
#	  	  }	
	  
#		  if(ntopwr < act_pwr) act_pwr <- ntopwr
#		  if(ntoalpha > act_alpha) act_alpha <- ntoalpha
		  if(first.good.n == 0) first.good.n <- val else if(first.good.n <= val*(1-drp_rate)) break  
		  val <- val + 1
		    }
      }
    }
  res <- val
  # Return options
   cat("Se:",Se,"   Sp:",Sp,"   pnull:",p_null,"    alt:",alt,"     ptrue:",p_alt,"    met:",met,"    drop out rate:"
       ,drp_rate, "    Sample size:",res, "\n", file = "interimresults.txt", append = TRUE)
  # cat("Se:",Se,"   Sp:",Sp,"   pnull:",p_null,"    alt:",alt,"     ptrue:",p_alt,"    met:",met,"    drop out rate:"
  #    ,drp_rate, "    Sample size:",res)
  # return(res)
 }
 
SampleSizetoPowerSeSp(p_null=.1, 
                      alt = "two", 
                      met = "Cl", 
                      alpha = 0.05, 
                      max_alpha=.06, 
                      pwr = 0.9, 
                      p_alt=0.1, 
                      min_n = 20, 
                      max_n=500, 
                      drp_rate = 0.3, 
                      Se = 1, 
                      Sp = 1)		  
