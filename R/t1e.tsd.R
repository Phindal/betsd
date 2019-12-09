# ------------------------------------------------------------------------------
# Function to adjust the significance levels for bioequivalence studies using
# adaptive two-stage 2x2 crossover designs.
#
# Authors: E. Molins, D. Labes, H. Schuetz, J. Ocaña
# ------------------------------------------------------------------------------
#ADJUST ALPHA1 AND ALPHA2 ON A TSD 2x2 CXO DESIGN, CONTROLLING T1E BELOW ALPHA
t1e.tsd <- function(n1, CV, GMR = 0.95, Nmax = 150, min.n2 = n1/2, type = 1,
                    alpha = 0.05, alpha1, alpha2, targetpower = 0.8,
                    setseed = TRUE, theta1, theta2,
                    details = TRUE, print = TRUE,
                    ...) {

  #Non-modified parameters => Nmax == Inf & min.n2 == 0
  #Modified CRITERIA ACCORDING TO DOI: 10.1002/sim.7452
  min.n2 <- sapply(min.n2, function(y) if (y %% 2 != 0) y+y%%2 else y)

  #DEBUG AT FIRST ITERATION
  if (missing(n1))
    stop("Number of subjects in stage 1 must be given")
  if (any(n1 < 12))
    stop("Number of subjects in stage 1 must be at least 12")
  if (missing(CV))
    stop("CV must be given")
  if (any(CV <= 0.05)) #This is 0.05 because the algorithm will look at the given CV +/- 0.05, and CV must be > 0
    stop("CV must be > 0.05 becuase the algorithm will look at CV +/- 0.05")
  if (length(alpha) > 1)
    stop("alpha must be of length 1")
  if (missing(alpha1) & missing(alpha2)) {
    alpha1 <- 0.0294
    alpha2 <- alpha1
  }
  if (missing(alpha1) & !missing(alpha2)) alpha1 <- alpha2
  if (!missing(alpha1) & missing(alpha2)) alpha2 <- alpha1
  if (alpha1 && alpha2 > 0.05) 
    stop("apha values must be below 0.05 (default at 0.0294)")
  if (length(type) > 1)
    stop("type must be of length 1")
  if (type!=1 & type!=2)
    stop("type must be 1 or 2")
  if (type == 1) {type <- "B"} else {type <- "C"} #type 1: Potvin "B"; type 2: Potvin "C or D"
  if (missing(theta1) & missing(theta2))  theta1 <- 0.8
  if (!missing(theta1) & missing(theta2)) theta2 <- 1/theta1
  if (missing(theta1) & !missing(theta2)) theta1 <- 1/theta2
  if (GMR <= theta1 | GMR >= theta2)
    stop("GMR must be within acceptance range")

  #SETSEED
  if (is.numeric(setseed)) {
    setseed <- TRUE
  } else if (is.character(setseed)) {
    stop("setseed should be TRUE for setseed = 1234567 or FALSE for a random setseed.")
  }
  if (setseed) {
    seed <- 1234567 #This corresponds to set.seed(1234567) used in power.tsd if setseed = TRUE
  } else {
    seed <- runif(1, max=1E7)
    set.seed(seed)
  }
  
  max_iter <- FALSE
  iter <- -1
  current_alpha <- matrix(c(alpha1, alpha2), nrow = 1, ncol = 2)
  # --------------------------------------------------------------------------
  # burn in: make an evaluation with 30 000 sims and choose the points
  # with T1E > 90% percentile
  # dataframe with the entire grid
  d <- cbind(expand.grid(n1 = n1, CV = CV, GMR = GMR, alpha1 = alpha1, alpha2 = alpha2), min.n2 = min.n2)
  res_T1E <- potvin(type = type, d = d, Nmax, targetpower, setseed, nsims = 30000, pmethod = "nct",
                    theta0 = theta2, theta1, theta2, ...)
  T1E  <- res_T1E["pBE",] #This gives the T1E
  m.n2 <- res_T1E["min.n2",] #This gives min.n2
  T1E  <- t(matrix(as.numeric(T1E), nrow = length(CV), ncol = length(n1), byrow = TRUE))
  m.n2 <- t(matrix(as.numeric(m.n2), nrow = length(CV), ncol = length(n1), byrow = TRUE))
  rownames(T1E) <- n1
  colnames(T1E) <- CV
  # now entries with T1E >= 90% percentile
  index <- which(T1E >= quantile(T1E, 0.90), arr.ind = TRUE)
  n1.quant <- as.numeric(rownames(T1E)[index[, 1]])
  m.n2 <- m.n2[index]
  CV.quant <- as.numeric(colnames(T1E)[index[, 2]])
  d90 <- data.frame(n1 = n1.quant,
                    min.n2 = m.n2,
                    CV = CV.quant,
                    GMR = GMR,
                    alpha1 = alpha1, alpha2 = alpha2)

  # LOOP until current_alpha and alpha_adj are different
  # --------------------------------------------------------------------------
  #browser()
  repeat{
    iter <- iter + 1
    # evaluate the T1E grid with 1E6 sims at theta0=theta2 with current_alpha
    d90$alpha1 <- current_alpha[1,1]
    d90$alpha2 <- current_alpha[1,2]
    res_T1E90 <- potvin(type = type, d = d90, Nmax, targetpower, setseed, nsims = 1E+06, pmethod = "nct",
                        theta0 = theta2, theta1, theta2, ...)
    T1E_high <- data.frame(cbind(d90[, "n1"], d90[, "CV"], d90[, "GMR"], d90[, "min.n2"],
                                 d90[, "alpha1"], d90[, "alpha2"],
                                 as.numeric(res_T1E90["pBE", ])))
    colnames(T1E_high) <- c("n1", "CV", "GMR", "min.n2", "alpha1", "alpha2", "pBE")
    # --------------------------------------------------------------------------
    # choose the max T1E
    max_T1E <- T1E_high[T1E_high[, "pBE"] == max(T1E_high["pBE"]), ]
    m <- which.max(max_T1E$CV)
    max_T1E <- max_T1E[m,]
    nperc90 <- cbind(n1 = res_T1E90["n1", ], CV = res_T1E90["CV", ], GMR = res_T1E90["GMR", ],
                     min.n2 = res_T1E90["min.n2", ], pBE_s1 = res_T1E90["pBE_s1", ],
                     pct_s2 = res_T1E90["pct_s2", ], Nperc = res_T1E90["nperc", ])
    max_T1E <- merge(max_T1E, nperc90, all.x = TRUE, by = c("n1", "CV", "GMR", "min.n2"))
    if (details) {
      if (iter < 1) {
        cat("Run-in with alphas", paste(as.numeric(current_alpha), collapse=", "), "\n")
      } else {
        cat("Iteration", iter, "with alphas", paste(as.numeric(floor(current_alpha*1e5)/1e5), collapse=", "), "\n")
      }
        cat("- max.TIE", max_T1E$pBE, "at n1 =", max_T1E$n1, "and loc.CV =" , max_T1E$CV, "\n")
    }
    # --------------------------------------------------------------------------
    # make a grid around the old_alphas
    new_alpha1 <- seq(current_alpha[1, 1] - 0.0005, current_alpha[1, 1] + 0.0005, length.out = 7)
    new_alpha2 <- seq(current_alpha[1, 2] - 0.0005, current_alpha[1, 2] + 0.0005, length.out = 7)
    step_size <- median(diff(seq(current_alpha[1, 2] - 0.0005, current_alpha[1, 2] + 0.0005, length.out = 7)))
     
    if (alpha1 == alpha2) {
      n_d <- cbind(n1 = max_T1E["n1"], CV = max_T1E["CV"], GMR = max_T1E["GMR"],
                   min.n2 = max_T1E["min.n2"], alpha1 = new_alpha1,
                   alpha2 = new_alpha2, row.names = NULL)
    } else {
      n_d <- cbind(n1 = max_T1E["n1"], CV = max_T1E["CV"], GMR = max_T1E["GMR"],
                   min.n2 = max_T1E["min.n2"],
                   expand.grid(alpha1 = alpha1, alpha2 = new_alpha2), row.names = NULL)
    }
    colnames(n_d) <- c("n1", "CV", "GMR", "min.n2", "alpha1", "alpha2")

    # evaluate and use inverse regression to obtain new alphas (alpha.adj)
    res_new_T1E <- potvin(type = type, d = n_d, Nmax, targetpower, setseed, nsims = 1E+06, pmethod = "nct",
                          theta0 = theta2, theta1, theta2, ...) #Takes much time if alpha1 <> alpha2
    res_new_d_T1E <- data.frame(CV = unlist(res_new_T1E["CV", ]),
                                n1 = unlist(res_new_T1E["n1", ]),
                                GMR = unlist(res_new_T1E["GMR", ]),
                                min.n2 = unlist(res_new_T1E["min.n2", ]),
                                alpha1 = sapply(1:ncol(res_new_T1E), function(x) res_new_T1E[["alpha", x]][1]),
                                alpha2 = sapply(1:ncol(res_new_T1E), function(x) res_new_T1E[["alpha", x]][2]),
                                T1E = unlist(res_new_T1E["pBE", ]))

    #INVERSE REGRESSION
    alpha_adj <- inv.reg(alpha, alpha1, alpha2, res_new_d_T1E)
    colnames(alpha_adj) <- c("alpha1", "alpha2")
    #browser()
    
    # Emergency brake
    if(iter > 10) {
      max_iter <- TRUE
      warning("Max. iterations reached.")
      current_alpha <- floor(current_alpha*1e4)/1e4
      repeat{
        max_T1E$alpha1 <- current_alpha[1,1]
        max_T1E$alpha2 <- current_alpha[1,2]
        check_T1E <- potvin(type = type, d = max_T1E, Nmax, targetpower, setseed, nsims = 1E+06, pmethod = "nct",
                            theta0 = theta2, theta1, theta2, ...) #Takes much time if alpha1 <> alpha2
        #Punctual T1E for CV +/- 0.05 is evaluated following the article of Xu et al.  
        lo.power.tsd <- power.tsd(alpha=c(max_T1E$alpha1, max_T1E$alpha2), CV = CV[1] - 0.05, n1 = max_T1E$n1, 
                                  GMR = GMR, min.n2 = max_T1E$min.n2, Nmax = Nmax, targetpower = targetpower, setseed = setseed, 
                                  pmethod = "nct", theta0 = 1.25, theta1 = theta1, theta2 = theta2, method = type) 
        hi.power.tsd <- power.tsd(alpha=c(max_T1E$alpha1, max_T1E$alpha2), CV = tail(CV, n=1) + 0.05, n1 = max_T1E$n1, 
                                  GMR = GMR, min.n2 = max_T1E$min.n2, Nmax = Nmax, targetpower = targetpower, setseed = setseed, 
                                  pmethod = "nct", theta0 = 1.25, theta1 = theta1, theta2 = theta2, method = type)
        if (unlist(check_T1E["pBE", ]) < alpha && lo.power.tsd$pBE < alpha && hi.power.tsd$pBE < alpha) { 
            max_T1E$pBE <- check_T1E["pBE",]
            break 
        } else if (alpha_adj[1] == alpha_adj[2]) { 
          current_alpha <- current_alpha - 1E-4
        } else {
          current_alpha[1] <- floor(alpha1*1e4)/1e4 #In case of alpha 1 different than alpha 2, alph1 is fixed
          current_alpha[2] <- current_alpha[2] - 1E-4
        }
      } #end repeat
      break
    } #end if

    # check that the significance level is below alpha, and the new alphas against the old_one
    #current_alpha[2] because this is the only one we adjust. If alpha1=alpha2, then the alpha2 adjusted
    #is copied in alpha1. If alpha1 and alpha are different, the only that we adjust is alpha2
    diff <- abs(current_alpha[2] - alpha_adj[2])
    #SCENARIO 1/4: T1E < 0.05 and Diff <= 2E-4
    if (max_T1E["pBE"] < alpha & diff <= 2E-4) { 
      current_alpha <- floor(current_alpha*1e4)/1e4
      #The following repeat serves to find a solution when the algorithm enters into a cycle 
      repeat{
        max_T1E$alpha1 <- current_alpha[1,1]
        max_T1E$alpha2 <- current_alpha[1,2]
        check_T1E <- potvin(type = type, d = max_T1E, Nmax, targetpower, setseed, nsims = 1E+06, pmethod = "nct",
                            theta0 = theta2, theta1, theta2, ...) #Takes much time if alpha1 <> alpha2
        #Punctual T1E for CV +/- 0.05 is evaluated following the article of Xu et al.  
        lo.power.tsd <- power.tsd(alpha=c(max_T1E$alpha1, max_T1E$alpha2), CV = CV[1] - 0.05, n1 = max_T1E$n1, 
                                  GMR = GMR, min.n2 = max_T1E$min.n2, Nmax = Nmax, targetpower = targetpower, setseed = setseed, 
                                  pmethod = "nct", theta0 = 1.25, theta1 = theta1, theta2 = theta2, method = type) 
        hi.power.tsd <- power.tsd(alpha=c(max_T1E$alpha1, max_T1E$alpha2), CV = tail(CV, n=1) + 0.05, n1 = max_T1E$n1, 
                                  GMR = GMR, min.n2 = max_T1E$min.n2, Nmax = Nmax, targetpower = targetpower, setseed = setseed, 
                                  pmethod = "nct", theta0 = 1.25, theta1 = theta1, theta2 = theta2, method = type)
        if (unlist(check_T1E["pBE", ]) < alpha && lo.power.tsd$pBE < alpha && hi.power.tsd$pBE < alpha) { 
          max_T1E$pBE <- check_T1E["pBE",]
          break 
        } else if (alpha_adj[1] == alpha_adj[2]) { 
          current_alpha <- current_alpha - 1E-4
        } else {
          current_alpha[1] <- floor(alpha1*1e4)/1e4 #In case of alpha 1 different than alpha 2, alph1 is fixed
          current_alpha[2] <- current_alpha[2] - 1E-4
        }
      } #end repeat
    break
    #SCENARIO 2/4: T1E > 0.05 and Diff <= 2E-4
    } else if (max_T1E["pBE"] > alpha & diff <= 2E-4) { #Scenario 2/4: T1E > 0.05 and Diff <= 1E-4
      if (alpha_adj[2] > current_alpha[2]) {
        if (alpha_adj[1] == alpha_adj[2]) {
          current_alpha <- current_alpha - step_size*(iter+1)
        } else {
          current_alpha[1] <- alpha1 #In case of alpha 1 different than alpha 2, alph1 is fixed
          current_alpha[2] <- current_alpha[2] - step_size*(iter+1)  
        } 
      } else {
        if (alpha_adj[1] == alpha_adj[2]) {
          current_alpha <- alpha_adj
        } else {
          current_alpha[1] <- alpha1 #In case of alpha 1 different than alpha 2, alph1 is fixed
          current_alpha[2] <- alpha_adj[2]  
        }
      }
    #SCENARIO 3/4: Diff > 2E-04 and T1E < 0.05  
    } else if (max_T1E["pBE"] < alpha & diff > 2E-4) { 
      if (max(current_alpha[2], alpha_adj[2]) > alpha) { # Scenario 3 (first option): T1E < 0.05 but either current_alpha or alpha_adj is > alpha
        if (alpha_adj[1] == alpha_adj[2]) {
          current_alpha <- matrix(c(min(current_alpha - step_size*(iter+1), alpha_adj), min(current_alpha - step_size*(iter+1), alpha_adj)), nrow = 1, ncol =2)
          if (current_alpha[2] < 0.025) current_alpha <- matrix(c(alpha/2 + step_size*(iter+1), alpha2/2+ step_size*(iter+1)), nrow = 1, ncol = 2)
        } else {
          current_alpha[1] <- alpha1
          current_alpha[2] <- min(current_alpha[2] - step_size*(iter+1), alpha_adj[2])
          if (current_alpha[2] < 0.025) current_alpha[2] <- alpha/2 + step_size*(iter+1)
        }         
      } else { # Secenario 3 (second option): T1E < 0.05 and current_alpha < alpha
        if (alpha_adj[1] == alpha_adj[2]) {
            current_alpha <- matrix(c(max(current_alpha + step_size*(iter+1), alpha_adj), max(current_alpha + step_size*(iter+1), alpha_adj)), nrow = 1, ncol =2)
        } else {
            current_alpha[1] <- alpha1 #In case of alpha 1 different than alpha 2, alph1 is fixed
            current_alpha[2] <- max(current_alpha[2] + step_size*(iter+1), alpha_adj[2])  
        }
      }
    #SCENARIO 4/4: Diff > 2E-04 and T1E > 0.05  
    } else { 
      if (min(current_alpha[2], alpha_adj[2]) > alpha) { #Both are > alpha
        if (alpha_adj[1] == alpha_adj[2]) {
          current_alpha <- matrix(c(alpha1 - step_size*(iter+1), alpha2 - step_size*(iter+1)), nrow = 1, ncol =2)
        } else {
          current_alpha[1] <- alpha1 #In case of alpha 1 different than alpha 2, alph1 is fixed
          current_alpha[2] <- alpha2 - step_size*(iter+1)  
        } 
      } else {
        if (alpha_adj[1] == alpha_adj[2]) {
          current_alpha <- matrix(c(min(current_alpha - step_size*(iter+1), alpha_adj), min(current_alpha - step_size*(iter+1), alpha_adj)), nrow = 1, ncol =2)
          if (current_alpha[2] < 0.025) current_alpha <- matrix(c(alpha/2 + step_size*(iter+1), alpha2/2+ step_size*(iter+1)), nrow = 1, ncol = 2)
        } else {
          current_alpha[1] <- alpha1 #In case of alpha 1 different than alpha 2, alph1 is fixed
          current_alpha[2] <- min(current_alpha[2] - step_size*(iter+1), alpha_adj[2])
          if (current_alpha[2] < 0.025) current_alpha[2] <- alpha/2 + step_size*(iter+1)
        }
      }  
    }
  } # end repeat loop current_alpha different from alpha_adj
  # ---------------------------------------------------------------------------
  # output
  mt <- if (type == "B") "1 (B)" else "2 (C/D)"
  #browser()
  p <- potvin(type = type, d = max_T1E, Nmax, targetpower, setseed, nsims = 1E+05, pmethod = "nct",
              theta0 = max_T1E[["GMR"]], theta1 = theta1, theta2 = theta2, ...)
  if (print) {
    cat("\nMethod type =", mt, "\n",
        "setseed =", setseed, "\n",
        "BE acceptance range (theta1, theta2) =",
        paste(c(round(unlist(res_T1E["theta1", ][1]), 4),
              round(unlist(res_T1E["theta2", ][1]), 4)), collapse=" ... "), "\n",
        "Power calculation method =", unlist(p["pmethod",][1]), "\n",
        "n1 =", n1, "\n",
        "GMR =", GMR, "\n",
        "CV =", paste(CV, collapse=", "), "\n",
        "Nmax =", Nmax, "\n",
        "min.n2 =", min.n2, "\n",
        "Adjusted alpha at stage 1 =",
        (max_T1E$alpha1*1e4)/1e4, "and alpha at stage 2 =",
        (max_T1E$alpha2*1e4)/1e4, "\n",
        "Maximum empirical type I error =", max_T1E$pBE[[1]], "\n",
        "  at loc.n1 =", max_T1E$n1, "and loc.CV =" , max_T1E$CV, "\n",
        "Power: Overall probability of ABE =", unlist(p["pBE", ]),"\n",
        "Power: Probability of ABE at stage 1 =", unlist(p["pBE_s1", ]),"\n",
        "Studies in stage 2 =", paste0(unlist(p["pct_s2", ]),"%\n"),
        names(p["nperc",][[1]]), "percentiles of N =",
        paste(unlist(p["nperc", ]), collapse=", "),"\n",
        "max.iter =", max_iter, "\n\n"
    )
  }
  if (print == FALSE) {
    pct <- as.numeric(unlist(p["nperc", ], use.names=FALSE))
    res <- list(type = mt, alpha = alpha, CV = CV, n1 = n1, GMR = GMR,
                min.n2 = min.n2, Nmax = Nmax, targetpower = targetpower,
                alpha1 = max_T1E$alpha1, alpha2 = max_T1E$alpha2,
                theta1 = round(unlist(res_T1E["theta1", ][1]), 4),
                theta2 = round(unlist(res_T1E["theta2", ][1]), 4),
                pmethod = unlist(p["pmethod", ]), 
                TIE = max_T1E$pBE[1],
                loc.CV = max_T1E$CV,loc.n1 = max_T1E$n1,
                pBE = as.numeric(unlist(p["pBE", ])),
                pBE_s1 = as.numeric(unlist(p["pBE_s1", ])),
                pct_s2 = as.numeric(unlist(p["pct_s2", ])),
                nperc5 = pct[1], median = pct[2], nperc95 = pct[3],
                max.iter = max_iter)
    return(res)
  }
} # end of t1e.tsd()

# Power.tsd function
# d is a dataframe with alpha1, alpha2, n1, GMR, CV, min.n2
potvin <- function(type, d, Nmax, targetpower, setseed, nsims, pmethod = "nct", theta0=theta0,
                   theta1=theta1, theta2=theta2, npct = c(0.05, 0.5, 0.95)) { # uses FUNCTION power.tsd
  # BY LABES D., LANG B., SCHUETZ H.
  sapply(1:nrow(d), function(x) {
    return(power.tsd(method = type,
                     alpha = c(d[x, "alpha1"], d[x, "alpha2"]),
                     n1  = d[x, "n1"],
                     GMR = d[x, "GMR"],
                     CV  = d[x, "CV"],
                     targetpower = targetpower,
                     pmethod = pmethod,
                     Nmax = Nmax,
                     min.n2 = d[x, "min.n2"],
                     # theta0: True unknown GMR (t_effect);
                     # theta0=theta2 for T1E; theta0=GMR for power
                     theta0 = theta0,
                     theta1 = theta1,
                     theta2 = theta2,
                     npct = npct,
                     nsims = nsims,
                     setseed = setseed))
  }
  )
}

#LOOKING FOR A ALPHA AT EACH STAGE WHOSE T1E < alpha (BASED ON MIN(AIC)) USING INVERSE REGRESSION
inv.reg <- function(alpha, alpha1, alpha2, res_new_d_T1E) {
  mod1 <- lm(T1E ~ alpha2, data = res_new_d_T1E)   # linear
  mod2 <- lm(T1E ~ alpha2 + I(alpha2^2), data = res_new_d_T1E)  #quadratic
  if (extractAIC(mod1, k=2)[2] <= extractAIC(mod2, k=2)[2]) { # select the better model
    #cat("B0 =", coef(mod1)[[1]], "B1 =", coef(mod1)[[2]], "\n")
    alpha.adj <- (alpha-coef(mod1)[[1]])/coef(mod1)[[2]]
  } else {
    #cat("B0 =", coef(mod2)[[1]], "B1 =", coef(mod2)[[2]], "B2 =", coef(mod2)[[3]], "\n")
    det <- (coef(mod2)[[2]]/2/coef(mod2)[[3]])^2-(coef(mod2)[[1]]-alpha)/coef(mod2)[[3]]
    #cat("det =", det, "sqrt.det =", sqrt(det), "\n")
    if (det > 0) {
      if (coef(mod2)[[3]] < 0) {
        alpha.adj <- -(coef(mod2)[[2]]/2/coef(mod2)[[3]]+sqrt(det))
      } else {
        alpha.adj <- -(coef(mod2)[[2]]/2/coef(mod2)[[3]]-sqrt(det))
      }
    } else {
      alpha.adj <- (alpha-coef(mod1)[[1]])/coef(mod1)[[2]]
    }
  } #end else
  #if (alpha.adj > alpha) alpha.adj <- alpha/2
  if (alpha.adj < 0) alpha.adj <- alpha/2 
  if (alpha1 == alpha2) {
    return(matrix(c(alpha.adj, alpha.adj), nrow = 1, ncol = 2))
  } else {
    return(matrix(c(alpha1, alpha.adj), nrow = 1, ncol = 2))
  }
}
