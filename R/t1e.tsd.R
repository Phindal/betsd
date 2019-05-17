# ------------------------------------------------------------------------------
# funnction to Significance level adjustment for bioequivalence studies using 
# two-stage adaptive 2x2 crossover designs.
#
# Authors: E. Molins, D. Labes, H. Schuetz, J. Ocaña
# ------------------------------------------------------------------------------
#PREDICTS ALPHA1 AND ALPHA2 ON A TSD 2x2 CXO DESIGN, CONTROLLING T1E BELOW ALPHA
t1e.tsd <- function(n1, CV, GMR = 0.95, Nmax = 150, min.n2 = n1/2, type = 1,
                    alpha = 0.05, alpha1 = 0.0294, alpha2 = alpha1,
                    targetpower = 0.8, setseed = TRUE, details = TRUE, print = TRUE,
                    ...) {

  #Non-modified parameters => Nmax == Inf & min.n2 == 0
  #Modified CRITERIA ACCORDING TO DOI: 10.1002/sim.7452
  min.n2 <- sapply(min.n2, function(y) if (y %% 2 != 0) y+y%%2 else y)

  #DEBUG AT FIRST ITERATION
  if (type!=1 & type!=2)
    stop("type must be 1 or 2")
  if (type == 1) {type <- "B"} else {type <- "C"} #type 1: Potvin "B"; type 2: Potvin "C or D"
  if (missing(n1))
    stop("Number of subjects in stage 1 must be given.")
  if (any(n1 < 12))
    stop("Number of subjects in stage 1 must be at least 12.")
  if (missing(CV))
    stop("CV must be given.")
  if (GMR <= 0.8 | GMR >= 1.25)
    stop("GMR must be within acceptance range.")
  if (any(CV < 0))
    stop("CV must be > 0.")
  if (length(type) > 1)
    stop("type must be of length 1.")

  #SETSEED
  if (is.numeric(setseed)) {
    setseed <- TRUE
  } else if (is.character(setseed)) {
    stop("setseed should be TRUE for setseed = 1234567 or FALSE for a random setseed.")
  }
  if (setseed) {
    setseed <- 1234567
  } else {
    setseed <- runif(1, max=1E7)
  }

  iter <- -1
  old_alpha <- matrix(c(alpha1, alpha2), nrow = 1, ncol = 2)
  # --------------------------------------------------------------------------
  # burn in: make an evaluation with 30 000 sims and choose the points
  # with T1E > 90% percentile

  # dataframe with the entire grid
  d <- cbind(expand.grid(n1 = n1, CV = CV, GMR = GMR, alpha1 = alpha1, alpha2 = alpha2), min.n2 = min.n2)
  res_T1E <- potvin(type = type, d = d, Nmax, targetpower, setseed, nsims = 30000, ...)
  T1E  <- res_T1E["pBE",] #This gives the T1E
  m.n2 <- res_T1E["min.n2",] #This gives min.n2
  T1E  <- t(matrix(as.numeric(T1E), nrow = length(CV), ncol = length(n1), byrow = TRUE))
  m.n2 <- t(matrix(as.numeric(m.n2), nrow = length(CV), ncol = length(n1), byrow = TRUE))
  rownames(T1E) <- n1
  colnames(T1E) <- CV
  # now entries with T1E >= 90% percentile
  index <- which(T1E >= quantile(T1E, 0.90), arr.ind = TRUE)
  n1.quant <- as.numeric(rownames(T1E)[index[,1]])
  m.n2 <- m.n2[index]
  CV.quant <- as.numeric(colnames(T1E)[index[,2]])
  d90 <- data.frame(n1 = n1.quant,
                    min.n2 = m.n2,
                    CV = CV.quant,
                    GMR = GMR,
                    alpha1 = alpha1, alpha2 = alpha2)

  # LOOP until old_alpha and alpha_adj are different
  # --------------------------------------------------------------------------
  #browser()
  repeat{
    iter <- iter + 1
    # evaluate the T1E grid with 1E6 sims at theta0=1.25 with old_alpha
    d90$alpha1 <- old_alpha[1,1]
    d90$alpha2 <- old_alpha[1,2]
    res_T1E90 <- potvin(type = type, d = d90, Nmax, targetpower, setseed,
                        nsims = 1E+06, ...)
    T1E_high <- data.frame(cbind(d90[,"n1"], d90[,"CV"], d90[,"GMR"], d90[,"min.n2"],
                                 d90[,"alpha1"], d90[,"alpha2"],
                                 as.numeric(res_T1E90["pBE",])))
    colnames(T1E_high) <- c("n1", "CV", "GMR", "min.n2", "alpha1", "alpha2", "pBE")
    # --------------------------------------------------------------------------
    # choose the max T1E
    max_T1E <- T1E_high[T1E_high[,"pBE"] == max(T1E_high["pBE"]), ]
    m <- which.max(max_T1E$CV)
    max_T1E <- max_T1E[m,]
    nperc90 <- cbind(n1 = res_T1E90["n1",], CV = res_T1E90["CV",], GMR = res_T1E90["GMR",],
                     min.n2 = res_T1E90["min.n2",], pBE_s1 = res_T1E90["pBE_s1",],
                     pct_s2 = res_T1E90["pct_s2",], Nperc = res_T1E90["nperc",])
    max_T1E <- merge(max_T1E, nperc90, all.x = TRUE, by = c("n1", "CV", "GMR", "min.n2"))
    if (details) {
      if (iter < 1) {
        cat("Run-in with alphas", paste(as.numeric(old_alpha), collapse=", "), "\n")
      } else {
        cat("Iteration", iter, "with alphas", paste(as.numeric(floor(old_alpha*1e5)/1e5), collapse=", "), "\n")
      }
      cat("- max.TIE", max_T1E$pBE, "at n1 =", max_T1E$n1, "and CV =" , max_T1E$CV, "\n")
    }
    # --------------------------------------------------------------------------
    # make a grid around the old_alphas
    new_alpha1 <- seq(old_alpha[1,1] - 0.0005, old_alpha[1,1] + 0.0005, length.out = 7)
    new_alpha2 <- seq(old_alpha[1,1] - 0.0005, old_alpha[1,1] + 0.0005, length.out = 7)
    if (alpha1 == alpha2) {
      n_d <- cbind(n1 = max_T1E["n1"], CV = max_T1E["CV"], GMR = max_T1E["GMR"],
                   min.n2 = max_T1E["min.n2"], alpha1 = new_alpha1,
                   alpha2 = new_alpha2, row.names = NULL)
    } else {
      n_d <- cbind(n1 = max_T1E["n1"], CV = max_T1E["CV"], GMR = max_T1E["GMR"],
                   min.n2 = max_T1E["min.n2"],
                   expand.grid(alpha1 = new_alpha1, alpha2 = new_alpha2), row.names = NULL)
    }
    colnames(n_d) <- c("n1", "CV", "GMR", "min.n2", "alpha1", "alpha2")

    # evaluate and use inverse regression to obtain new alphas (alpha.adj)
    res_new_T1E <- potvin(type = type, d = n_d, Nmax, targetpower, setseed, nsims = 1E+06, ...) #Takes much time if alpha1 <> alpha2
    res_new_d_T1E <- data.frame(CV = unlist(res_new_T1E["CV",]),
                                n1 = unlist(res_new_T1E["n1",]),
                                GMR = unlist(res_new_T1E["GMR",]),
                                min.n2 = unlist(res_new_T1E["min.n2",]),
                                alpha1 = sapply(1:ncol(res_new_T1E), function(x) res_new_T1E[["alpha",x]][1]),
                                alpha2 = sapply(1:ncol(res_new_T1E), function(x) res_new_T1E[["alpha",x]][2]),
                                T1E = unlist(res_new_T1E["pBE",]))

    #INVERSE REGRESSION
    alpha_adj <- inv.reg(alpha, alpha1, alpha2, res_new_d_T1E)
    colnames(alpha_adj) <- c("alpha1","alpha2")

    # check that the significance level is below alpha, and the new alphas against the old_one
    diff <- sum(abs(old_alpha - alpha_adj))
    if (max_T1E["pBE"] <= alpha & diff <= 1E-4) break

    # Emergency brake
    if(iter > 10) {
      warning("Max. iterations reached.")
      break
    }

    # start over with new alphas
    if (all(old_alpha == alpha_adj)) {
      old_alpha <- old_alpha - 0.0001
    } else {
      old_alpha <- alpha_adj
    }

  } # end repeat loop old_alpha different from alpha_adj
  # ---------------------------------------------------------------------------
  # output of whatever comes here
  # and/or a return structure as value
  mt <- if (type == "B") "1 (B)" else "2 (C/D)"
  #browser()
  p_d <- cbind(n1 = max_T1E["n1"], CV = max_T1E["CV"], GMR = max_T1E["GMR"],
               min.n2 = max_T1E["min.n2"], alpha1 = alpha_adj[1],
               alpha2 = alpha_adj[2], row.names = NULL)
  p <- potvin(type = type, d = p_d, Nmax, targetpower, setseed, nsims = 1E+05, theta0 = 0.95, ...)
  if (print) {
    cat("\nMethod type =", mt, "\n",
      "setseed =", setseed, "\n",
      "LBL , UBL =", unlist(res_T1E["theta1",][1]),",", unlist(res_T1E["theta2",][1]), "\n",
      "Power calculation method:", unlist(p["pmethod",][1]), "\n",
      "n1 =", n1, "\n",
      "GMR =", GMR, "\n",
      "CV =", CV, "\n",
      "Nmax =", Nmax, "\n",
      "min.n2 =", min.n2, "\n",
      "Adjusted alpha at stage 1 =", floor(max_T1E$alpha1*1e4)/1e4, "and alpha at stage 2 =", floor(max_T1E$alpha2*1e4)/1e4, "\n",
      "Maximum empirical type I error =", max_T1E$pBE[1], "\n",
      "  at n1=", max_T1E$n1, "and CV=" , max_T1E$CV, "\n",
      "Power: Overall probability of ABE =", unlist(p["pBE",]),"\n",
      "Power: Probability of ABE at stage 1 =", unlist(p["pBE_s1",]),"\n",
      "Studies (%) in stage 2 =", unlist(p["pct_s2",]),"\n",
      "Percentiles",names(p["nperc",][[1]]),"of N =", unlist(p["nperc",]),"\n\n"
      )
  }
  if (print == FALSE) {
    pct <- as.numeric(unlist(p["nperc",], use.names=FALSE))
    res <- list(type=mt,alpha=alpha,CV=CV,n1=n1,GMR=GMR,min.n2=min.n2,Nmax=Nmax,
                targetpower=targetpower,
                alpha1=max_T1E$alpha1,alpha2=max_T1E$alpha2,
                TIE=max_T1E$pBE[1],
                loc.CV=max_T1E$CV,loc.n1=max_T1E$n1,
                pBE=as.numeric(unlist(p["pBE",])),
                pBE_s1=as.numeric(unlist(p["pBE_s1",])),
                pct_s2=as.numeric(unlist(p["pct_s2",])),
                nperc5=pct[1],median=pct[2],nperc95=pct[3])
    return(res)
  }
} # end of t1e.tsd()

# Power.tsd function
# d is a dataframe with alpha1, alpha2, n1, GMR, CV, min.n2
potvin <- function(type, d, Nmax, targetpower, setseed, nsims, pmethod = "nct", theta0 = 1.25, theta1 = 0.8, theta2 = 1.25, npct = c(0.05, 0.5, 0.95)) { # uses FUNCTION power.tsd
    # BY LABES D., LANG B., SCHUETZ H.
    sapply(1:nrow(d), function(x) {
      return(power.tsd(method = type,
                       alpha = c(d[x,"alpha1"], d[x,"alpha2"]),
                       n1  = d[x,"n1"],
                       GMR = d[x, "GMR"],
                       CV  = d[x,"CV"],
                       targetpower = targetpower,
                       pmethod = pmethod,
                       Nmax = Nmax,
                       min.n2 = d[x, "min.n2"],
                       # theta0: True unknown GMR (t_effect);
                       # theta0=1.25 for T1E; theta0=0.95 for power
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
    if (det < 0) {
      det <- 0
    }
    if (coef(mod2)[[3]] < 0) {
      alpha.adj <- -(coef(mod2)[[2]]/2/coef(mod2)[[3]]+sqrt(det))
    } else {
      alpha.adj <- -(coef(mod2)[[2]]/2/coef(mod2)[[3]]-sqrt(det))
    }
  }
  #To avoid an alpha.adj < 0 (and the resulting error when running the following example:
  #t1e.tsd(n1 = 12, CV = 0.6, GMR = 0.9, type = 2), the next two rows are added:
  if (alpha.adj > alpha) alpha.adj <- alpha #Rarely happens, maybe with CV = 0.6
  if (alpha.adj < 0) alpha.adj <- 0.025 #Rarely happens, maybe with CV = 0.6
  if (alpha1 == alpha2) {
    return(matrix(c(alpha.adj, alpha.adj), nrow = 1, ncol = 2))
  } else {
    return(matrix(c(alpha1, alpha.adj), nrow = 1, ncol = 2))
  }
}
