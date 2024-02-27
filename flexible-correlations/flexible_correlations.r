#
# Dace Apšvalka, 2022, www.dcdace.net
#
# =======================================================
# INSTALL AND LOAD THE REQUIRED PACKAGES
# =======================================================
# Define a function that checks and installs missing packages
install_packages <- function(packages) {
  lapply(packages,
    FUN = function(x)(
        if (length(find.package(x, quiet = TRUE)) == 0) {
          install.packages(x, dependencies = TRUE)
        }))
}
# A list of required packages
packages_required <- c(
  "ggplot2", # for plotting
  "Hmisc", # for several stats things
  "aplpack", # for bi-variate outlier, bgplot
  "mvnTest", # for multivariate-normality test
  "WRS2", # for percentage bend correlation coefficient. See https://link-springer-com.ezp.lib.cam.ac.uk/content/pdf/10.1007/BF02294395.pdf
# WRS for Wilcox correlation functions (skipped correlations)
# WRS requires several dependent packages, as specified here https://github.com/nicebread/WRS
  "MASS", "akima", "robustbase", # MASS's select function conflicts with dplyr select!
  "cobs", "robust", "mgcv", "scatterplot3d", "quantreg", "rrcov", "lars", "pwr", "trimcluster", "mc2d", "psych", "Rfit", "DepthProc", "class", "fda", "rankFD",
  "devtools" # to install WRS from the github site
  )
# Install the missing required packages
install_packages(packages_required)
# Install WRS package from GitHub
if (length(find.package("WRS", quiet = TRUE)) == 0) {
  library(devtools)
  install_github("nicebread/WRS", subdir = "pkg")
}

# Not all required packages need to be loaded. Only load the ones that are needed
# A list of packages that need to be loaded
packages_to_load <- c("ggplot2", "Hmisc", "aplpack", "mvnTest", "WRS2", "WRS", "MASS")
# Load the packages:
invisible(lapply(packages_to_load, library, character.only = TRUE))

# =======================================================
# FUNCTION CHECK DATA (for outliers and normality)
# =======================================================
check_data <- function(var1, var2, # required
  var1name = "var1", var2name = "var2",
  var1ylab = "", var2ylab = "",
  disp = FALSE) {
  # put the 2 univariate, 1 bivariate, and Q-Q plots together
  par(mfrow = c(1, 4))

  # UNIVARIATE, boxplot method
  # ------------------------------------
  bpVar1 <- boxplot(var1, main = var1name, ylab = var1ylab, plot = disp)
  bpVar2 <- boxplot(var2, main = var2name, ylab = var2ylab, plot = disp)
  # oVar1
  ifelse(length(bpVar1$out) == 0,
      oVar1 <- vector(mode = "numeric", length = 0),
      oVar1 <- which(var1 %in% bpVar1$out)
    )
  # oVar2
  ifelse(length(bpVar2$out) == 0,
      oVar2 <- vector(mode = "numeric", length = 0),
      oVar2 <- which(var2 %in% bpVar2$out)
    )

  # BI-VARIATE, bagplot method
  # ------------------------------------
  # A bagplot is a bivariate generalization of the  boxplot. A bivariate boxplot.
  bpl <- compute.bagplot(var1, var2)
  # oBiv
  ifelse(is.null(bpl$pxy.outlier),
      oBiv <- vector(mode = "numeric", length = 0),
      oBiv <- which(var1 %in% bpl$pxy.outlier[, 1])
    )
  bp <-
      bagplot(
        var1, var2,
        ylab = var2name, xlab = var1name,
        main = "bi-variate",
        cex = 2,
        create.plot = disp
      )

  # BI-VARIATE NORMALITY
  # ------------------------------------
  # using Henze-Zirkler Test for Multivariate Normality
  bvn <- HZ.test(cbind(var1, var2), qqplot = disp)
  bvn_result <- ifelse(bvn@p.value > 0.05, "Data are bi-variate normal", "Data are not bi-variate normal")

  # Report results
  # ------------------------------------
  if (disp) {
    # var1
    cat(sprintf("\n%s outlier cases(n=%d): %s \n", var1name, length(bpVar1$out), paste(oVar1, collapse = ",")))
    # var2
    cat(sprintf("%s outlier cases(n=%d): %s \n", var2name, length(bpVar2$out), paste(oVar2, collapse = ",")))
    # bi-var
    cat(sprintf("Bi-variate outlier cases(n=%d): %s \n", length(oBiv), paste(oBiv, collapse = ",")))
    # normality
    cat(sprintf("\nHenze-Zirkler test for Multivariate Normality:\n%s(HZ = %.2f, p = %.2f)\n",
      bvn_result, bvn@HZ, bvn@p.value))
  }

  # Return results
  # ------------------------------------
  data_check <- list()
  data_check$outliers_all <- c(oVar1, oVar2, oBiv)
  data_check$outliers_uni <- c(oVar1, oVar2)
  data_check$outliers_bi <- oBiv
  data_check$normality <- bvn@p.value > 0.05

  return(data_check)
}

# =======================================================
# FUNCTION DO CORRELATIONS
# =======================================================
# Depending on the data, 3 types of correlations possible. Following recommendations by Pernet et al.(2013)
do_correlation <- function(var1, var2, data_check = NULL) {
  corr_results <- list()
  datainfo <- list()
  i <- 0
  # if data check not provided, get them
  if (is.null(data_check)) {
    data_check <- check_data(var1, var2)
  }
  # WHICH CORRELATION
  # ------------------------------------
  isOutliers <- length(data_check$outliers_all) > 0
  isUnivariate <- length(data_check$outliers_uni) > 0
  isBivariate <- length(data_check$outliers_bi) > 0
  isNormal <- data_check$normality

  # if is Bivariate do Spearman skipped, using the minimum covariance determinant (MCD) estimator
  if (isBivariate) {
    corr_results <-
      mscor(data.frame(var1, var2), corfun = spear, cop = 3)
    r <- corr_results$cor[2]
    p <- 2 * pt(-abs(corr_results$test.stat[2]), df = length(var1) - 1)
    f <- "Spearman skipped"
    subs <- "ss"
    i <- i + 1
    datainfo[i] <- "bi-variate outlier/s"
  }

  # if is not Bivar but is Univar or is not Normal, do 20% Bend
  if (!isBivariate & (isUnivariate | !isNormal)) {
    corr_results <-
      pbcor(var1, var2, beta = 0.2)
    r <- corr_results$cor
    p <- corr_results$p.value
    f <- "Percentage-bend"
    subs <- "pb"
  }
  if (isUnivariate) {
    i <- i + 1
    datainfo[i] <- "univariate outlier/s"
  }
  if (!isNormal) {
    i <- i + 1
    datainfo[i] <- "not bi-variate normality"
  }

  # if no outliers and is normal do Pearson
  if (!isOutliers & isNormal) {
    corr_results <-
      rcorr(var1, var2)
    r <- corr_results$r[2]
    p <- corr_results$P[2]
    f <- "Pearson"
    subs <- ""
    i <- i + 1
    datainfo[i] <- "no outliers and has bi-variate normality"
  }

  pval <- ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  corr_results$txt <- bquote(.(f) ~ "correlation" ~ r[.(subs)] == .(sprintf("%.3f, %s", r, pval)))
  corr_results$p <- p
  corr_results$datainfo <- sprintf("Data: %s", paste(datainfo, collapse = "; "))

  return(corr_results)
}

# =======================================================
# FUNCTION PLOT CORRELATIONS (with 95%CI)
# =======================================================
plot_correlation <- function(var1, var2, #required
  var1name = "var1", var2name = "var2", # axis lables
  corr_results = NULL,
  data_check = NULL,
  pointsize = 1.8, txtsize = 11, # default point and font size
  plotoutliers = FALSE,
  pthreshold = NULL,
  datainfo = TRUE) {
  # If data_check not given, get them
  if (is.null(data_check)) {
    data_check <- check_data(var1, var2)
  }
  # If correlation results not given, get them
  if (is.null(corr_results)) {
    corr_results <- do_correlation(var1, var2, data_check)
  }
  # Format the output
  if (length(data_check$outliers_all) > 0) {
    out1 <- var1[data_check$outliers_all]
    out2 <- var2[data_check$outliers_all]
    outdata <- data.frame(out1, out2)

    var1 <- var1[-c(data_check$outliers_all)]
    var2 <- var2[-c(data_check$outliers_all)]
    ifelse(plotoutliers == TRUE,
           addTxt <- sprintf("\n Outliers (n=%d) displayed in red", length(unique(data_check$outliers_all))),
           addTxt <- sprintf("\n Outliers (n=%d) not displayed", length(unique(data_check$outliers_all)))
           )
    corr_results$txt <- bquote(atop(.(corr_results$txt), .(addTxt)))
  }

  # If alpha not defined, set it to 0 (to ignore it)
  if (is.null(pthreshold)) {
    pthreshold <- 0
  }
  if (corr_results$p < pthreshold) {
    titleface <- "bold"
    framecolor <- "red"
  } else {
    titleface <- "plain"
    framecolor <- "lightgrey"
  }

  dataset <- data.frame(var1, var2)

  corplot <- ggplot(dataset, aes(var1, var2)) +
    geom_smooth(
      method = lm, formula = y ~ x, level = 0.95,
      color = "black", fill = "grey", size = 0.3
    ) +
    geom_point(
      colour = "black", alpha = .8, fill = "orange",
      size = pointsize, stroke = 0.2, shape = 21
    ) +
    labs(x = var1name, y = var2name, subtitle = corr_results$txt) +
    theme_minimal() +
    theme(text = element_text(size = txtsize),
          plot.title = element_text(hjust = 0.5, size = txtsize + 2, face = titleface, color = "black"),
          plot.subtitle = element_text(hjust = 0.5),
          plot.background = element_rect(colour = framecolor, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
          )
  # if asked to display, add outliers on the plot
  if (length(data_check$outliers_all) > 0 & plotoutliers == TRUE) {
    corplot <- corplot + geom_point(data = outdata, aes(out1, out2), color = "red")
  }
  # if asked to show data info, add it to the plot caption
  if (datainfo) {
    corplot <- corplot + labs(caption = corr_results$datainfo)
  }
  return(corplot)
}