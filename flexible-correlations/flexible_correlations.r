# Dace Apšvalka, 2022, www.dcdace.net
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
  "mvnormalTest", # for multivariate-normality test
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
packages_to_load <- c("ggplot2", "Hmisc", "aplpack", "mvnormalTest", "WRS2", "WRS", "MASS")
# Load the packages:
invisible(lapply(packages_to_load, library, character.only = TRUE))

# =======================================================
# FUNCTION GET OUTLIERS (from a variable pair, for correlations)
# =======================================================
get_outliers <-
  function(var1, var2, # required
  var1name = "var1", var2name = "var2",
  var1ylab = "", var2ylab = "",
  disp = FALSE) {
    # put the 2 univariate and 1 bivariate plots together
    par(mfrow = c(1, 3))

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
    # Report results
    # ------------------------------------
    if (disp) {
      # var1
      print(sprintf("%s outliers: %d", var1name, length(bpVar1$out)))
      print(oVar1)
      # var2
      print(sprintf("%s outliers: %d", var2name, length(bpVar2$out)))
      print(oVar2)
      # bi-var
      print(sprintf("Bi-variate outliers: %d", length(oBiv)))
      print(oBiv)
    }
    # Return results
    # ------------------------------------
    all <- c(oVar1, oVar2, oBiv)
    univariate <- c(oVar1, oVar2)
    bivariate <- oBiv
    return(list(all, univariate, bivariate))
  }

# =======================================================
# FUNCTION DO CORRELATIONS
# =======================================================
# Depending on the data, 3 types of correlations possible. Following recommendations by Pernet et al.(2013)
do_correlation <- function(var1, var2, outliers = NULL) {
  # if outliers not provided, get them
  if (is.null(outliers)) {
    outliers <- get_outliers(var1, var2)
  }
  # WHICH CORRELATION
  # ------------------------------------
  isOutliers <- length(outliers[[1]]) > 0
  isUnivariate <- length(outliers[[2]]) > 0
  isBivariate <- length(outliers[[3]]) > 0
  # using Henze-Zirkler Test for Multivariate Normality
  isNormal <- mhz(cbind(var1, var2))$mv.test["p-value"] > 0.05

  # if is Bivariate do Spearman skipped, using the minimum covariance determinant (MCD) estimator
  if (isBivariate) {
    corRes <-
      mscor(data.frame(var1, var2), corfun = spear, cop = 3)
    r <- corRes$cor[2]
    p <- 2 * pt(-abs(corRes$test.stat[2]), df = length(var1) - 1)
    f <- "Spearman skipped"
    subs <- "ss"
  }

  # if is not Bivar but is Univar or is not Normal, do 20% Bend
  if (!isBivariate & (isUnivariate | !isNormal)) {
    corRes <-
      pbcor(var1, var2, beta = 0.2)
    r <- corRes$cor
    p <- corRes$p.value
    f <- "Percentage-bend"
    subs <- "pb"
  }

  # if no outliers and is normal do Pearson
  if (!isOutliers & isNormal) {
    corRes <-
      rcorr(var1, var2)
    r <- corRes$r[2]
    p <- corRes$P[2]
    f <- "Pearson"
    subs <- ""
  }

  pval <- ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  corResTxt <- bquote(.(f) ~ "correlation" ~ r[.(subs)] == .(sprintf("%.3f, %s", r, pval)))

  return(list(corResTxt, p))
}

# =======================================================
# FUNCTION PLOT CORRELATIONS (with 95%CI)
# =======================================================
plot_correlation <- function(var1, var2, #required
  var1name = "var1", var2name = "var2", # axis lables
  corRes = NULL,
  pointsize = 1.8, txtsize = 11, # default point and font size
  outliers = NULL,
  plotoutliers = FALSE,
  pthreshold = NULL) {
  # If outliers not given, get them
  if (is.null(outliers)) {
    outliers <- get_outliers(var1, var2)
  }
  # If correlation results not given, get them
  if (is.null(corRes)) {
    corRes <- do_correlation(var1, var2, outliers)
  }
  # Format the output
  if (length(outliers[[1]]) > 0) {
    out1 <- var1[outliers[[1]]]
    out2 <- var2[outliers[[1]]]
    outdata <- data.frame(out1, out2)

    var1 <- var1[-c(outliers[[1]])]
    var2 <- var2[-c(outliers[[1]])]
    ifelse(plotoutliers == TRUE,
           addTxt <- sprintf("\n Outliers (n=%d) displayed in red", length(unique(outliers[[1]]))),
           addTxt <- sprintf("\n Outliers (n=%d) not displayed", length(unique(outliers[[1]])))
           )
    resTXT <- bquote(atop(.(corRes[[1]]), .(addTxt)))
  }
  else {
    resTXT <- corRes[[1]]
  }

  # If alpha not defined, set it to 0 (to ignore it)
  if (is.null(pthreshold)) {
    pthreshold <- 0
  }
  if (corRes[[2]] < pthreshold) {
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
    labs(x = var1name, y = var2name, subtitle = resTXT) +
    theme_minimal() +
    theme(text = element_text(size = txtsize),
          plot.title = element_text(hjust = 0.5, size = txtsize + 2, face = titleface, color = "black"),
          plot.subtitle = element_text(hjust = 0.5),
          plot.background = element_rect(colour = framecolor, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
          )
  # if asked to display, add outliers on the plot
  if (length(outliers[[1]]) > 0 & plotoutliers == TRUE) {
    corplot <- corplot + geom_point(data = outdata, aes(out1, out2), color = "red")
  }
  return(corplot)
}