# Various R-scripts

## 2×2 repeated measures ANOVA
[2×2 repeated measures ANOVA R-script](https://github.com/dcdace/R_functions/tree/main/RM-2by2-ANOVA) `rm_2by2_anova.r` outputs assumption checks (outliers and normality), interaction and main effect results, pairwise comparisons, and produces a result plot with within-subject error bars (SD, SE or 95% CI) and significance stars added to the plot.

## Flexible correlation
[Flexible correlations R-script](https://github.com/dcdace/R_functions/tree/main/flexible-correlations) `flexible_correlations.r` first inspects the data for outliers and normality and then chooses the most appropriate from three correlation methods. Based on [Pernet et al. (2013)](https://doi.org/10.3389/fpsyg.2012.00606 "Pernet, C. R., Wilcox, R. R., & Rousselet, G. A. (2013). Robust correlation analyses: false positive and power validation using a new open source matlab toolbox. Frontiers in psychology, 606."), I follow three simple rules for selecting the correlation method:
* **Pearson's correlation**: Data is normally distributed and has no outliers
* **Spearman skipped correlation**: Data has bi-variate outliers
    * Using the minimum covariance determinant (MCD) estimator
* **(20%) Percentage-bend correlation**: Data has no bi-variate outliers but is not normally distributed or has univariate outliers
