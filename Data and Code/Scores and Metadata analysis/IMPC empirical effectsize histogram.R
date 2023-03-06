library(DRrequiredAgeing)
df = read.csv(tail(
  UnzipAndfilePath(
    'WithQvaluesAndMPtermsdata_point_of_type_unidimensional.zip'
  ),
  1
), stringsAsFactors = TRUE)
hist(
  abs(df$N_Genotype.standardised.effect.size),
  breaks = 20,
  xlab = 'Standardised Effect Size (SES)',
  main = '',
  probability = TRUE,
  ylim = c(0, 1),
  cex.lab = 2
)
mses = mean(abs(df$N_Genotype.standardised.effect.size), trim = .1)
abline(v = mses,
       lwd = 5,
       col = 2,
       lty = 2)
abline(v = .5,
       lwd = 5,
       col = 3,
       lty = 2)
legend(
  'topright',
  legend = c(
    paste0(
      'Empirical standardised mean effect size from the data (',
      round(mses, 2),
      ')'
    ),
    '0.5'
  ),
  xpd = FALSE,
  horiz = FALSE,
  col = 2:3,
  lwd = 5,
  lty = 2,
  cex = 1.2
)
