load('Rawdata_IMPC_HEM_onlyControls.Rdata')
df = subset(df, df$biological_sample_group == 'control')
# Filter for the date of experiment
df = df[df$date_of_experiment >= '2018-01-01T00:00:00Z' &
          df$date_of_experiment <= '2020-12-31T00:00:00Z',]
df = as.data.frame(df)

centres = sort(unique(df$phenotyping_center))
parameters = unique(df$parameter_name)
df$data_point = as.numeric(df$data_point)
df = df[!is.na(df$data_point),]
df$date_of_experiment = as.Date(df$date_of_experiment)
df = df[order(df$phenotyping_center),]
pdf('VariationOverTime.pdf', width = 15, height = 13)


for (parameter in parameters) {
  df1 = subset(df, df$parameter_name %in% parameter)
  tbldf1 = table(df1$production_center)
  nm = length(tbldf1[tbldf1 > 0])
  if (nm < 1)
    next
  par(mar = c(10, 5, 3, 3),
      mfrow = c(4, ceiling(nm / 4)))
  for (centre in centres) {
    df2 = subset(df,
                 df$parameter_name %in% parameter &
                   df$phenotyping_center %in% centre)
    if (nrow(df2) < 1 || length(unique(df2$date_of_experiment)) < 3)
      next
    
    colm = corrplot::COL2('RdYlBu', length(unique(df2$date_of_experiment)))
    colf = corrplot::COL2('PRGn', length(unique(df2$date_of_experiment)))
    col2 = c()
    for (i in 1:length(colm)) {
      j = 2 * i
      col2[j] = colm[i]
      col2[j - 1] = colm[i]
    }
    
    
    boxplot(
      df2$data_point ~ df2$date_of_experiment,
      ylab = parameter,
      main = paste0('IMPC centre: ', centre),
      las = 3,
      xlab = '',
      col = col2,
      outline = FALSE
    )
    abline(
      h = mean(df2$data_point),
      lty = 1,
      lwd = 3,
      col = 4
    )
  }
}
graphics.off()
