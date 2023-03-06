library(DRrequiredAgeing)
# Install the package above via:
# library(devtools)
# install_github(
#   repo = 'mpi2/impc_stats_pipeline/Late adults stats pipeline/DRrequiredAgeing/DRrequiredAgeingPackage',
#   dependencies = FALSE,
#   upgrade = 'always',
#   force = TRUE,
#   build = TRUE,
#   quiet = FALSE
# )
# Load meta analysis R package
library(metafor)

# Load the raw results
df = read.csv(tail(
  UnzipAndfilePath(
    'ForScoresWithQvaluesAndMPtermsdata_point_of_type_unidimensional.zip'
  ),
  1
), stringsAsFactors = TRUE)

# Filter the unnecessary results (for late adult (LA) mice)
# We do not use the LA in this paper
dfearly = droplevels(subset(
  df,!grepl(
    pattern = 'LA_',
    x = df$parameter_stable_id,
    fixed = TRUE
  )
))


# The Driver function
# The Driver function is written is a way that is more general to this paper
# The default set of parameters in this function should reproduce exact results
# in the paper
MainFunction = function(df,
                  pattern = 'HEM',
                  centres = NULL,
                  late = FALSE,
                  minRequiredPerSex = 50,
                  qvalThresh = 0.05,
                  expES = .5,
                  reverse = FALSE) {
  if (!is.null(centres))
    centres = sort(centres)
  df = droplevels(subset(
    df,
    grepl(
      pattern = pattern,
      x = df$parameter_stable_id,
      fixed = TRUE
    )
  ))
  # Remove tests with less than xx animals
  dfraw = df
  df =  subset(df, df$status %in% 'Successful')
  df = df[df$Total.females >= minRequiredPerSex,]
  df = df[df$Total.males >= minRequiredPerSex,]
  df = droplevels(df)
  
  dfraw$parnameid = paste(dfraw$parameter_name, dfraw$parameter_stable_id, sep = '.')
  df$parnameid = paste(df$parameter_name, df$parameter_stable_id, sep = '.')
  
  
  ####################################
  # For the meta analysis
  ####################################
  agESMeta = aggregate(df[, c('Sex.standardised.effect.size',
                              'Sex.estimate',
                              'Sex.standard.error')], by = list(df$parnameid, df$phenotyping_center), function(x) {
                                mean(x)
                              })
  # From the supplementary materials in https://academic.oup.com/ije/article/47/4/1343/5042988?login=false
  for (iparm in seq_along(unique(agESMeta$Group.1))) {
    parm = unique(agESMeta$Group.1)[iparm]
    meta <-
      rma(
        yi = agESMeta$Sex.estimate[agESMeta$Group.1 %in% parm],
        sei = agESMeta$Sex.standard.error[agESMeta$Group.1 %in% parm],
        method = 'REML'
      )
    eff <- unlist(meta[c("b", "ci.lb", "ci.ub", "se", "pval")])
    eff <- t(as.matrix(eff))
    colnames(eff) <-
      c(
        "MetaAnalysis.Estimate",
        "MetaAnalysis.low_95%",
        "MetaAnalysis.upp_95%",
        "MetaAnalysis.S.E.",
        "MetaAnalysis.p-value"
      )
    eff = as.data.frame(eff)
    eff = cbind('parameter' = parm, eff)
    
    if (iparm <= 1) {
      effinal <- eff
    } else{
      effinal = rbind(effinal, eff)
    }
  }
  write.csv(effinal, file = 'MetaAnalysisResults.csv', row.names = FALSE)
  ####################################
  
  
  ####################################
  # For scoring system
  ####################################
  agES = aggregate(df$Sex.standardised.effect.size, by = list(df$parnameid, df$phenotyping_center), function(x) {
    mean(x)
  })
  agQv = aggregate(df$Sex.q.value, by = list(df$parnameid, df$phenotyping_center), function(x) {
    if (reverse)
      r = mean(x) > qvalThresh
    else
      r =  mean(x) < qvalThresh
    return (r)
  })
  agQv2 = aggregate(df$Sex.q.value, by = list(df$parnameid, df$phenotyping_center), function(x) {
    mean(x)
  })
  agES$qval = agQv2$x
  agES$hit = agQv$x
  names(agES) = c('parameter',
                  'centre',
                  'sex.standardised.effect.size',
                  'qvalue',
                  'is_hit')
  output = NULL
  counter = 1
  # Prepare to calculate scores
  for (par in unique(agES$parameter)) {
    r = NULL
    score = 0
    ESs = c()
    testcounter = 0
    for (cen in centres) {
      v = subset(agES, agES$parameter == par & agES$centre == cen)
      v2 = subset(dfraw,
                  dfraw$parnameid == par &
                    dfraw$phenotyping_center == cen)
      
      
      if (nrow(v) > 0) {
        if (v$is_hit) {
          v100 = v$sex.standardised.effect.size
        } else{
          v100 = ifelse(reverse,
                        '(reverse) Not_significant',
                        'Not_significant')
        }
        score = score + sqrt(abs(v$sex.standardised.effect.size)) *
          ifelse(reverse,
                 1 - max(v$qvalue, 10 ^ -16),
                 max(v$qvalue, 10 ^ -16))
        ESs = c(ESs, v$sex.standardised.effect.size)
        testcounter = testcounter + 1
      } else{
        if (nrow(v2) > 0) {
          v100 = paste(
            'Less than ',
            minRequiredPerSex,
            ' participants per sex. ',
            ifelse(nrow(v2) > 1, paste0('X', nrow(v2), collapse = ''), ''),
            sep = '',
            collapse = ''
          )
        } else{
          v100 = 'Not_tested'
        }
      }
      r = c(r, v100)
    }
    
    # Calculate scores and create a table of results
    SumSigES = sum(sign(ESs))
    if (abs(SumSigES) > 3) {
      fscore = score / ((testcounter ^ 2) * sqrt(expES) * ifelse(reverse, 1 - qvalThresh, qvalThresh)) *
        max(testcounter, ifelse(late, 8, length(centres)) /  2)
    } else{
      fscore = 1
    }
    r = c(r, fscore, mean(ESs), sum(sign(ESs)))
    names(r) = c(as.character(centres), 'score', 'mean ES', 'Direction')
    if (counter <= 1) {
      output = r
    } else{
      output = rbind(output, r)
    }
    counter = counter + 1
  }
  
  
  output = as.data.frame(output)
  output$parameter = unique(agES$parameter)
  output = merge(x = output, y = effinal, by = 'parameter')
  return(output)
}

# Execute the driver function and write the output to a CSV file
write.csv(x = MainFunction(dfearly,
                     'HEM',
                     centres = unique(df$phenotyping_center)),
          file = 'SexualDimorphism_HEM_early.csv')
