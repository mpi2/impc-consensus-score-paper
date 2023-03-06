load(file = 'Rawdata_IMPC_HEM_onlyControls.Rdata')
df = subset(df, df$biological_sample_group == 'control')
# Filter for the date of experiment
df = df[df$date_of_experiment >= '2018-01-01T00:00:00Z' &
          df$date_of_experiment <= '2020-12-31T00:00:00Z', ]

# for the scoring paper
tblss=table(df$genetic_background,df$phenotyping_center)
tblss[tblss>0]=1
tblss[tblss<=0]='-'
tblss[tblss>0]='✓'
tblss

write.csv(tblss,file='Supplementry table X. BackgroundStrainVariationInIMPC.csv')
