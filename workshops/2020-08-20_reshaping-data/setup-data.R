library(tidyverse) # install.packages('tidyverse')

data <- read_csv('https://raw.githubusercontent.com/wvictor14/TOG/master/data/GSE98224.csv')
data

pdat <- data %>% select(expr_geo_id:ga_days)

expr <- data %>% select(contains('transcript'))

meth <- data %>% select(contains('cg'))

saveRDS(pdat, 'pdat.rds')
saveRDS(expr, 'expr.rds')
saveRDS(meth, 'meth.rds')