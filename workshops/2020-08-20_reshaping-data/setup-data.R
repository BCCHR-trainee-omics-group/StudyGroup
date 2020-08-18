library(tidyverse) # install.packages('tidyverse')

data <- read_csv('https://raw.githubusercontent.com/wvictor14/TOG/master/data/GSE98224.csv')
data

pdat <- data %>% select(expr_geo_id:ga_days)

expr <- data %>% select(expr_geo_id, contains('transcript'))

meth <- data %>% select(meth_geo_id, contains('cg'))
here::here()
write_csv(pdat, here::here('workshops', '2020-08-20_reshaping-data', 'pdat.csv'))
write_csv(expr, here::here('workshops', '2020-08-20_reshaping-data', 'expr.csv'))
write_csv(meth, here::here('workshops', '2020-08-20_reshaping-data', 'meth.csv'))
