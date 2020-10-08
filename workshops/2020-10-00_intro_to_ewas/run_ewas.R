## @knitr load_packages
library(GEOquery)
library(Biobase)
library(limma)
library(ggplot2)

## @knitr load_data

eset <- getGEO(file="GSE100197_series_matrix.txt.gz",)
metadata <- pData(eset)
colnames(metadata)<- gsub(":ch1| ","",colnames(metadata))
metadata <- metadata[!metadata$pathologygroup == "REPLICATE",]
metadata$gestationalage <- as.numeric(metadata$gestationalage)
metadata$int_grp <- paste0(metadata$pathologygroup,metadata$fetalsex)
metadata$status <- sapply( metadata$pathologygroup, function(x)if(x== "Term" | x == "PreT"){return("CONTROL")} else{return(x)})
methy <- exprs(eset)
matched <- match(rownames(metadata),colnames(methy))
methy <- methy[,matched]


## @knitr run_fit
model <- model.matrix(~ 0+ pathologygroup + fetalsex, data=metadata)
fit <- lmFit(methy, model)
fit <- eBayes(fit)
contrasts <- makeContrasts(
preTEOPE=pathologygroupPreT-pathologygroupEOPE,
preTLOPE=pathologygroupPreT-pathologygroupLOPE,
preTIUGR=pathologygroupPreT-pathologygroupIUGR,
termEOPE=pathologygroupTerm-pathologygroupEOPE,
termLOPE=pathologygroupTerm-pathologygroupLOPE,
termIUGR=pathologygroupTerm-pathologygroupIUGR,
sex=fetalsexMALE,levels=model)
fitCont <- contrasts.fit(fit,contrasts)
fitCont <- eBayes(fitCont)

## @knitr run_age_fit

age_model <- model.matrix(~ 0+ status + fetalsex + gestationalage,  data=metadata)
age_fit <- lmFit(methy, age_model)
age_fit <- eBayes(age_fit)
contrasts <- makeContrasts(
  EOPE=statusCONTROL-statusEOPE,
  LOPE=statusCONTROL-statusLOPE,
  IUGR=statusCONTROL-statusIUGR,
  sex=fetalsexMALE,levels=age_model)
age_fitCont <- contrasts.fit(age_fit,contrasts)
age_fitCont <- eBayes(age_fitCont)

## @knitr run_interact_fit

interact_model <- model.matrix(~ 0+int_grp,  data=metadata)
interact_fit <- lmFit(methy, interact_model)
interact_fit <- eBayes(interact_fit)
contrasts <- makeContrasts(
  EOPE_int = (int_grpPreTMALE- int_grpEOPEMALE) - (int_grpPreTFEMALE - int_grpEOPEFEMALE),
  EOPE_male = (int_grpPreTMALE- int_grpEOPEMALE),
  EOPE_female = (int_grpPreTMALE- int_grpEOPEMALE),
  levels=interact_model)
interact_fitCont <- contrasts.fit(interact_fit,contrasts)
interact_fitCont <- eBayes(interact_fitCont)



