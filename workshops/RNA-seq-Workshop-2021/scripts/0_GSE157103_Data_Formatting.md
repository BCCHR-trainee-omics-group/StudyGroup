---
title: "TOG RNA-seq Workshop" 
author: Nikita Telkar 
date: July 2021
output: 
  html_document: 
    keep_md: yes 
    theme: flatly  
    highlight: pygments 
---  

***  

### Raw data formatting  


Here are the steps which I used to make the _*formatted.txt* expression count and clinical/phenotype information tab-delimited files. I've added detailed annotation for all the functions/commands I've used.    

> The raw files can be downloaded from [GEO Accession GSE157103](https://www-ncbi-nlm-nih-gov.ezproxy.library.ubc.ca/geo/query/acc.cgi?acc=GSE157103).  





```r
# required packages
library(tidyverse)
library(here)

# importing the count matrix (gene expression) and clinical data file
eDat <- read.delim(here::here("data", "GSE157103_genes.ec.tsv"))
cDat <- read.delim(here::here("data", "GSE157103_series_matrix.txt"), sep = "\t")


# FORMATTING THE CLINCAL INFORMATION FILE

# removing the ! character at the start of every value in the X.Sample_title column of the cDat
# dataframe
cDat$Sample_title <- str_sub(cDat$Sample_title, 2)

# transposing rows to columns
cDat <- as.data.frame(t(cDat))

# setting the column names of cDat to be the values in the first row
colnames(cDat) <- cDat[1, ]

# removing first row of cDat
cDat <- cDat[-1, ]

# converting the rownames to a column called Sample
cDat <- cDat %>%
    rownames_to_column(var = "Sample")

# displaying the cDat dataframe as a tibble
as_tibble(cDat)
```

```
## # A tibble: 126 x 65
##    Sample    Sample_geo_acces… Sample_status   Sample_submissi… Sample_last_upd…
##    <chr>     <chr>             <chr>           <chr>            <chr>           
##  1 COVID_01… GSM4753021        Public on Aug … Aug 28 2020      Aug 29 2020     
##  2 COVID_02… GSM4753022        Public on Aug … Aug 28 2020      Aug 29 2020     
##  3 COVID_03… GSM4753023        Public on Aug … Aug 28 2020      Aug 29 2020     
##  4 COVID_04… GSM4753024        Public on Aug … Aug 28 2020      Aug 29 2020     
##  5 COVID_05… GSM4753025        Public on Aug … Aug 28 2020      Aug 29 2020     
##  6 COVID_06… GSM4753026        Public on Aug … Aug 28 2020      Aug 29 2020     
##  7 COVID_07… GSM4753027        Public on Aug … Aug 28 2020      Aug 29 2020     
##  8 COVID_08… GSM4753028        Public on Aug … Aug 28 2020      Aug 29 2020     
##  9 COVID_09… GSM4753029        Public on Aug … Aug 28 2020      Aug 29 2020     
## 10 COVID_10… GSM4753030        Public on Aug … Aug 28 2020      Aug 29 2020     
## # … with 116 more rows, and 60 more variables: Sample_type <chr>,
## #   Sample_channel_count <chr>, Sample_source_name_ch1 <chr>,
## #   Sample_organism_ch1 <chr>, Sample_characteristics_ch1 <chr>,
## #   Sample_characteristics_ch1.1 <chr>, Sample_characteristics_ch1.2 <chr>,
## #   Sample_characteristics_ch1.3 <chr>, Sample_characteristics_ch1.4 <chr>,
## #   Sample_characteristics_ch1.5 <chr>, Sample_characteristics_ch1.6 <chr>,
## #   Sample_characteristics_ch1.7 <chr>, Sample_characteristics_ch1.8 <chr>,
## #   Sample_characteristics_ch1.9 <chr>, Sample_characteristics_ch1.10 <chr>,
## #   Sample_characteristics_ch1.11 <chr>, Sample_characteristics_ch1.12 <chr>,
## #   Sample_characteristics_ch1.13 <chr>, Sample_characteristics_ch1.14 <chr>,
## #   Sample_characteristics_ch1.15 <chr>, Sample_characteristics_ch1.16 <chr>,
## #   Sample_characteristics_ch1.17 <chr>, Sample_characteristics_ch1.18 <chr>,
## #   Sample_characteristics_ch1.19 <chr>, Sample_characteristics_ch1.20 <chr>,
## #   Sample_treatment_protocol_ch1 <chr>, Sample_growth_protocol_ch1 <chr>,
## #   Sample_molecule_ch1 <chr>, Sample_extract_protocol_ch1 <chr>,
## #   Sample_extract_protocol_ch1.1 <chr>, Sample_taxid_ch1 <chr>,
## #   Sample_description <chr>, Sample_description.1 <chr>,
## #   Sample_description.2 <chr>, Sample_data_processing <chr>,
## #   Sample_data_processing.1 <chr>, Sample_data_processing.2 <chr>,
## #   Sample_data_processing.3 <chr>, Sample_data_processing.4 <chr>,
## #   Sample_platform_id <chr>, Sample_contact_name <chr>,
## #   Sample_contact_email <chr>, Sample_contact_department <chr>,
## #   Sample_contact_institute <chr>, Sample_contact_address <chr>,
## #   Sample_contact_city <chr>, Sample_contact_state <chr>,
## #   Sample_contact_zip/postal_code <chr>, Sample_contact_country <chr>,
## #   Sample_data_row_count <chr>, Sample_instrument_model <chr>,
## #   Sample_library_selection <chr>, Sample_library_source <chr>,
## #   Sample_library_strategy <chr>, Sample_relation <chr>,
## #   Sample_relation.1 <chr>, Sample_supplementary_file_1 <chr>,
## #   series_matrix_table_begin <chr>, D_REF <chr>, series_matrix_table_end <chr>
```

```r
# displaying the column names of cDat
names(cDat)
```

```
##  [1] "Sample"                         "Sample_geo_accession"          
##  [3] "Sample_status"                  "Sample_submission_date"        
##  [5] "Sample_last_update_date"        "Sample_type"                   
##  [7] "Sample_channel_count"           "Sample_source_name_ch1"        
##  [9] "Sample_organism_ch1"            "Sample_characteristics_ch1"    
## [11] "Sample_characteristics_ch1.1"   "Sample_characteristics_ch1.2"  
## [13] "Sample_characteristics_ch1.3"   "Sample_characteristics_ch1.4"  
## [15] "Sample_characteristics_ch1.5"   "Sample_characteristics_ch1.6"  
## [17] "Sample_characteristics_ch1.7"   "Sample_characteristics_ch1.8"  
## [19] "Sample_characteristics_ch1.9"   "Sample_characteristics_ch1.10" 
## [21] "Sample_characteristics_ch1.11"  "Sample_characteristics_ch1.12" 
## [23] "Sample_characteristics_ch1.13"  "Sample_characteristics_ch1.14" 
## [25] "Sample_characteristics_ch1.15"  "Sample_characteristics_ch1.16" 
## [27] "Sample_characteristics_ch1.17"  "Sample_characteristics_ch1.18" 
## [29] "Sample_characteristics_ch1.19"  "Sample_characteristics_ch1.20" 
## [31] "Sample_treatment_protocol_ch1"  "Sample_growth_protocol_ch1"    
## [33] "Sample_molecule_ch1"            "Sample_extract_protocol_ch1"   
## [35] "Sample_extract_protocol_ch1.1"  "Sample_taxid_ch1"              
## [37] "Sample_description"             "Sample_description.1"          
## [39] "Sample_description.2"           "Sample_data_processing"        
## [41] "Sample_data_processing.1"       "Sample_data_processing.2"      
## [43] "Sample_data_processing.3"       "Sample_data_processing.4"      
## [45] "Sample_platform_id"             "Sample_contact_name"           
## [47] "Sample_contact_email"           "Sample_contact_department"     
## [49] "Sample_contact_institute"       "Sample_contact_address"        
## [51] "Sample_contact_city"            "Sample_contact_state"          
## [53] "Sample_contact_zip/postal_code" "Sample_contact_country"        
## [55] "Sample_data_row_count"          "Sample_instrument_model"       
## [57] "Sample_library_selection"       "Sample_library_source"         
## [59] "Sample_library_strategy"        "Sample_relation"               
## [61] "Sample_relation.1"              "Sample_supplementary_file_1"   
## [63] "series_matrix_table_begin"      "D_REF"                         
## [65] "series_matrix_table_end"
```

```r
# making a new dataframe that contains only the clinical variables of our interest called pDat
pDat <- cDat %>%
    dplyr::select(Sample, Sample_characteristics_ch1.1, Sample_characteristics_ch1.2, Sample_characteristics_ch1.3,
        Sample_characteristics_ch1.4, Sample_characteristics_ch1.5, Sample_characteristics_ch1.6, Sample_characteristics_ch1.7,
        Sample_characteristics_ch1.9, Sample_characteristics_ch1.10, Sample_characteristics_ch1.11, Sample_characteristics_ch1.13,
        Sample_characteristics_ch1.14, Sample_characteristics_ch1.15)

# dimensions (no.of rows and columns) in our new pDat dataframe
dim(pDat)
```

```
## [1] 126  14
```

```r
# renaming columns 2 to 4 of our new pDat dataframe
colnames(pDat)[2:14] <- c("Age", "Sex", "ICU", "APACHEII_Score", "Charlson_Score", "Mechanical_Ventilation",
    "Ventilator_free_days", "Hospital_free_days_post_45_days", "Ferritin_ng.ml", "CRP_mg.l", "Procalcitonin_ng.ml",
    "Lactate_mmol.l", "Fibrinogen_mg.dL")

# displaying first 6 rows of pDat
head(pDat)
```

```
##                     Sample             Age       Sex     ICU    APACHEII_Score
## 1 COVID_01_39y_male_NonICU age (years): 39 Sex: male icu: no      apacheii: 15
## 2 COVID_02_63y_male_NonICU age (years): 63 Sex: male icu: no apacheii: unknown
## 3 COVID_03_33y_male_NonICU age (years): 33 Sex: male icu: no apacheii: unknown
## 4 COVID_04_49y_male_NonICU age (years): 49 Sex: male icu: no apacheii: unknown
## 5 COVID_05_49y_male_NonICU age (years): 49 Sex: male icu: no      apacheii: 19
## 6  COVID_06_.y_male_NonICU  age (years): : Sex: male icu: no apacheii: unknown
##      Charlson_Score      Mechanical_Ventilation     Ventilator_free_days
## 1 charlson score: 0 mechanical ventilation: yes  ventilator-free days: 0
## 2 charlson score: 2  mechanical ventilation: no ventilator-free days: 28
## 3 charlson score: 2  mechanical ventilation: no ventilator-free days: 28
## 4 charlson score: 1  mechanical ventilation: no ventilator-free days: 28
## 5 charlson score: 1 mechanical ventilation: yes ventilator-free days: 23
## 6 charlson score: 1  mechanical ventilation: no ventilator-free days: 28
##                      Hospital_free_days_post_45_days         Ferritin_ng.ml
## 1  hospital-free days post 45 day followup (days): 0  ferritin (ng/ml): 946
## 2 hospital-free days post 45 day followup (days): 39 ferritin (ng/ml): 1060
## 3 hospital-free days post 45 day followup (days): 18 ferritin (ng/ml): 1335
## 4 hospital-free days post 45 day followup (days): 39  ferritin (ng/ml): 583
## 5 hospital-free days post 45 day followup (days): 27  ferritin (ng/ml): 800
## 6 hospital-free days post 45 day followup (days): 36  ferritin (ng/ml): 563
##              CRP_mg.l         Procalcitonin_ng.ml            Lactate_mmol.l
## 1    crp (mg/l): 73.1   procalcitonin (ng/ml): 36     lactate (mmol/l): 0.9
## 2 crp (mg/l): unknown procalcitonin (ng/ml): 0.37 lactate (mmol/l): unknown
## 3    crp (mg/l): 53.2 procalcitonin (ng/ml): 0.07 lactate (mmol/l): unknown
## 4   crp (mg/l): 251.1 procalcitonin (ng/ml): 0.98    lactate (mmol/l): 0.87
## 5   crp (mg/l): 355.8 procalcitonin (ng/ml): 4.92    lactate (mmol/l): 1.48
## 6   crp (mg/l): 129.1 procalcitonin (ng/ml): 0.67    lactate (mmol/l): 0.86
##      Fibrinogen_mg.dL
## 1     fibrinogen: 513
## 2 fibrinogen: unknown
## 3     fibrinogen: 513
## 4     fibrinogen: 949
## 5     fibrinogen: 929
## 6     fibrinogen: 769
```

```r
# only keeping the last two characters of all values in the Age column which is the age of the
# patients and removing all alphabetical charaters (age (years): )
pDat$Age <- str_sub(pDat$Age, -2)

# removing first 6 characters of all values in the Sex column (sex: ), to keep only keep male or
# female
pDat$Sex <- str_sub(pDat$Sex, 6)

# removing first 6 characters of all values in the ICU column (icu: ), to keep only yes or no
pDat$ICU <- str_sub(pDat$ICU, 6)

# removing first 11 characters of all values in the APACHEII_Score column (apacheii: ), to keep
# only the score
pDat$APACHEII_Score <- str_sub(pDat$APACHEII_Score, 11)

# removing the first 17 characters of all values in the Charlson_Score(charlson score: ), to only
# keep the score
pDat$Charlson_Score <- str_sub(pDat$Charlson_Score, 17)

# removing the first 25 characters of all values in the Mechanical_Ventilation(mechanical
# ventilation: ), to only keep yes or no
pDat$Mechanical_Ventilation <- str_sub(pDat$Mechanical_Ventilation, 25)

# removing the first 23 characters of all values in the Ventilator_free_days(ventilator-free days:
# ), to only keep number
pDat$Ventilator_free_days <- str_sub(pDat$Ventilator_free_days, 23)

# finding out length of string (i.e. characters present) in hospital-free days post 45 day followup
# (days): in the Hospital_free_days_post_45_days variables as it is too many to count make new
# value capturing the 1st value int the Hospital_free_days_post_45_days column
value <- pDat$Hospital_free_days_post_45_days[1]
# find number of characters
nchar(value)
```

```
## [1] 49
```

```r
# removing to keep only the number
pDat$Hospital_free_days_post_45_days <- str_sub(pDat$Hospital_free_days_post_45_days, 49)

# removing the first 19 characters of all values in the Ferritin_ng.ml to only keep value
value <- pDat$Ferritin_ng.ml[1]
nchar(value)
```

```
## [1] 21
```

```r
pDat$Ferritin_ng.ml <- str_sub(pDat$Ferritin_ng.ml, 19)

# removing the first 13 characters of all values in the CRP_mg.l to only keep value
value <- pDat$CRP_mg.l[1]
nchar(value)
```

```
## [1] 16
```

```r
pDat$CRP_mg.l <- str_sub(pDat$CRP_mg.l, 13)

# removing the first 24 characters of all values in the Procalcitonin_ng.ml to only keep value
value <- pDat$Procalcitonin_ng.ml[1]
nchar(value)
```

```
## [1] 25
```

```r
pDat$Procalcitonin_ng.ml <- str_sub(pDat$Procalcitonin_ng.ml, 24)

# removing the first 19 characters of all values in the Lactate_mmol.l to only keep value
value <- pDat$Lactate_mmol.l[1]
nchar(value)
```

```
## [1] 21
```

```r
pDat$Lactate_mmol.l <- str_sub(pDat$Lactate_mmol.l, 19)

# removing the first 13 characters of all values in the Fibrinogen_mg.dL to only keep value
value <- pDat$Fibrinogen_mg.dL[1]
nchar(value)
```

```
## [1] 15
```

```r
pDat$Fibrinogen_mg.dL <- str_sub(pDat$Fibrinogen_mg.dL, 13)

str(pDat)
```

```
## 'data.frame':	126 obs. of  14 variables:
##  $ Sample                         : chr  "COVID_01_39y_male_NonICU" "COVID_02_63y_male_NonICU" "COVID_03_33y_male_NonICU" "COVID_04_49y_male_NonICU" ...
##  $ Age                            : chr  "39" "63" "33" "49" ...
##  $ Sex                            : chr  "male" "male" "male" "male" ...
##  $ ICU                            : chr  "no" "no" "no" "no" ...
##  $ APACHEII_Score                 : chr  "15" "unknown" "unknown" "unknown" ...
##  $ Charlson_Score                 : chr  "0" "2" "2" "1" ...
##  $ Mechanical_Ventilation         : chr  "yes" "no" "no" "no" ...
##  $ Ventilator_free_days           : chr  "0" "28" "28" "28" ...
##  $ Hospital_free_days_post_45_days: chr  "0" "39" "18" "39" ...
##  $ Ferritin_ng.ml                 : chr  "946" "1060" "1335" "583" ...
##  $ CRP_mg.l                       : chr  "73.1" "unknown" "53.2" "251.1" ...
##  $ Procalcitonin_ng.ml            : chr  "36" "0.37" "0.07" "0.98" ...
##  $ Lactate_mmol.l                 : chr  "0.9" "unknown" "unknown" "0.87" ...
##  $ Fibrinogen_mg.dL               : chr  "513" "unknown" "513" "949" ...
```

```r
# converting all measurement values to numeric
pDat <- pDat %>%
    mutate_at(vars(-Sample, -Sex, -ICU, -Mechanical_Ventilation), as.numeric)

str(pDat)
```

```
## 'data.frame':	126 obs. of  14 variables:
##  $ Sample                         : chr  "COVID_01_39y_male_NonICU" "COVID_02_63y_male_NonICU" "COVID_03_33y_male_NonICU" "COVID_04_49y_male_NonICU" ...
##  $ Age                            : num  39 63 33 49 49 NA 38 78 64 62 ...
##  $ Sex                            : chr  "male" "male" "male" "male" ...
##  $ ICU                            : chr  "no" "no" "no" "no" ...
##  $ APACHEII_Score                 : num  15 NA NA NA 19 NA NA 43 31 34 ...
##  $ Charlson_Score                 : num  0 2 2 1 1 1 7 7 2 1 ...
##  $ Mechanical_Ventilation         : chr  "yes" "no" "no" "no" ...
##  $ Ventilator_free_days           : num  0 28 28 28 23 28 28 0 0 2 ...
##  $ Hospital_free_days_post_45_days: num  0 39 18 39 27 36 42 0 0 0 ...
##  $ Ferritin_ng.ml                 : num  946 1060 1335 583 800 ...
##  $ CRP_mg.l                       : num  73.1 NA 53.2 251.1 355.8 ...
##  $ Procalcitonin_ng.ml            : num  36 0.37 0.07 0.98 4.92 0.67 0.06 4.2 0.25 1.1 ...
##  $ Lactate_mmol.l                 : num  0.9 NA NA 0.87 1.48 0.86 1.17 1.65 1.14 1.45 ...
##  $ Fibrinogen_mg.dL               : num  513 NA 513 949 929 769 478 780 317 893 ...
```

```r
# converting all NA values to 0
pDat <- pDat %>%
    replace(is.na(.), 0)

# changing the name of column 1 to gene
colnames(eDat)[1] <- "gene"

# converting the column titled gene to rownames of the eDat dataframe
eDat <- eDat %>%
    column_to_rownames(var = "gene")

# making a new column called ID in the pDat dataframe and adding the column names as sample IDs
# (columns are in order of row)
pDat$ID <- colnames(eDat)

# making ID as column 1
pDat <- pDat %>%
    dplyr::select(ID, everything())

# making a new column called COVID by coopying the first letter of the ID column
pDat <- pDat %>%
    mutate(COVID = substr(ID, 1L, 1L))

# recoding the COVID column: if C, change value to yes; if N, change value to no
pDat <- pDat %>%
    mutate(COVID = recode(COVID, C = "yes", N = "no"))

# making COVID as column 5
pDat <- pDat %>%
    dplyr::select(ID, Sample, Age, Sex, COVID, everything())

# converting the rownames to a column titled gene - we will be saving this as a delimited tab file,
# and so saving the rownames is not an available option.
eDat <- eDat %>%
    rownames_to_column(var = "gene")

# saving both eDat and pDat as tab delimted text files
write_delim(pDat, file = "GSE157103_formatted_pDat.txt", delim = "\t")
write_delim(eDat, file = "GSE157103_formatted_eDat.txt", delim = "\t")
```
