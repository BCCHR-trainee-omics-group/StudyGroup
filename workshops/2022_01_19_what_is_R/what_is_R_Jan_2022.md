---
title: "TOG: What is R? Workshop"
date: "Jan 19th 2022"
author: "Divya Anantsri"
output: 
  html_document: 
    keep_md: yes #keep the .md version of this .rmd
    theme: cosmo  #the themes your rmd-to-html document
    highlight: kate #your code highlight theme
---



***  

### 0.0 Introduction  

Welcome to the BCCHR Trainee 'Omics Group's `"What is R?"` workshop! 

In this **tutorial**, we are covering: 

1. installing R studio 
2. layout of R studio 
3. introduction to R markdown 
4. setting up - installing packages and loading libraries 
4. exploratory data analysis (EDA) using a sample dataset 


Our **goals** for the end of today are: 

1. for you to feel comfortable navigating R studio and working with R markdown  
2. to give you the ability to self-guide online tutorials or attend future workshops  
  
  
*Please also refer to the accompanying powerpoint slides for a list of useful resources!*    

***  

### 1.0 Getting help in R  


```r
help(mean)

?mean

vignette()
```

***  

### 2.0 Installing packages 


```r
install.packages("tidyverse")
#this a collection of packages that we will use for data manipulation and plotting 

install.packages("palmerpenguins")
#this package has the sample dataset we will use - Palmerpenguins dataset 
#this dataset is great for learning exploratory data analysis 
```

***  

### 3.0 Load the respective libraries


```r
library(tidyverse)
```

```
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
```

```
## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
## ✓ tibble  3.1.4     ✓ dplyr   1.0.7
## ✓ tidyr   1.1.3     ✓ stringr 1.4.0
## ✓ readr   2.0.1     ✓ forcats 0.5.1
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## x dplyr::filter() masks stats::filter()
## x dplyr::lag()    masks stats::lag()
```

```r
library(palmerpenguins)
```

***  

### 4.0 Look at the dataset (penguins)


```r
#making a new object/vector using the `penguins` dataset that is pre-built into the palmerpenguins package
penguins <- palmerpenguins::penguins

head(penguins)
```

```
## # A tibble: 6 × 8
##   species island bill_length_mm bill_depth_mm flipper_length_… body_mass_g sex  
##   <fct>   <fct>           <dbl>         <dbl>            <int>       <int> <fct>
## 1 Adelie  Torge…           39.1          18.7              181        3750 male 
## 2 Adelie  Torge…           39.5          17.4              186        3800 fema…
## 3 Adelie  Torge…           40.3          18                195        3250 fema…
## 4 Adelie  Torge…           NA            NA                 NA          NA <NA> 
## 5 Adelie  Torge…           36.7          19.3              193        3450 fema…
## 6 Adelie  Torge…           39.3          20.6              190        3650 male 
## # … with 1 more variable: year <int>
```

```r
#view the first few entries in your dataset

tail(penguins)
```

```
## # A tibble: 6 × 8
##   species island bill_length_mm bill_depth_mm flipper_length_… body_mass_g sex  
##   <fct>   <fct>           <dbl>         <dbl>            <int>       <int> <fct>
## 1 Chinst… Dream            45.7          17                195        3650 fema…
## 2 Chinst… Dream            55.8          19.8              207        4000 male 
## 3 Chinst… Dream            43.5          18.1              202        3400 fema…
## 4 Chinst… Dream            49.6          18.2              193        3775 male 
## 5 Chinst… Dream            50.8          19                210        4100 male 
## 6 Chinst… Dream            50.2          18.7              198        3775 fema…
## # … with 1 more variable: year <int>
```

```r
#view the last few entries

glimpse(penguins)
```

```
## Rows: 344
## Columns: 8
## $ species           <fct> Adelie, Adelie, Adelie, Adelie, Adelie, Adelie, Adel…
## $ island            <fct> Torgersen, Torgersen, Torgersen, Torgersen, Torgerse…
## $ bill_length_mm    <dbl> 39.1, 39.5, 40.3, NA, 36.7, 39.3, 38.9, 39.2, 34.1, …
## $ bill_depth_mm     <dbl> 18.7, 17.4, 18.0, NA, 19.3, 20.6, 17.8, 19.6, 18.1, …
## $ flipper_length_mm <int> 181, 186, 195, NA, 193, 190, 181, 195, 193, 190, 186…
## $ body_mass_g       <int> 3750, 3800, 3250, NA, 3450, 3650, 3625, 4675, 3475, …
## $ sex               <fct> male, female, female, NA, female, male, female, male…
## $ year              <int> 2007, 2007, 2007, 2007, 2007, 2007, 2007, 2007, 2007…
```

```r
#view the structure of each variable/vector in the dataset
```

***  

### 5.0 What is the structure of this dataset?


```r
str(penguins)
```

```
## tibble [344 × 8] (S3: tbl_df/tbl/data.frame)
##  $ species          : Factor w/ 3 levels "Adelie","Chinstrap",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ island           : Factor w/ 3 levels "Biscoe","Dream",..: 3 3 3 3 3 3 3 3 3 3 ...
##  $ bill_length_mm   : num [1:344] 39.1 39.5 40.3 NA 36.7 39.3 38.9 39.2 34.1 42 ...
##  $ bill_depth_mm    : num [1:344] 18.7 17.4 18 NA 19.3 20.6 17.8 19.6 18.1 20.2 ...
##  $ flipper_length_mm: int [1:344] 181 186 195 NA 193 190 181 195 193 190 ...
##  $ body_mass_g      : int [1:344] 3750 3800 3250 NA 3450 3650 3625 4675 3475 4250 ...
##  $ sex              : Factor w/ 2 levels "female","male": 2 1 1 NA 1 2 1 2 NA NA ...
##  $ year             : int [1:344] 2007 2007 2007 2007 2007 2007 2007 2007 2007 2007 ...
```

```r
#str is the baseR fucntion of the glimpse fucntion we used above
```

***  

### 6.0 Using ggplot2 for exploratory analysis of the `palmerpenguins` dataset 

Code from https://www.reed.edu/data-at-reed/resources/R/scatterplots.html

ggplot2 syntax includes:  

1. data = dataset being plotted
2. aes = aesthetic mappings - specifiying the variables for the X and Y axes
3. geom = genometric objects - points, lines, bars 


```r
ggplot(data = penguins, 
       mapping = aes(x = bill_length_mm, y = body_mass_g)) +
  geom_point()
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

![](what_is_R_Jan_2022_files/figure-html/scatterplot-1.png)<!-- -->

```r
# using the tidyverse `pipe` to select our dataset and the ggplot `+` to add layers and specifications of plotting
# penguins %>% 
#   ggplot(aes(x = bill_length_mm, y = body_mass_g)) +
#   geom_point()
```

Add the visual variable of color to show patterns within each island 


```r
ggplot(data = penguins, 
       mapping = aes(x = bill_length_mm, y = body_mass_g, color = island)) +
  geom_point()
```

```
## Warning: Removed 2 rows containing missing values (geom_point).
```

![](what_is_R_Jan_2022_files/figure-html/colour-scatterplot-1.png)<!-- -->

```r
# using the tidyverse `pipe` to select our dataset and the ggplot `+` to add layers and specifications of plotting
# penguins %>% 
#   ggplot(aes(x = bill_length_mm, y = body_mass_g, colour = island)) +
#   geom_point()
```

***  

### 7.0 Store your project information in .Rproj   

Whenever you start a project or analysis, it is best practice to make a separate directory and sub-directories containing all of your files. An example would be:  

`project_folder`  
|  
|__ `data` (where you store all your data files; you can also divide the folder further into `raw`, `processed`, `external`)  
|  
|__ `scripts` (all your .RMD, .MD, and HTML/Word?PDF knitted files)  

You can then store all of the directory information in a `.Rproj` file (and you can open that .Rproj file whenever working on a project), by clicking on the blue cube with the R and green plus sign in the option bar above --> Existing Directory --> navigate to your `project_folder`. The top level, i.e. your `project_folder` will be set as your working directory for that specific .Rproj

***  

### 8.0 Knit your file!  


```r
#rmarkdown::render("what_is_R_Jan_2022.Rmd", "html_document")
```

Easier way to Knit:  

- Click on the Knit symbol dropdown arrow in the option bar above --> Knit Directory --> Project Directory (to ensure that you set the Knit parameter to your .Rproj folder directory)
- Click on the Knit symbol  

OR  

- Ctrl/Cmd + Shift + K  

***  

