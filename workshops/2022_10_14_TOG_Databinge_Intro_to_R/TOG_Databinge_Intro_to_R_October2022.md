---
title: "TOG + Databinge: Intro to R"
date: "14 October 2022"
output:
  html_document:
    keep_md: TRUE
    toc: true
    toc_float: true
    code_download: true
    theme: flatly
    highlight: haddock
---



***

## 0.0 Introduction  

This script is adapted from from the British Columbia Children's Hospital Research Institute Trainee 'Omics Group's [Intro to R](https://github.com/BCCHR-trainee-omics-group/StudyGroup/blob/master/workshops/2021-07-24_intro-to-R/intro-R-workshop.Rmd) and [Jordan Sicherman's](https://github.com/jsicherman/WorkshopR/blob/master/workshops/01_getting-started.Rmd) and has been reworked by [Will Cassaza](https://github.com/wilcas) and [Nikita Telkar](https://github.com/nikita-telkar/) for the TOG + Databinge Intro to Programming Workshops 2022 

***  

## 1.0 Outline and Goals of the Workshop  

Today, we are going to:  

- Introduce R and RStudio  
- Basic R syntax  
- Data Types within R  
- Working with Your Data  

***  

## 2.0 Introducing R and RStudio  

### 2.1 Downlaod R and RStudio  

> 1.  Download R from [CRAN](https://cran.r-project.org/)  
> 2.  Download [RStudio](https://rstudio.com/products/rstudio/download/preview/)  

### 2.1  What is R?  

You've probably heard people talking about R before, but *what really is R*?  

R is a *programming language and environment for statistical computing and graphics*. The R language was originally designed for statistics, and is now also frequently used for data analysis and visualization. R is free and open source, and the scientific community has designed many great freely available R tools/packages that are useful for all types of data analyses, ranging from finance to biological and genomic analyses.  


### 2.2 What is RStudio?  

R is a *programming language*, which means it is a language you can use to give instructions to your computer. R can be used in the *command line interface (CLI)*, but it is easier to use the R language in the RStudio application, which is an *integrated development environment (IDE)* that you can use to write code, run it line-by-line, view your figures, and save code for later use.   

The RStudio IDE now also supports not just the R language, but also Python and Julia!  

**In almost all cases, people who use 'R', actually work in RStudio the majority of the time, however you can very well use the CLI/Terminal, too.**

Notice that when you open RStudio, it is divided into different 'panes'. For now we will focus on the *console* (left side).  


### 2.3 The R console  


The console is your view into the R engine. You can use it to give commands to R and immediately see the output.  The R Console is synonymous to the Terminal within a CLI.  

The `>` symbol in the console is called the *prompt*. It is inviting you to type R commands, such as `2 * 2` or `mean(x)`. When you press Enter(Windows) / Return(Mac), you can see the output printed by R.  

> Try this out in your own console!  

Sometimes, you may accidentally press <kbd>Return</kbd> before finishing a complete command --> This often happens because you forgot to close parentheses.  

For example, you may have written: `mean(c(1,2,3)`; this command is incomplete because it is missing one `)` at the end. When you give the console an incomplete command, R expects you to complete it. Therefore, the prompt (`>`) is replaced by a plus sign (`+`). This is R's way of telling you to continue writing. When this happens, you have two alternatives. Either you complete the command (e.g., write the missing `)`) and then press <kbd>Return</kbd> *or* you can press <kbd>Esc</kbd>, to cancel the command and show the prompt.  

*A word of caution: auto-saving of the R environment*  

It is also best practice to disable RStudio's environment auto-saving option. Do the following:  

  - Click `Tools` > `Global Options`  
  - Where it says `Save workspace to .RData on exit`, select the option `Never`  
  - Click `OK` and close the dialog  

Now R will never save or prompt you to save your environment from RStudio.  

***  

### 2.4 R Scripts / RMarkdown  

You'll see that the extension of this file that is open is `.Rmd`, meaning a RMarkdown script. When you have a lot of code that you want to save, you can save all of your code within an R Script (`.R` extension) or a RMarkdown script.  

The difference between an R Script and a RMarkdown:  

- A R script --> used to run code --> output is a text file  
- RMarkdown --> used to run and document code --> several different outputs --> RMarkdown allows you to 'knit' your script and export your code as an interactive HTML document or write your thesis in RMarkdown and export it as a Word or PDF document  

***  

## 3.0 Basic R Syntax  

You write your R code within a code chunk in a RMarkdown. Use the shortcut Ctrl+Shift+Enter/Command+Shift+Enter to make a new chunk. It's always best practice to name your chunk as well  

You **assign values** to a certain variable or object using the `<-` assignment operator  


```r
a <- 3
a
```

```
## [1] 3
```

To perform **operations** on the value stored in that variable: 


```r
a+2
```

```
## [1] 5
```

```r
a*2
```

```
## [1] 6
```

```r
a^2
```

```
## [1] 9
```

```r
log(a)
```

```
## [1] 1.098612
```

But what if you want to assign multiple values to a single variable? We can use the `c` argument to combine values  


```r
b <- c(2, 5.5) 
b
```

```
## [1] 2.0 5.5
```

Now, you might have noticed that the output contains a [1] as the starting index - which means that the first value is designated at position 1. This is in contrast to Python, which uses 0-based indexing, meaning that the number 2 in the above example will be designated at position 0. *This is very important when you want to subset or split your data according to the cell, row, or column number.*  

Think of it like this: in North America, the entry-level floor of a building is called the First Floor or the Ground Floor and the floor directly above that one is called the Second Floor. Whereas in Europe or the U.K., that entry-level floor is called the Ground Floor, and the floor above is called the First Floor. Depending on which country you are in, going to the Second Floor means something different.  

<br> 

Now, imagine you're grading a class test, and calculating the percentage each student got. If you just have some code written and give it off to your colleague, they might not understand what you've done. For this reason, to ensure clarity and reproducibility, it is essential to **comment or document your code**  

You can do this within a code chunk using the `#` symbol  


```r
# grades of 5 students in the class
c <- c(32, 48, 45, 38, 41)

# highest score achievable = 50
d <- 50

# percentage grades for each student
(c/d)*100
```

```
## [1] 64 96 90 76 82
```

**Syntax Table**  

| Description | Symbol | Example |  
|:--|:-|:-----|  
| Add | `+` | `2 + 2` |  
| Subtract | `-` | `4 - 1` |  
| Multiply | `*` | `4 * 2` |  
| Power | `^` | `2^8` |  
| Divide | `/` | `5 / 3` |  
| Modulus | `%%` | `6 %% 2` |  
| Equality | `==` | `6 + 2 == 8` |  
| Inequality | `!=` | `TRUE != FALSE` |  
| Less than | `<` | `2 < 9` |  
| Greater than | `>` | `9 > 2` |  
| Less than or equal | `<=`  | `pi <= 5` |  
| Greater than or equal | `>=` | `5 * 2 >= 2` |  
| Short-circuiting AND | `&&` | `TRUE && FALSE` |  
| Short-circuiting OR | `||` | `FALSE || TRUE` |  
| Element-wise AND | `&` | `c(T, T, F, T) & c(T, F, F, T)` |  
| Element-wise OR | `|` | `c(T, F, F, T) | c(F, T, T, F)` |  
| Assignment | `<-` or `->` | `something <- 10` or `10 -> something` |  
| Vector indexing | `[]` or `$` | `c("A", "B", "C", "D")[1]` or `data.frame(A = 1:5, B = 6:10)$A` |  
| List indexing | `[[]]`  | `list(A = c(1, 2, 3), B = c(4, 5, 6))[[1]]]` |  
| Checking help | `?` | `?function_name` |   


***  

### 4.0 Data Types  

R has a few data types:  

- Numeric  
- Character  
- Factor  
- Logical  
- Vectors  

You can check the data type by using the `class` argument  

At the end, we'll also combine a few of the vectors into a `data.frame`  

### 4.1 Numeric  


```r
# whole numbers = integers
e <- c(5, 7.2)

class(e)
```

```
## [1] "numeric"
```

### 4.2 Character 

To make an object that stores letters and words, we enclose them in double(`"`) or single(`'`) inverted commas  


```r
pizzas <- c("pepperoni", "pesto chicken", "margherita")
class(pizzas)
```

```
## [1] "character"
```

More often than not, when you load data in external data (e.g., from Excel) numbers get assigned as characters, and hence, numeric operations cannot be perfromed on them. In this case, you will have to convert them into numeric variables  


```r
number_as_char <- c("4", "8", "16")
class(number_as_char)
```

```
## [1] "character"
```

```r
# 2*number_as_char

char_to_number <- as.numeric(number_as_char)
2*char_to_number
```

```
## [1]  8 16 32
```

### 4.3 Factors

A sequence of characters, given a specific order


```r
planets <- factor(c("Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"))
planets
```

```
## [1] Mercury Venus   Earth   Mars    Jupiter Saturn  Uranus  Neptune
## Levels: Earth Jupiter Mars Mercury Neptune Saturn Uranus Venus
```

Now, here, R automatically encodes the levels or the order of the character values by ascending alphabetical order. We specific Mercury as the first value, but the output shows that Earth was the assigned as the first level. It is possible to relevel the order of the factors using the `fct_relvel` argument, and you can explore that in your own time!  

### 4.4 Logical  

A TRUE or FALSE value  


```r
f <- 10

f > 5
```

```
## [1] TRUE
```

```r
f < 5
```

```
## [1] FALSE
```

### 4.5 Vectors  

A sequence of elements of the same class type 


```r
pets <- c("dog", "cat", "bird", "fish")

legs <- c(4, 4, 2, 0)

mammal <- c(TRUE, TRUE, FALSE, FALSE)
```

### 4.6 Data Frame  

Now, let's say we want to combine the above vectors we have in a single object. We can do so two ways:  

1. Make the vectore separately as above  
2. Or make them within the `data.frame` argument


```r
# Option 1  
animals <- data.frame(pets, legs, mammal)

# Option 2 
animals <- data.frame("pets" = c("dog", "cat", "bird", "fish"),
                      "legs" = c(4, 4, 2, 0),
                      "mammal" = c(TRUE, TRUE, FALSE, FALSE))
```

Note however that R only allows to add vectors to a dataframe that have the same number of entries, i.e., you cannot have a vector with 4 values and another with 3 values  

***

## 5.0 Working with Your Data

The following are some of the best practices when starting with your own data:  

1. Make a R Project  

> To make a new project, click on `File` --> `New Project`, and either choose an existing directory or make a new directory  


2. Make separate directories/folders for your scripts, data, etc.  

`project_folder.Rproj`   
|  
|__ `data`  
|  
|__ `scripts`


3. Use relative rather than absolute paths 

In a Mac, the path to a file might look like:  

```
/home/user/intro_to_r/data/Letter.txt
```

and in Windows, it could look like:  

```
C:\user\intro_to_r\data\Letter.txt
```

And so to ensure compatibility and reproducibility, it's always best to set the root of your project folder as your *working directory*. 

You can get the current working directory with `getwd()`, and set the working directory with `setwd("path/to/folder")` but we're going to use the `{here}` package. When you open your `.RProj` Project file, the `{here}` package automatically assignes the directory or folder the Project file lives in becomes your root directory 

Let's say you have the `Letter.txt` file within your data folder, you would access it as:


```r
letter <- read.delim(here::here("data", "Letter.txt")
```

But, given that we don't have such a folder at the moment, let's load in some data available online


```r
# Link to data in csv format
path <- "https://raw.githubusercontent.com/EDUCE-UBC/educer/main/data-raw/data_intro_ws.csv" 

# Read csv file
dat <- read.csv(path) 

class(dat)
```

```
## [1] "data.frame"
```

### 5.1 Viewing Data  

To check some specifics of this data frame:  


```r
# top 6 entries
head(dat)
```

```
##   Season Depth_m   O2_uM Add_data
## 1   Fall      10 203.533     TRUE
## 2   Fall      20 183.787    FALSE
## 3   Fall      40 130.579    FALSE
## 4   Fall      60  91.115     TRUE
## 5   Fall      75  69.828    FALSE
## 6   Fall      85  26.972    FALSE
```

```r
# assign these to a new object
top_dat <- head(dat)

# bottom 6 entries
tail(dat)
```

```
##    Season Depth_m  O2_uM Add_data
## 27 Summer     120 32.354     TRUE
## 28 Summer     135 20.446    FALSE
## 29 Summer     150  0.000    FALSE
## 30 Summer     165  0.000    FALSE
## 31 Summer     185  0.000    FALSE
## 32 Summer     200  0.000     TRUE
```

```r
# value within the cell at the second row and third column 
dat[2,3]
```

```
## [1] 183.787
```

```r
# column names of the dataframe
names(dat)
```

```
## [1] "Season"   "Depth_m"  "O2_uM"    "Add_data"
```

```r
# classes of each of the columns/vectors
str(dat)
```

```
## 'data.frame':	32 obs. of  4 variables:
##  $ Season  : chr  "Fall" "Fall" "Fall" "Fall" ...
##  $ Depth_m : int  10 20 40 60 75 85 90 97 100 110 ...
##  $ O2_uM   : num  203.5 183.8 130.6 91.1 69.8 ...
##  $ Add_data: logi  TRUE FALSE FALSE TRUE FALSE FALSE ...
```

### 5.2 Subsetting Data  

We recommend using the `tidyverse` package to perform all of your data transformations rather than base R, as the syntax of the functions within this package are more intuitive to understand when you're just starting to learn R.  


```r
# install.packages(tidyverse)
library(tidyverse)
```

```
## Warning: package 'tidyverse' was built under R version 4.2.1
```

```
## -- Attaching packages --------------------------------------- tidyverse 1.3.2 --
## v ggplot2 3.3.6      v purrr   0.3.4 
## v tibble  3.1.8      v dplyr   1.0.10
## v tidyr   1.2.1      v stringr 1.4.1 
## v readr   2.1.3      v forcats 0.5.2
```

```
## Warning: package 'tibble' was built under R version 4.2.1
```

```
## Warning: package 'tidyr' was built under R version 4.2.1
```

```
## Warning: package 'readr' was built under R version 4.2.1
```

```
## Warning: package 'dplyr' was built under R version 4.2.1
```

```
## Warning: package 'stringr' was built under R version 4.2.1
```

```
## Warning: package 'forcats' was built under R version 4.2.1
```

```
## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
## x dplyr::filter() masks stats::filter()
## x dplyr::lag()    masks stats::lag()
```

The message here stating that `dplyr::filter() masks stats::filter()` means that there are two functions/arguments of the same name from different packages. So, in order to ensure that we are using the function that we want, we can specify which package ot comes from using the name of the package followed by the `::` notation

A few ways to transform or subset your data

**dplyr::select**: to select certain columns  


```r
dat1 <- dat %>% 
  # the pipe operator %>% signifies consecutive steps that we're applying to our data
  dplyr::select(Season, Depth_m) # by name

head(dat1)
```

```
##   Season Depth_m
## 1   Fall      10
## 2   Fall      20
## 3   Fall      40
## 4   Fall      60
## 5   Fall      75
## 6   Fall      85
```

```r
dat1 <- dat %>% 
  dplyr::select(1,2) # by column number
```

**dplyr::filter**: to select rows containing a specific value 


```r
dat2 <- dat %>% 
  filter(Season == "Summer") 
# in R, to signify equivalence, we use double equal to signs ==
head(dat2)
```

```
##   Season Depth_m   O2_uM Add_data
## 1 Summer      10 216.667     TRUE
## 2 Summer      20 159.672    FALSE
## 3 Summer      40 141.778    FALSE
## 4 Summer      60  97.894     TRUE
## 5 Summer      75  44.978    FALSE
## 6 Summer      85  25.807    FALSE
```

**dplyr::mutate**: to make a new column  


```r
dat3 <- dat %>% 
  mutate(Depth_cm = Depth_m*100)
head(dat3)
```

```
##   Season Depth_m   O2_uM Add_data Depth_cm
## 1   Fall      10 203.533     TRUE     1000
## 2   Fall      20 183.787    FALSE     2000
## 3   Fall      40 130.579    FALSE     4000
## 4   Fall      60  91.115     TRUE     6000
## 5   Fall      75  69.828    FALSE     7500
## 6   Fall      85  26.972    FALSE     8500
```

### 5.3 Making Plots  

Making a very simple plot of O2 vs Depth using the `ggplot2` package:  


```r
dat %>% 
  ggplot(aes(x = O2_uM, y = Depth_m)) + #defining the x and y axis under the aes (aesthetics) parameter
  geom_point() # dot plot / scatter plot
```

![](TOG_Databinge_Intro_to_R_October2022_files/figure-html/plot-1.png)<!-- -->






