---
title: "Introduction to data wrangling in R"
output:
  html_document:
    keep_md: true
    toc: true
    toc_depth: 3
    toc_float:
      collapsed: true
      smooth_scroll: true
    number_sections:  true
    theme: paper
    self_contained: yes
editor_options: 
  chunk_output_type: inline
---

Adapted from Tutorial by: Victor Yuan

Who: Almas Khan, [Trainee Omics Group (TOG)](https://bcchr-trainee-omics-group.github.io/)

What: Introduction to data wrangling in R

Where: BC Children's Hospital Research Institute

When: July 23,2020

# Introduction 

## Who is this tutorial for?

- Those interested using R for data analysis
- Beginners in R, and those curious about `dplyr` and **tidyverse**

## Setup

- Please have the latest version of R and Rstudio installed. 
- Download the [r markdown file](workshops/2020-07-23_data_wrangling_ak/Completed_instructor_data_wrangle_workshop.Rmd) for this workshop and open it in Rstudio.

- Install the following R packages now, if you haven't already:


```r
install.packages(c('gapminder', 'tidyverse'))
```


## Overall Learning outcomes

- Familiarity with tidyverse syntax
- Understanding how to use dplyr to manipulate data in R
- "Piping"

<!---The following chunk allows errors when knitting--->





To learn data manipulation in R, some understanding of how to use R is required. A comprehensive lesson on R is out of the scope of this tutorial, but here I go over some of the basics that will be helpful to follow along the data wrangling portion of this workshop.


## Resources

- [R Swirl](https://swirlstats.com/) for interactive lesson on programming with R
- *R for data science* [tibbles chapter](http://r4ds.had.co.nz/tibbles.html) to learn more about tibble/data.frames
- [Learnr package](https://rstudio.github.io/learnr/) another interactive lesson within R on programming


I recommend using **r markdown file (recommended)** for all of your analysis scripts. Think of it as a like a notebook for your code that includes comments and I find it's easier to keep track of your analysis. 

In `.rmd` files, code needs to be in code chunks (insert code chunk shortcut is Windows: `ctrl` `+` `alt` `+` `i`, Mac: `cmd` `+` `alt` `+` `i`) in order to run.


## Error messages (4 min)

As you work with R you will encounter error messages along the way, here are some steps I recommend you follow: 

1. Inspect your code for simple mistakes. Examples: are you missing a bracket or quotation mark? Did you mean to use `=` instead of `==`? Did you spell the name of the function or object wrong?
2. Are using the function correctly? Use the `?` command to pull up the help page. 
3. Make sure your input is correct. Does the function expect a `data.frame` but you gave it a `vector`?
4. Copy and paste the error and the name of the function you are using into **Google**.

# Understand the data in base R: Intro to data types and initial exploration

## Topics

- Basic data types
- `data.frame` and `tibble` objects
- Basic functions for exploring dataframes

Before we begin to manipulate data we should understand what type of data we have:

In R there are multiple types of data and I have listed a couple of basic ones below: 


```r
3.14 # This is a numeric 
```

```
## [1] 3.14
```

```r
"This is a character" # A character type 
```

```
## [1] "This is a character"
```

```r
"3.14" # This is also a character 
```

```
## [1] "3.14"
```

```r
TRUE # logical
```

```
## [1] TRUE
```

```r
FALSE # These are logicals
```

```
## [1] FALSE
```

These different types of data can be combined into objects like data frames or tibbles:

## Data frames and tibbles (5 min)

`data.frame` objects are by far the most common and useful way to work with your data in R.

A data.frame is basically just a table, it has a certain number of rows, and a certain number of columns.

"trees" is a built-in dataset in R available as a `data.frame`.


```r
trees
```

```
##    Girth Height Volume
## 1    8.3     70   10.3
## 2    8.6     65   10.3
## 3    8.8     63   10.2
## 4   10.5     72   16.4
## 5   10.7     81   18.8
## 6   10.8     83   19.7
## 7   11.0     66   15.6
## 8   11.0     75   18.2
## 9   11.1     80   22.6
## 10  11.2     75   19.9
## 11  11.3     79   24.2
## 12  11.4     76   21.0
## 13  11.4     76   21.4
## 14  11.7     69   21.3
## 15  12.0     75   19.1
## 16  12.9     74   22.2
## 17  12.9     85   33.8
## 18  13.3     86   27.4
## 19  13.7     71   25.7
## 20  13.8     64   24.9
## 21  14.0     78   34.5
## 22  14.2     80   31.7
## 23  14.5     74   36.3
## 24  16.0     72   38.3
## 25  16.3     77   42.6
## 26  17.3     81   55.4
## 27  17.5     82   55.7
## 28  17.9     80   58.3
## 29  18.0     80   51.5
## 30  18.0     80   51.0
## 31  20.6     87   77.0
```

The columns of the trees `data.frame` object are individual `vector` objects. So trees has 3 columns/vectors that are each 31 elements long. The dbl is just a type of numeric class of data. 

### Some basic functions to help understand your `data.frame` objects are:


```r
# number of rows
nrow(trees)
```

```
## [1] 31
```

```r
# number of columns
ncol(trees)
```

```
## [1] 3
```

```r
# row x columns
dim(trees)
```

```
## [1] 31  3
```

```r
# some basic info on the "structure" of the data.frame
str(trees)
```

```
## 'data.frame':	31 obs. of  3 variables:
##  $ Girth : num  8.3 8.6 8.8 10.5 10.7 10.8 11 11 11.1 11.2 ...
##  $ Height: num  70 65 63 72 81 83 66 75 80 75 ...
##  $ Volume: num  10.3 10.3 10.2 16.4 18.8 19.7 15.6 18.2 22.6 19.9 ...
```

```r
# calculates some summary statistics on each column
summary(trees)
```

```
##      Girth           Height       Volume     
##  Min.   : 8.30   Min.   :63   Min.   :10.20  
##  1st Qu.:11.05   1st Qu.:72   1st Qu.:19.40  
##  Median :12.90   Median :76   Median :24.20  
##  Mean   :13.25   Mean   :76   Mean   :30.17  
##  3rd Qu.:15.25   3rd Qu.:80   3rd Qu.:37.30  
##  Max.   :20.60   Max.   :87   Max.   :77.00
```

```r
# print first 6 rows
head(trees)
```

```
##   Girth Height Volume
## 1   8.3     70   10.3
## 2   8.6     65   10.3
## 3   8.8     63   10.2
## 4  10.5     72   16.4
## 5  10.7     81   18.8
## 6  10.8     83   19.7
```

```r
# print last 6 rows
tail(trees)
```

```
##    Girth Height Volume
## 26  17.3     81   55.4
## 27  17.5     82   55.7
## 28  17.9     80   58.3
## 29  18.0     80   51.5
## 30  18.0     80   51.0
## 31  20.6     87   77.0
```

### Tibbles

We are going to use the "gapminder" dataset today, which is stored as a special type of `data.frame`, called a [`tibble`](https://tibble.tidyverse.org/index.html). 

- tibbles have a special printing output
- tibbles never have row names
- any function that works with `data.frame` objects (e.g. `str`, `summary`), will also work with `tibble` objects 
- use functions `as_tibble` and `as.data.frame` to convert between `tibble` and `data.frame`


```r
library(gapminder) #load gapminder
```


```r
# tibble
gapminder
```

```
## # A tibble: 1,704 x 6
##    country     continent  year lifeExp      pop gdpPercap
##    <fct>       <fct>     <int>   <dbl>    <int>     <dbl>
##  1 Afghanistan Asia       1952    28.8  8425333      779.
##  2 Afghanistan Asia       1957    30.3  9240934      821.
##  3 Afghanistan Asia       1962    32.0 10267083      853.
##  4 Afghanistan Asia       1967    34.0 11537966      836.
##  5 Afghanistan Asia       1972    36.1 13079460      740.
##  6 Afghanistan Asia       1977    38.4 14880372      786.
##  7 Afghanistan Asia       1982    39.9 12881816      978.
##  8 Afghanistan Asia       1987    40.8 13867957      852.
##  9 Afghanistan Asia       1992    41.7 16317921      649.
## 10 Afghanistan Asia       1997    41.8 22227415      635.
## # … with 1,694 more rows
```

```r
# data.frame 
head(as.data.frame(gapminder), n = 50)
```

```
##        country continent year lifeExp      pop gdpPercap
## 1  Afghanistan      Asia 1952  28.801  8425333  779.4453
## 2  Afghanistan      Asia 1957  30.332  9240934  820.8530
## 3  Afghanistan      Asia 1962  31.997 10267083  853.1007
## 4  Afghanistan      Asia 1967  34.020 11537966  836.1971
## 5  Afghanistan      Asia 1972  36.088 13079460  739.9811
## 6  Afghanistan      Asia 1977  38.438 14880372  786.1134
## 7  Afghanistan      Asia 1982  39.854 12881816  978.0114
## 8  Afghanistan      Asia 1987  40.822 13867957  852.3959
## 9  Afghanistan      Asia 1992  41.674 16317921  649.3414
## 10 Afghanistan      Asia 1997  41.763 22227415  635.3414
## 11 Afghanistan      Asia 2002  42.129 25268405  726.7341
## 12 Afghanistan      Asia 2007  43.828 31889923  974.5803
## 13     Albania    Europe 1952  55.230  1282697 1601.0561
## 14     Albania    Europe 1957  59.280  1476505 1942.2842
## 15     Albania    Europe 1962  64.820  1728137 2312.8890
## 16     Albania    Europe 1967  66.220  1984060 2760.1969
## 17     Albania    Europe 1972  67.690  2263554 3313.4222
## 18     Albania    Europe 1977  68.930  2509048 3533.0039
## 19     Albania    Europe 1982  70.420  2780097 3630.8807
## 20     Albania    Europe 1987  72.000  3075321 3738.9327
## 21     Albania    Europe 1992  71.581  3326498 2497.4379
## 22     Albania    Europe 1997  72.950  3428038 3193.0546
## 23     Albania    Europe 2002  75.651  3508512 4604.2117
## 24     Albania    Europe 2007  76.423  3600523 5937.0295
## 25     Algeria    Africa 1952  43.077  9279525 2449.0082
## 26     Algeria    Africa 1957  45.685 10270856 3013.9760
## 27     Algeria    Africa 1962  48.303 11000948 2550.8169
## 28     Algeria    Africa 1967  51.407 12760499 3246.9918
## 29     Algeria    Africa 1972  54.518 14760787 4182.6638
## 30     Algeria    Africa 1977  58.014 17152804 4910.4168
## 31     Algeria    Africa 1982  61.368 20033753 5745.1602
## 32     Algeria    Africa 1987  65.799 23254956 5681.3585
## 33     Algeria    Africa 1992  67.744 26298373 5023.2166
## 34     Algeria    Africa 1997  69.152 29072015 4797.2951
## 35     Algeria    Africa 2002  70.994 31287142 5288.0404
## 36     Algeria    Africa 2007  72.301 33333216 6223.3675
## 37      Angola    Africa 1952  30.015  4232095 3520.6103
## 38      Angola    Africa 1957  31.999  4561361 3827.9405
## 39      Angola    Africa 1962  34.000  4826015 4269.2767
## 40      Angola    Africa 1967  35.985  5247469 5522.7764
## 41      Angola    Africa 1972  37.928  5894858 5473.2880
## 42      Angola    Africa 1977  39.483  6162675 3008.6474
## 43      Angola    Africa 1982  39.942  7016384 2756.9537
## 44      Angola    Africa 1987  39.906  7874230 2430.2083
## 45      Angola    Africa 1992  40.647  8735988 2627.8457
## 46      Angola    Africa 1997  40.963  9875024 2277.1409
## 47      Angola    Africa 2002  41.003 10866106 2773.2873
## 48      Angola    Africa 2007  42.731 12420476 4797.2313
## 49   Argentina  Americas 1952  62.485 17876956 5911.3151
## 50   Argentina  Americas 1957  64.399 19610538 6856.8562
```

Note that I printed only the first 50 rows of the `data.frame` version of gapminder. By default, printing `tibble` objects will display the top 10 rows to avoid over-printing.

# Let's wrangle:  Intro to tidyverse

Now that we've learned some basics about data types in R. Let's wrangle some data using the tidyverse package. 

## Topics

* Relational/comparison and logical operators
* "Base R" vs tidyverse
* dplyr basic functions


## Resources

Most of this material is inspired from the [STAT 545 `dplyr` lecture](https://stat545guidebook.netlify.com/intro-to-data-wrangling-part-i.html)

Other resources for `dplyr`:

* [Old STAT 545 `dplyr` tutorial](http://stat545.com/block009_dplyr-intro.html)
* [The `dplyr` vignette](https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html)

### R Operators (4 min)

Before diving into the tidyverse, let's take a look at a list of built-in operators in R that can help in wrangling our data and using within tidyverse functions. 

See knitted html for pretty tables.

**Arithmetic** operators allow us to carry out mathematical operations:

| Operator | Description |
|------|:---------|
| + | Add |
| - | Subtract |
| * | Multiply |
| / | Divide |
| ^ | Exponent |
| %% | Modulus (remainder from division) |

**Relational** operators allow us to compare values:

| Operator | Description |
|------|:---------|
| < | Less than |
| > | Greater than |
| <= | Less than or equal to |
| >= | Greater than or equal to |
| == | Equal to |
| != | Not equal to |

* Arithmetic and relational operators work on vectors.

**Logical** operators allow us to carry out boolean operations:

| Operator | Description |
|---|:---|
| ! | Not |
| \| | Or (element_wise) |
| & | And (element-wise) |
| \|\| | Or |
| && | And |

## Base-R vs tidyverse

The [tidyverse](https://www.tidyverse.org/) is a collection of packages that were designed under a specific philosophy of making data analysis in R easy, reproducible, and consistent. The tidyverse packages:

* has an easy-to-read syntax
* can be much quicker (to write)  
* are comprehensive in covering a wide variety of applications (e.g. [`ggplot2`](https://ggplot2.tidyverse.org/) for plotting, [`dplyr`](https://dplyr.tidyverse.org/) for data manipulation, ..., etc.)

Most base-r operations can be replaced with tidyverse code. However, in genomics research, this is less true because many genomics analysis workflows rely on [**Bioconductor**](https://www.bioconductor.org/).

[Bioconductor](https://www.bioconductor.org/) is an online R package repository specifically for analyzing and storing genomic data. If you are in genomics, it is highly likely you will use these packages. Bioconductor is basically all based in base-R code, so learning *some* base-R is going to be a requirement for any genomics analysis.


## Tidyverse syntax (2 min)

- Almost all functions from a **tidyverse** package will accept a `tibble` or `data.frame` as the first argument.

```
function(a_tibble, operations_involving_columns)
```

- **tidyverse** functions tend to be [named](https://principles.tidyverse.org/function-names.html) as verbs. (e.g. `select()`, `pivot_longer()`, `str_extract()`)

- More information on the design principles behind tidyverse packages can be found [here](https://principles.tidyverse.org/index.html).

# Basic `dplyr` functions

We're going to be working with the gapminder dataframe we loaded before. Additionally we need to we need to load in the `dplyr` R package. This is a package within the tidyverse that stands for data frame applyer and contains functions to help wrangle data in data.frames or tibbles. 



```r
# load your package here:
library(dplyr)
```

We're going to go over these dplyr functions:

1. select 
2. arrange
3. filter 
4. mutate
5. grouped operations in dplyr:
    a. group_by
    b. summarize

### `select()` (8 min)

`select()` allows you to subset to the columns(or variables) that you specify. 

1. Make a data frame containing the columns `year`, `lifeExp`, `country` from the gapminder data, in that order.


```r
select(gapminder, year,lifeExp,country)
```

```
## # A tibble: 1,704 x 3
##     year lifeExp country    
##    <int>   <dbl> <fct>      
##  1  1952    28.8 Afghanistan
##  2  1957    30.3 Afghanistan
##  3  1962    32.0 Afghanistan
##  4  1967    34.0 Afghanistan
##  5  1972    36.1 Afghanistan
##  6  1977    38.4 Afghanistan
##  7  1982    39.9 Afghanistan
##  8  1987    40.8 Afghanistan
##  9  1992    41.7 Afghanistan
## 10  1997    41.8 Afghanistan
## # … with 1,694 more rows
```


2. Select all variables, from `country` to `lifeExp`.


```r
# This will work:
select(gapminder, country, continent, year, lifeExp)
```

```
## # A tibble: 1,704 x 4
##    country     continent  year lifeExp
##    <fct>       <fct>     <int>   <dbl>
##  1 Afghanistan Asia       1952    28.8
##  2 Afghanistan Asia       1957    30.3
##  3 Afghanistan Asia       1962    32.0
##  4 Afghanistan Asia       1967    34.0
##  5 Afghanistan Asia       1972    36.1
##  6 Afghanistan Asia       1977    38.4
##  7 Afghanistan Asia       1982    39.9
##  8 Afghanistan Asia       1987    40.8
##  9 Afghanistan Asia       1992    41.7
## 10 Afghanistan Asia       1997    41.8
## # … with 1,694 more rows
```

```r
# Better way:
select(gapminder, country:lifeExp)
```

```
## # A tibble: 1,704 x 4
##    country     continent  year lifeExp
##    <fct>       <fct>     <int>   <dbl>
##  1 Afghanistan Asia       1952    28.8
##  2 Afghanistan Asia       1957    30.3
##  3 Afghanistan Asia       1962    32.0
##  4 Afghanistan Asia       1967    34.0
##  5 Afghanistan Asia       1972    36.1
##  6 Afghanistan Asia       1977    38.4
##  7 Afghanistan Asia       1982    39.9
##  8 Afghanistan Asia       1987    40.8
##  9 Afghanistan Asia       1992    41.7
## 10 Afghanistan Asia       1997    41.8
## # … with 1,694 more rows
```


3. Select all variables, except `lifeExp`.


```r
select(gapminder, country:year)
```

```
## # A tibble: 1,704 x 3
##    country     continent  year
##    <fct>       <fct>     <int>
##  1 Afghanistan Asia       1952
##  2 Afghanistan Asia       1957
##  3 Afghanistan Asia       1962
##  4 Afghanistan Asia       1967
##  5 Afghanistan Asia       1972
##  6 Afghanistan Asia       1977
##  7 Afghanistan Asia       1982
##  8 Afghanistan Asia       1987
##  9 Afghanistan Asia       1992
## 10 Afghanistan Asia       1997
## # … with 1,694 more rows
```

4. Put `continent` first. Hint: use the `everything()` function.


```r
select(gapminder, continent,everything())
```

```
## # A tibble: 1,704 x 6
##    continent country      year lifeExp      pop gdpPercap
##    <fct>     <fct>       <int>   <dbl>    <int>     <dbl>
##  1 Asia      Afghanistan  1952    28.8  8425333      779.
##  2 Asia      Afghanistan  1957    30.3  9240934      821.
##  3 Asia      Afghanistan  1962    32.0 10267083      853.
##  4 Asia      Afghanistan  1967    34.0 11537966      836.
##  5 Asia      Afghanistan  1972    36.1 13079460      740.
##  6 Asia      Afghanistan  1977    38.4 14880372      786.
##  7 Asia      Afghanistan  1982    39.9 12881816      978.
##  8 Asia      Afghanistan  1987    40.8 13867957      852.
##  9 Asia      Afghanistan  1992    41.7 16317921      649.
## 10 Asia      Afghanistan  1997    41.8 22227415      635.
## # … with 1,694 more rows
```

### `arrange()` (8 min)

1. Order by year.


```r
arrange(gapminder, year)
```

```
## # A tibble: 1,704 x 6
##    country     continent  year lifeExp      pop gdpPercap
##    <fct>       <fct>     <int>   <dbl>    <int>     <dbl>
##  1 Afghanistan Asia       1952    28.8  8425333      779.
##  2 Albania     Europe     1952    55.2  1282697     1601.
##  3 Algeria     Africa     1952    43.1  9279525     2449.
##  4 Angola      Africa     1952    30.0  4232095     3521.
##  5 Argentina   Americas   1952    62.5 17876956     5911.
##  6 Australia   Oceania    1952    69.1  8691212    10040.
##  7 Austria     Europe     1952    66.8  6927772     6137.
##  8 Bahrain     Asia       1952    50.9   120447     9867.
##  9 Bangladesh  Asia       1952    37.5 46886859      684.
## 10 Belgium     Europe     1952    68    8730405     8343.
## # … with 1,694 more rows
```

2. Order by year, in descending order.


```r
arrange(gapminder, desc(year))
```

```
## # A tibble: 1,704 x 6
##    country     continent  year lifeExp       pop gdpPercap
##    <fct>       <fct>     <int>   <dbl>     <int>     <dbl>
##  1 Afghanistan Asia       2007    43.8  31889923      975.
##  2 Albania     Europe     2007    76.4   3600523     5937.
##  3 Algeria     Africa     2007    72.3  33333216     6223.
##  4 Angola      Africa     2007    42.7  12420476     4797.
##  5 Argentina   Americas   2007    75.3  40301927    12779.
##  6 Australia   Oceania    2007    81.2  20434176    34435.
##  7 Austria     Europe     2007    79.8   8199783    36126.
##  8 Bahrain     Asia       2007    75.6    708573    29796.
##  9 Bangladesh  Asia       2007    64.1 150448339     1391.
## 10 Belgium     Europe     2007    79.4  10392226    33693.
## # … with 1,694 more rows
```

3. Order by year, then by life expectancy.


```r
arrange(gapminder,year, lifeExp)
```

```
## # A tibble: 1,704 x 6
##    country       continent  year lifeExp     pop gdpPercap
##    <fct>         <fct>     <int>   <dbl>   <int>     <dbl>
##  1 Afghanistan   Asia       1952    28.8 8425333      779.
##  2 Gambia        Africa     1952    30    284320      485.
##  3 Angola        Africa     1952    30.0 4232095     3521.
##  4 Sierra Leone  Africa     1952    30.3 2143249      880.
##  5 Mozambique    Africa     1952    31.3 6446316      469.
##  6 Burkina Faso  Africa     1952    32.0 4469979      543.
##  7 Guinea-Bissau Africa     1952    32.5  580653      300.
##  8 Yemen, Rep.   Asia       1952    32.5 4963829      782.
##  9 Somalia       Africa     1952    33.0 2526994     1136.
## 10 Guinea        Africa     1952    33.6 2664249      510.
## # … with 1,694 more rows
```

## Piping, `%>%` (8 min)

*Piping* refers to using the `%>%` operator to write nested function calls in a more readable manner.

- Takes an output as the input for the first argument of the next function.

- Think of `%>%` as the word "then"!

**Demonstration:** Here I want to combine `select()` Task 1 with `arrange()` Task 3.

This is how I could do it by *nesting* the two function calls, **without piping**:


```r
# Nesting function calls can be hard to read
arrange(select(gapminder, year, lifeExp, country), year, lifeExp)
```

Now using **with piping**:


```r
# alter the above to include 2 "pipes"
gapminder %>%
  select(year, lifeExp, country) %>%
  arrange(year, lifeExp)
```

```
## # A tibble: 1,704 x 3
##     year lifeExp country      
##    <int>   <dbl> <fct>        
##  1  1952    28.8 Afghanistan  
##  2  1952    30   Gambia       
##  3  1952    30.0 Angola       
##  4  1952    30.3 Sierra Leone 
##  5  1952    31.3 Mozambique   
##  6  1952    32.0 Burkina Faso 
##  7  1952    32.5 Guinea-Bissau
##  8  1952    32.5 Yemen, Rep.  
##  9  1952    33.0 Somalia      
## 10  1952    33.6 Guinea       
## # … with 1,694 more rows
```


## `filter()` (6 min)

Use `filter()` to subset to rows within your data where the condition you specify is TRUE. 

1. Only take data with population greater than 100 million.


```r
gapminder %>%
  filter(pop>100000000)
```

```
## # A tibble: 77 x 6
##    country    continent  year lifeExp       pop gdpPercap
##    <fct>      <fct>     <int>   <dbl>     <int>     <dbl>
##  1 Bangladesh Asia       1987    52.8 103764241      752.
##  2 Bangladesh Asia       1992    56.0 113704579      838.
##  3 Bangladesh Asia       1997    59.4 123315288      973.
##  4 Bangladesh Asia       2002    62.0 135656790     1136.
##  5 Bangladesh Asia       2007    64.1 150448339     1391.
##  6 Brazil     Americas   1972    59.5 100840058     4986.
##  7 Brazil     Americas   1977    61.5 114313951     6660.
##  8 Brazil     Americas   1982    63.3 128962939     7031.
##  9 Brazil     Americas   1987    65.2 142938076     7807.
## 10 Brazil     Americas   1992    67.1 155975974     6950.
## # … with 67 more rows
```

2. Your turn: of those rows filtered from step 1., only take data from Asia.


```r
gapminder %>%
  filter(pop>100000000 & continent=='Asia')
```

```
## # A tibble: 52 x 6
##    country    continent  year lifeExp       pop gdpPercap
##    <fct>      <fct>     <int>   <dbl>     <int>     <dbl>
##  1 Bangladesh Asia       1987    52.8 103764241      752.
##  2 Bangladesh Asia       1992    56.0 113704579      838.
##  3 Bangladesh Asia       1997    59.4 123315288      973.
##  4 Bangladesh Asia       2002    62.0 135656790     1136.
##  5 Bangladesh Asia       2007    64.1 150448339     1391.
##  6 China      Asia       1952    44   556263527      400.
##  7 China      Asia       1957    50.5 637408000      576.
##  8 China      Asia       1962    44.5 665770000      488.
##  9 China      Asia       1967    58.4 754550000      613.
## 10 China      Asia       1972    63.1 862030000      677.
## # … with 42 more rows
```

3. Only take data from countries Brazil and China. 


```r
gapminder %>%
  filter(pop>100000000,country == "Brazil"| country == "China")
```

```
## # A tibble: 20 x 6
##    country continent  year lifeExp        pop gdpPercap
##    <fct>   <fct>     <int>   <dbl>      <int>     <dbl>
##  1 Brazil  Americas   1972    59.5  100840058     4986.
##  2 Brazil  Americas   1977    61.5  114313951     6660.
##  3 Brazil  Americas   1982    63.3  128962939     7031.
##  4 Brazil  Americas   1987    65.2  142938076     7807.
##  5 Brazil  Americas   1992    67.1  155975974     6950.
##  6 Brazil  Americas   1997    69.4  168546719     7958.
##  7 Brazil  Americas   2002    71.0  179914212     8131.
##  8 Brazil  Americas   2007    72.4  190010647     9066.
##  9 China   Asia       1952    44    556263527      400.
## 10 China   Asia       1957    50.5  637408000      576.
## 11 China   Asia       1962    44.5  665770000      488.
## 12 China   Asia       1967    58.4  754550000      613.
## 13 China   Asia       1972    63.1  862030000      677.
## 14 China   Asia       1977    64.0  943455000      741.
## 15 China   Asia       1982    65.5 1000281000      962.
## 16 China   Asia       1987    67.3 1084035000     1379.
## 17 China   Asia       1992    68.7 1164970000     1656.
## 18 China   Asia       1997    70.4 1230075000     2289.
## 19 China   Asia       2002    72.0 1280400000     3119.
## 20 China   Asia       2007    73.0 1318683096     4959.
```

## `mutate()` (8 min)

The `mutate()` function _creates_ new columns in the tibble by transforming other variables. Like `select()`, `filter()`, and `arrange()`, the `mutate()` function also takes a tibble as its first argument, and returns a tibble. 

The general syntax is:

```
mutate(tibble, NEW_COLUMN_NAME = CALCULATION)
```

1. Make a new column named `GDP` that equals to multiplying GPD per capita with population.


```r
gapminder %>%
  mutate(GDP=gdpPercap*pop)
```

```
## # A tibble: 1,704 x 7
##    country     continent  year lifeExp      pop gdpPercap          GDP
##    <fct>       <fct>     <int>   <dbl>    <int>     <dbl>        <dbl>
##  1 Afghanistan Asia       1952    28.8  8425333      779.  6567086330.
##  2 Afghanistan Asia       1957    30.3  9240934      821.  7585448670.
##  3 Afghanistan Asia       1962    32.0 10267083      853.  8758855797.
##  4 Afghanistan Asia       1967    34.0 11537966      836.  9648014150.
##  5 Afghanistan Asia       1972    36.1 13079460      740.  9678553274.
##  6 Afghanistan Asia       1977    38.4 14880372      786. 11697659231.
##  7 Afghanistan Asia       1982    39.9 12881816      978. 12598563401.
##  8 Afghanistan Asia       1987    40.8 13867957      852. 11820990309.
##  9 Afghanistan Asia       1992    41.7 16317921      649. 10595901589.
## 10 Afghanistan Asia       1997    41.8 22227415      635. 14121995875.
## # … with 1,694 more rows
```

2. Make a new column named `GDP_bill`, that is GDP in billions.


```r
gapminder %>%
 mutate(gdpBill = round(gdpPercap*pop/1E9,digits=2))
```

```
## # A tibble: 1,704 x 7
##    country     continent  year lifeExp      pop gdpPercap gdpBill
##    <fct>       <fct>     <int>   <dbl>    <int>     <dbl>   <dbl>
##  1 Afghanistan Asia       1952    28.8  8425333      779.    6.57
##  2 Afghanistan Asia       1957    30.3  9240934      821.    7.59
##  3 Afghanistan Asia       1962    32.0 10267083      853.    8.76
##  4 Afghanistan Asia       1967    34.0 11537966      836.    9.65
##  5 Afghanistan Asia       1972    36.1 13079460      740.    9.68
##  6 Afghanistan Asia       1977    38.4 14880372      786.   11.7 
##  7 Afghanistan Asia       1982    39.9 12881816      978.   12.6 
##  8 Afghanistan Asia       1987    40.8 13867957      852.   11.8 
##  9 Afghanistan Asia       1992    41.7 16317921      649.   10.6 
## 10 Afghanistan Asia       1997    41.8 22227415      635.   14.1 
## # … with 1,694 more rows
```

Your turn: Make a new column called `cc` that pastes the country name followed by the continent, separated by a comma. (Hint: use the `paste` function with the `sep=", "` argument).


```r
gapminder %>%
  mutate(cc = paste(country,continent,sep = ","))
```

```
## # A tibble: 1,704 x 7
##    country     continent  year lifeExp      pop gdpPercap cc              
##    <fct>       <fct>     <int>   <dbl>    <int>     <dbl> <chr>           
##  1 Afghanistan Asia       1952    28.8  8425333      779. Afghanistan,Asia
##  2 Afghanistan Asia       1957    30.3  9240934      821. Afghanistan,Asia
##  3 Afghanistan Asia       1962    32.0 10267083      853. Afghanistan,Asia
##  4 Afghanistan Asia       1967    34.0 11537966      836. Afghanistan,Asia
##  5 Afghanistan Asia       1972    36.1 13079460      740. Afghanistan,Asia
##  6 Afghanistan Asia       1977    38.4 14880372      786. Afghanistan,Asia
##  7 Afghanistan Asia       1982    39.9 12881816      978. Afghanistan,Asia
##  8 Afghanistan Asia       1987    40.8 13867957      852. Afghanistan,Asia
##  9 Afghanistan Asia       1992    41.7 16317921      649. Afghanistan,Asia
## 10 Afghanistan Asia       1997    41.8 22227415      635. Afghanistan,Asia
## # … with 1,694 more rows
```

# `dplyr`: grouped operations

*Most of this material was taken/modified from the [STAT545 lecture: Intro to data wrangling, Part II](https://stat545guidebook.netlify.com/intro-to-data-wrangling-part-ii.html)*

## Resources

- [STAT545 lecture: Intro to data wrangling, Part II](https://stat545guidebook.netlify.com/intro-to-data-wrangling-part-ii.html)
- [r4ds: transform](http://r4ds.had.co.nz/transform.html)
- [Intro to `dplyr` vignette](https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html)

## `summarize()` (3 min)

Like `mutate()`, the `summarize()` function also creates new columns, but the calculations that make the new columns must reduce down to a single number.

For example, let’s compute the mean and standard deviation of life expectancy in the gapminder data set:


```r
gapminder %>% 
  summarize(mu    = mean(lifeExp),
            sigma = sd(lifeExp))
```

```
## # A tibble: 1 x 2
##      mu sigma
##   <dbl> <dbl>
## 1  59.5  12.9
```

Notice that all other columns were dropped. This is necessary, because there’s no obvious way to compress the other columns down to a single row. This is unlike `mutate(`), which keeps all columns.

As it is, this is hardly useful. But that’s outside of the context of *grouping*, coming up next...

## `group_by()` (15 min)

* `group_by()` allows you to apply functions to separate chunks of your data frame, where the chunks are defined by a specified grouping variable.

Let’s group the gapminder dataset by continent and year:


```r
gapminder %>% 
  group_by(continent, year)
```

```
## # A tibble: 1,704 x 6
## # Groups:   continent, year [60]
##    country     continent  year lifeExp      pop gdpPercap
##    <fct>       <fct>     <int>   <dbl>    <int>     <dbl>
##  1 Afghanistan Asia       1952    28.8  8425333      779.
##  2 Afghanistan Asia       1957    30.3  9240934      821.
##  3 Afghanistan Asia       1962    32.0 10267083      853.
##  4 Afghanistan Asia       1967    34.0 11537966      836.
##  5 Afghanistan Asia       1972    36.1 13079460      740.
##  6 Afghanistan Asia       1977    38.4 14880372      786.
##  7 Afghanistan Asia       1982    39.9 12881816      978.
##  8 Afghanistan Asia       1987    40.8 13867957      852.
##  9 Afghanistan Asia       1992    41.7 16317921      649.
## 10 Afghanistan Asia       1997    41.8 22227415      635.
## # … with 1,694 more rows
```

The only thing different from a regular tibble is the indication of grouping variables above the tibble. This means that the tibble is recognized as having “chunks” defined by unique combinations of continent and year:

- Asia in 1952 is one chunk.
- Asia in 1957 is another chunk.
- Europe in 1952 is another chunk.
- etc…

Notice that the data frame isn’t rearranged by chunk!

Now that the tibble is grouped, operations that you do on a grouped tibble will be done independently within each chunk, as if no other chunks exist.

1. What is the mean and standard deviation of life expectancy for each year for every continent?


```r
gapminder %>% 
  group_by(continent, year) %>% 
  summarize(mu    = mean(lifeExp),
            sigma = sd(lifeExp))
```

```
## # A tibble: 60 x 4
## # Groups:   continent [5]
##    continent  year    mu sigma
##    <fct>     <int> <dbl> <dbl>
##  1 Africa     1952  39.1  5.15
##  2 Africa     1957  41.3  5.62
##  3 Africa     1962  43.3  5.88
##  4 Africa     1967  45.3  6.08
##  5 Africa     1972  47.5  6.42
##  6 Africa     1977  49.6  6.81
##  7 Africa     1982  51.6  7.38
##  8 Africa     1987  53.3  7.86
##  9 Africa     1992  53.6  9.46
## 10 Africa     1997  53.6  9.10
## # … with 50 more rows
```

2. In the gapminder dataset, how many rows are there for each continent? Hint: use the convenience function `dplyr::n()`


```r
# solution 1
gapminder %>%
  group_by(continent) %>%
 summarise(n())
```

```
## # A tibble: 5 x 2
##   continent `n()`
##   <fct>     <int>
## 1 Africa      624
## 2 Americas    300
## 3 Asia        396
## 4 Europe      360
## 5 Oceania      24
```

```r
# solution 2: use dplyr::count()
gapminder %>%
  count(continent)
```

```
## # A tibble: 5 x 2
##   continent     n
##   <fct>     <int>
## 1 Africa      624
## 2 Americas    300
## 3 Asia        396
## 4 Europe      360
## 5 Oceania      24
```

2. (a) What's the minimum life expectancy for each continent and each year? (b) Arrange by min life expectancy.


```r
gapminder %>% 
  group_by(continent,year) %>% 
 summarize(min_life = min(lifeExp)) %>% 
  arrange (min_life)
```

```
## # A tibble: 60 x 3
## # Groups:   continent [5]
##    continent  year min_life
##    <fct>     <int>    <dbl>
##  1 Africa     1992     23.6
##  2 Asia       1952     28.8
##  3 Africa     1952     30  
##  4 Asia       1957     30.3
##  5 Asia       1977     31.2
##  6 Africa     1957     31.6
##  7 Asia       1962     32.0
##  8 Africa     1962     32.8
##  9 Asia       1967     34.0
## 10 Africa     1967     34.1
## # … with 50 more rows
```

3. Calculate the growth in population since the first year on record _for each country_. Here's another convenience function for you: `dplyr::first()`. 


```r
gapminder %>% 
group_by(country) %>% 
arrange(year) %>% 
mutate(rel_growth = pop-first(pop)) 
```

```
## # A tibble: 1,704 x 7
## # Groups:   country [142]
##    country     continent  year lifeExp      pop gdpPercap rel_growth
##    <fct>       <fct>     <int>   <dbl>    <int>     <dbl>      <int>
##  1 Afghanistan Asia       1952    28.8  8425333      779.          0
##  2 Albania     Europe     1952    55.2  1282697     1601.          0
##  3 Algeria     Africa     1952    43.1  9279525     2449.          0
##  4 Angola      Africa     1952    30.0  4232095     3521.          0
##  5 Argentina   Americas   1952    62.5 17876956     5911.          0
##  6 Australia   Oceania    1952    69.1  8691212    10040.          0
##  7 Austria     Europe     1952    66.8  6927772     6137.          0
##  8 Bahrain     Asia       1952    50.9   120447     9867.          0
##  9 Bangladesh  Asia       1952    37.5 46886859      684.          0
## 10 Belgium     Europe     1952    68    8730405     8343.          0
## # … with 1,694 more rows
```

