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

Who: Victor Yuan, [Trainee Omics Group (TOG)](https://bcchr-trainee-omics-group.github.io/)

What: Introduction to data wrangling in R

Where: BC Children's Hospital Research Institute

When: October 3rd, 2019

# Introduction 

## Who is this tutorial for?

- Those interested using R for data analysis
- Beginners in R, and those curious about `dplyr` and **tidyverse**

## Setup

- Please have the latest version of R and Rstudio installed. 
- Download the [r markdown file](workshops/2019-10-03_intro_to_data_manipulation/intro_to_data_manipulation.Rmd) for this workshop and open it in Rstudio.

- Install the following R packages now, if you haven't already:


```r
install.packages(c('gapminder', 'tidyverse'))
```

## Learning outcomes

- Familiarity with tidyverse syntax
- Understanding how to use dplyr to manipulate data in R
- "Piping"

# Brief intro to R and Rstudio

To learn data manipulation in R, some understanding of how to use R is required. A comprehensive lesson on R is out of the scope of this tutorial, but here I go over some of the basics that will be helpful to follow along the data wrangling portion of this workshop.

## Topics

- Running code from scripts / rmd. Line-by-line, run all, run code chunk
- Error messages
- `data.frame` and `tibble` objects

## Resources

- [R Swirl](https://swirlstats.com/) for interactive lesson on programming with R
- *R for data science* [tibbles chapter](http://r4ds.had.co.nz/tibbles.html) to learn more about tibble/data.frames

## Running code (4 min)

Two ways to run code in Rstudio:

1. Run code in console 
2. Run code from scripts (several options):

* Running code line-by-line (Windows: `ctrl` `+` `enter`, Mac: `cmd` `+` `enter`)
* Run code chunk (for *rmarkdown* files)
* Run all

I recommend using **r markdown file (recommended)** for all of your analysis scripts.

Demonstrate different ways to run code:


```r
1 + 2 + 3 + 4 + 5
```

```
## [1] 15
```

```r
"Code on multiple lines:"
```

```
## [1] "Code on multiple lines:"
```

```r
sum(1, 2,
    3,
    4, 5)
```

```
## [1] 15
```

In `.rmd` files, code needs to be in code chunks (insert code chunk shortcut is Windows: `ctrl` `+` `alt` `+` `i`, Mac: `cmd` `+` `alt` `+` `i`)

## Error messages (4 min)

When you encounter an error, this is what you should do:

1. Inspect your code for simple mistakes. Examples: are you missing a bracket or quotation mark? Did you mean to use `=` instead of `==`? Did you spell the name of the function or object wrong?
2. Make sure you are using the function correctly. Use the `?` command to pull up the help page.
3. Make sure your input is correct. Does the function expect a `data.frame` but you gave it a `vector`?
4. Copy and paste the error and the name of the function you are using into **Google**.


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

The columns of the trees `data.frame` object are individual `vector` objects. So trees has 3 columns/vectors that are each 31 elements long.

Some basic functions to help understand your `data.frame` objects.


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

We are going to use the "gapminder" dataset today, which is stored as a special type of `data.frame`, called a [`tibble`](https://tibble.tidyverse.org/index.html). 

- tibbles have a special printing output
- tibbles never have row names
- any function that works with `data.frame` objects (e.g. `str`, `summary`), will also work with `tibble` objects 
- use functions `as_tibble` and `as.data.frame` to convert between `tibble` and `data.frame`


```r
library(gapminder)
```

```
## Warning: package 'gapminder' was built under R version 3.6.1
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
## # ... with 1,694 more rows
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

# Intro to tidyverse

<!---The following chunk allows errors when knitting--->



## Topics

* "Base R" vs tidyverse
* dplyr basic functions
* Relational/comparison and logical operators

## Resources

Most of this material is inspired from the [STAT 545 `dplyr` lecture](https://stat545guidebook.netlify.com/intro-to-data-wrangling-part-i.html)

Other resources for `dplyr`:

* [Old STAT 545 `dplyr` tutorial](http://stat545.com/block009_dplyr-intro.html)
* [The `dplyr` vignette](https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html)

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

**Let's continue the tutorial in the raw r markdown file and work on the exercises.**

Before filling out the exercises,  we need to load in the `dplyr` and `gapminder` R packages.


```r
# load your packages here:
library(FILL_THIS_IN)
```

```
## Error in library(FILL_THIS_IN): there is no package called 'FILL_THIS_IN'
```

```r
library(FILL_THIS_IN)
```

```
## Error in library(FILL_THIS_IN): there is no package called 'FILL_THIS_IN'
```

### `select()` (8 min)

1. Make a data frame containing the columns `year`, `lifeExp`, `country` from the gapminder data, in that order.


```r
select(gapminder, FILL_THIS_IN)
```

```
## Error in select(gapminder, FILL_THIS_IN): could not find function "select"
```


2. Select all variables, from `country` to `lifeExp`.


```r
# This will work:
select(gapminder, country, continent, year, lifeExp)
```

```
## Error in select(gapminder, country, continent, year, lifeExp): could not find function "select"
```

```r
# Better way:
select(gapminder, FILL_THIS_IN)
```

```
## Error in select(gapminder, FILL_THIS_IN): could not find function "select"
```


3. Select all variables, except `lifeExp`.


```r
select(gapminder, FILL_THIS_IN)
```

```
## Error in select(gapminder, FILL_THIS_IN): could not find function "select"
```

4. Put `continent` first. Hint: use the `everything()` function.


```r
select(gapminder, FILL_THIS_IN, FILL_THIS_IN)
```

```
## Error in select(gapminder, FILL_THIS_IN, FILL_THIS_IN): could not find function "select"
```

### `arrange()` (8 min)

1. Order by year.


```r
arrange(gapminder, FILL_THIS_IN)
```

```
## Error in arrange(gapminder, FILL_THIS_IN): could not find function "arrange"
```

2. Order by year, in descending order.


```r
arrange(gapminder, FILL_THIS_IN)
```

```
## Error in arrange(gapminder, FILL_THIS_IN): could not find function "arrange"
```

3. Order by year, then by life expectancy.


```r
arrange(gapminder, FILL_THIS_IN, FILL_THIS_IN)
```

```
## Error in arrange(gapminder, FILL_THIS_IN, FILL_THIS_IN): could not find function "arrange"
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
# alter the below to include 2 "pipes"
gapminder %>%
  select(year, lifeExp, country) %>%
  arrange(year, lifeExp)
```

```
## Error in gapminder %>% select(year, lifeExp, country) %>% arrange(year, : could not find function "%>%"
```

### R Operators (4 min)

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


## `filter()` (6 min)

1. Only take data with population greater than 100 million.


```r
gapminder %>%
  filter(FILL_THIS_IN)
```

```
## Error in gapminder %>% filter(FILL_THIS_IN): could not find function "%>%"
```

2. Your turn: of those rows filtered from step 1., only take data from Asia.


```r
gapminder %>%
  filter(FILL_THIS_IN)
```

```
## Error in gapminder %>% filter(FILL_THIS_IN): could not find function "%>%"
```

3. Only take data from countries Brazil and China. 


```r
gapminder %>%
  filter(FILL_THIS_IN)
```

```
## Error in gapminder %>% filter(FILL_THIS_IN): could not find function "%>%"
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
  mutate(FILL_THIS_IN)
```

```
## Error in gapminder %>% mutate(FILL_THIS_IN): could not find function "%>%"
```

2. Make a new column named `GDP_bill`, that is GDP in billions.


```r
gapminder %>%
  mutate(FILL_THIS_IN)
```

```
## Error in gapminder %>% mutate(FILL_THIS_IN): could not find function "%>%"
```

Your turn: Make a new column called `cc` that pastes the country name followed by the continent, separated by a comma. (Hint: use the `paste` function with the `sep=", "` argument).

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
## Error in gapminder %>% summarize(mu = mean(lifeExp), sigma = sd(lifeExp)): could not find function "%>%"
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
## Error in gapminder %>% group_by(continent, year): could not find function "%>%"
```

The only thing different from a regular tibble is the indication of grouping variables above the tibble. This means that the tibble is recognized as having “chunks” defined by unique combinations of continent and year:

- Asia in 1952 is one chunk.
- Asia in 1957 is another chunk.
- Europe in 1952 is another chunk.
- etc…

Notice that the data frame isn’t rearranged by chunk!

Now that the tibble is grouped, operations that you do on a grouped tibble will be done independently within each chunk, as if no other chunks exist.

1. What is the mean and standard deviation for each year for every continent?


```r
gapminder %>% 
  group_by(continent, year) %>% 
  FILL_THIS_IN(mu    = FILL_THIS_IN,
            sigma = FILL_THIS_IN)
```

```
## Error in gapminder %>% group_by(continent, year) %>% FILL_THIS_IN(mu = FILL_THIS_IN, : could not find function "%>%"
```

2. In the gapminder dataset, how many rows are there for each continent? Hint: use the convenience function `dplyr::n()`


```r
# solution 1
gapminder %>%
  group_by(FILL_THIS_IN) %>%
  FILL_THIS_IN
```

```
## Error in gapminder %>% group_by(FILL_THIS_IN) %>% FILL_THIS_IN: could not find function "%>%"
```

```r
# solution 2: use dplyr::count()
```

2. (a) What's the minimum life expectancy for each continent and each year? (b) Arrange by min life expectancy.


```r
gapminder %>% 
  group_by(FILL_THIS_IN, FILL_THIS_IN) %>% 
  FILL_THIS_IN(min_life = min(FILL_THIS_IN)) %>%
  arrange(FILL_THIS_IN)
```

```
## Error in gapminder %>% group_by(FILL_THIS_IN, FILL_THIS_IN) %>% FILL_THIS_IN(min_life = min(FILL_THIS_IN)) %>% : could not find function "%>%"
```

3. Calculate the growth in population since the first year on record _for each country_. Here's another convenience function for you: `dplyr::first()`. 



# Bonus Exercises

If there's time remaining, we can try these more difficult exercises.

1. Take all countries in Europe that have a GDP per capita greater than 10000, and select all variables except `gdpPercap`.

2. Take the first three columns, and extract the variable names.

3. Of the `iris` data frame, take all columns that start with the word "Petal". 
    - Hint: take a look at the "Select helpers" documentation by running the following code: `?tidyselect::select_helpers`.
    
4. Convert the population to a number in billions.

5. Filter the rows of the iris dataset for Sepal.Length >= 4.6 and Petal.Width >= 0.5.

Exercises 3. and 5. are from [r-exercises](https://www.r-exercises.com/2017/10/19/dplyr-basic-functions-exercises/).
