
# Introduction
For a nicer version of this tutorial go to [http://rpubs.com/NCI-biostats/microbiome-fixed-reference](http://rpubs.com/NCI-biostats/microbiome-fixed-reference).

In this tutorial we describe and show examples of R code implementing methods to use standard microbiome reference groups to simplify beta-diversity analyses and faclitate independent validation as described in **Using standard microbiome reference groups to simplify beta-diversity analyses and facilitate independent validation** by Marlena Maziarz, Ruth M. Pfeiffer, Yunhu Wan and Mitchell H. Gail. The code, test files, HMP reference sets and any other necessary files are on [https://github.com/NCI-biostats/microbiome-fixed-reference](https://github.com/NCI-biostats/microbiome-fixed-reference).

This tutorial describes the use of two functions: **HMPdistance** and **GHotelling**. HMPdistance calculates distance to the HMP stool and nasal reference sets, GHotelling is a generalized Hotelling test, which accounts for correlation between samples from the same individual.




# Calculating distance to HMP reference sets
## Load required packages and source the needed code from GitHub

```r
require(RCurl)  # for sourcing code from GitHub
require(repmis) # for sourcing rdata from GitHub 
```





```r
script <- getURL("https://raw.githubusercontent.com/NCI-biostats/microbiome-fixed-reference/master/microbiome-fixed-reference-functions.r", ssl.verifypeer = FALSE)
eval(parse(text = script))
```

## **HMPdistance** function definition

```r
# details below, in short:
# input
# tax.level = one of 'l2.phylum', 'l3.class', 'l4.order', 'l5.family', 'l6.genus'
# d.new.filename = your realative abundace data in csv or rdata format (it should be a local file, include the path in the filename, unless in current directory)
# d.new.ix.col.not.rel.abu = index of columns that are not relative abundnce, such as metadata, these will be excluded
# measure = 'bc' for Bray-Curtis distance, 'corr' for 1-Pearson correlation
#
# output
# n x 2 matrix with distances to HMP stool and HMP nasal references
HMPdistance <- function(tax.level = NULL,                 
                        d.new.filename = NULL,            
                        d.new.ix.col.not.rel.abu = NULL,  
                        measure = 'bc',                   
                        print.details = F)
```
This function computes the mean distance dstool of a microbiome sample (or samples) to a reference set of 92 stool samples from HMP (Human Microbiome Project) and the mean distance dnasal to a reference set of 74 HMP nasal samples.

### Input to **HMPdistance**
**tax.level** = defines the taxonomic level for which the distances are to be computed, should be one of: 'l2.phylum', 'l3.class', 'l4.order', 'l5.family', 'l6.genus'.

**d.new.filename** = filename of your csv or rdata (including path, ready to load) dataset, samples/subjects in rows, relative abundance (+ possibly metadata) in columns (more details below). 

**d.new.ix.col.not.rel.abu** = index of columns that are not relative abunance columns (ie. index of metadata columnsm if any), default = NULL.

**measure** = ‘bc’ for Bray-Curtis (default) or ‘corr’ for D_Pearson. For a given two $1 \times p$ composition vectors, $X = (X_1, \ldots, X_p)$ and $Y = (Y_1, \ldots, Y_p)$, the Bray-Curtis distance is $$BC = 1 - \sum_{i = 1}^{p}min(X_i, Y_i).$$ The $D_{Pearson}$ distance is 
$$D_{Pearson} = \frac{\sum_{i = 1}^{p}(X_i - \bar{X})(Y_i - \bar{Y})}{\{\sum_{i = 1}^{p}(X_i - \bar{X})^2\sum_{i = 1}^{p}(Y_i - \bar{Y})^2\}^{1/2}} = 1 - \frac{\sum_{i = 1}^{p}X_iY_i-(1/p)}{[\{\sum_{i = 1}^{p}X_i^2-(1/p)\}\{\sum_{i = 1}^{p}Y_i^2 - (1/p)\}]^{1/2}}$$ since $\bar{X} = 1/p$ and $\bar{Y} = 1/p$. Note, these are not metrics, dissimilarity measures would be a more accurate term to describe them.  

**print.details** = T/F, if TRUE prints details of the taxonomy matching.

Your datset should have $n$ rows (samples) and $k + p$ columns, where $p$ = number of relative abundance columns, and $k$ metadata columns such as sample ID, subject ID. The sum of relative abundances in a given row over the $p$ columns should equal 1 (up to the 6th decimal place). 

Before the distance to the reference set can be calculated, we need to match the taxonomy in your dataset to the reference set. For example, at the genus level, **Greengenes 13.8 (GG13.8)** has $p=2929$ different genus sequences (including missing ones). An example taxonomic string is $$k__Bacteria.p__Actinobacteria.c__Acidimicrobiia.o__Acidimicrobiales.f__Acidimicrobiaceae.g__Ferrimicrobium$$ which specifies the kingdom, phylum, class, order, family and genus. Another example of a taxonomic string at the $genus$ level, one that includes missing data on order, family and genus, might be $$GG13.8 is k__Bacteria.p__Acidobacteria.c__TM1.o_.f_.g_$$ which has missing information at the order, family and genus levels but is a unique genus sequence in GG13.8. At the phylum level, GG13.8 includes p=91 kingdom/phyla. A list of the 2929 sequences at the genus level and the sequences at the family (p=1116), order (p=664), class(p=319) and kingdom/phylum (p=91) are given in order in the attached file taxonomy-GG13-8-and-HMP-reference-sets.xlsx.  

### Output from **HMPdistance**
An $n$ x $2$ matrix of distances to stool and nasal HMP references (ie. with columns $dist.to.hmp.stool$, $dist.to.hmp.nasal$).

## Other comments about the inner workings of the HMPdistance function
If you decide not to use Greengenes13.8 for your microbiome analysis pipeline, use the following conventions:

1.	If you are only interested in the relative abundances at the kingdom/phylum level and want to compute distances at that level, then enter the sequences using the spelling as shown in the kingdom/phylum column of GG_13_8.xlsx (the spelling is important as it will be used for matching the relative abundance columns in the dataset to those in the reference sets. **GG_13_8_taxonomy-used-in-reference.xlsx** shows the names used in the reference sets from phylum to genus levels). Note, the kingdom/phylum names should be automatically generated by the analysis pipeline used (such as QIIME) when using GG13.8 closed reference at the phylum level. No information on order, class, family and genus should be included, only information up to the level you are working on. If the kingdom is known but the kingdom/phylum is not among the 91 string patterns in GG13.8, then the non-recognized phylum will be set to missing and its relative abundance will be added to K.\_. 

2.	In order to use the **HMPdistance** function to compute distances at a given level, you need to:

2a.	Use the naming conventions as in GG13.8 (also shown in GG_13_8_taxonomy-used-in-reference.xlsx) for each level. 

2b. At a given level, any name that is not in GG13.8 will be assigned missing values up to the lowest level at which the sequence does match a sequence in GG13.8. For example, let K.P.O.C.F.G (kingdom.phylum.order.class.family.genus) denote a particular sequence in GG13.8 and let K.P.O.C.f.g denote a sequence to be classified at the genus level but one that is not in GG13.8. The string will match the one in GG13.8 at the kingdom, phylum, order and class levels, but not at family or genus. It will be assigned K.P.O.C.f._ first, a match will be searched for and not found (since family *f* is also not in GG13.8), it will then be assigned K.P.O.C.\_.\_ and at that point a match will be found in GG13.8. This effectively means that any *leaves$ not found in the reference set will be collapsed and the relative abundances for those summed up and saved at the next level up that is found in GG13.8. Note, not all missingness patterns are in GG13.8. As a rule, if a given name string is not found in GG13.8, the lowest non-missing level will be assigned as missing and an attempt at matching will be made, and this process will repeat until a match is found or K.\_.\_.\_.\_.\_ is reached. 
2c. Be connected to the internet, since **HMPdistance** function needs to download **taxonomy-GG13-8-and-HMP-reference-sets-L2-L6.rdata** from GitHub


## Example of **HMPdistance**
This example calculates the Bray-Curtis distance from all the samples in HMP to the HMP stool and nasal reference sets: 

```r
d.l2.hmp <- HMPdistance(tax.level ='l2.phylum',
                    d.new.filename = "HMP_L2_030317_phylum-example-input-file.csv", 
                    d.new.ix.col.not.rel.abu = 1:4, 
                    measure = 'bc',
                    print.details = F)
dim(d.l2.hmp)
```

```
## [1] 681   2
```

```r
head(d.l2.hmp)
```

```
##   dist.to.hmp.stool dist.to.hmp.nasal
## 1         0.2205655         0.7363895
## 2         0.1886405         0.7537370
## 3         0.4300298         0.6236202
## 4         0.8891262         0.5393743
## 5         0.8002390         0.2495285
## 6         0.1935302         0.7418604
```

```r
summary(d.l2.hmp)
```

```
##  dist.to.hmp.stool dist.to.hmp.nasal
##  Min.   :0.1830    Min.   :0.2293   
##  1st Qu.:0.3430    1st Qu.:0.3906   
##  Median :0.5885    Median :0.6412   
##  Mean   :0.5519    Mean   :0.5841   
##  3rd Qu.:0.7179    3rd Qu.:0.7169   
##  Max.   :0.9992    Max.   :0.9984
```

Now, using an external dataset, with 1-Pearson correlation as the dissimilarity measure:

```r
d.l2.baxter <- HMPdistance(tax.level ='l2.phylum',
                    d.new.filename = "Baxter_L2_rel_abundance.csv", 
                    d.new.ix.col.not.rel.abu = 1:4, 
                    measure = 'corr',
                    print.details = F)
head(d.l2.baxter)
```

```
##   dist.to.hmp.stool dist.to.hmp.nasal
## 1         0.2014029         0.6639146
## 2         0.2353937         0.6156050
## 3         0.4965863         0.6233213
## 4         0.5338425         0.5765289
## 5         0.1386081         0.6831418
## 6         0.5358267         0.4076855
```

```r
summary(d.l2.baxter)
```

```
##  dist.to.hmp.stool dist.to.hmp.nasal
##  Min.   :0.06727   Min.   :0.1750   
##  1st Qu.:0.20084   1st Qu.:0.5693   
##  Median :0.36042   Median :0.6044   
##  Mean   :0.35411   Mean   :0.6066   
##  3rd Qu.:0.51661   3rd Qu.:0.6532   
##  Max.   :0.67852   Max.   :0.8229
```



# Testing between two (possibly correlated) groups
The *GHotelling* function calculates the generalized $Hotelling\; T^2$ statistic and resampling-based p-value for differences between two groups, A and B. A unique feature is that an individual can contribute data to both groups, as for example if an individual provides samples from two body sites, say skin sample (group A) and a saliva sample (group B). Correlations from such samples need to be taken into account, which is done by a bootstrap procedure to compute p-values.

## **GHotelling** function definition

```r
GHotelling(data, nBoot = 10000, print.details = F, seed = 1)
```
### Input to **GHotelling**
*data* = data frame or a matrix with $n$ rows that correspond to individuals who contribute data.  For a given row, columns 1 and 2 of *data* are $d_{A,stool}$ and $d_{A,nasal}$ respectively, which are the mean distances from the group A sample to the Human Microbiome Project (HMP) reference sets of 92 stool samples and 74 nasal samples, respectively. Columns 3 and 4 are $d_{B,stool}$ and $d_{B,nasal}$ respectively, which are the mean distances from the group B sample to the HMP stool and nasal reference sets. Missing data are indicated by NA. 

Only the following missingness patterns are permitted and any of these groups can have 0 or more individuals.

1. $d_{A,stool}$, $d_{A,nasal}$,  NA,  NA

2. NA, NA, $d_{B,stool}$, $d_{B,nasal}$

3. $d_{A,stool}$, $d_{A,nasal}$, $d_{B,stool}$, $d_{B,nasal}$


*nBoot* = number of bootstrap resamples to calculate the p-value. Defaul = 10000, if NULL, no p-value will be provided. Resampling is done withing groups 1, 2, and 3.

*print.details* = T/F, if T outputs intermediate matrices used in calculating the $T^2$ statistic: A, D and $\Sigma$. See the paper for details.
 
### Output from **GHotelling**
A list with two elements: $t2$ and $pval$, where $t2$ is the $T^2$ test statistic and $pval$ is the p-value based on nBoot bootstrap resamples. 

## Example of **GHotelling**
### Load test dataset (site A = nasal, site B = saliva)
Here, rows = individuals, columns = distances to HMP stool and nasal references, where columns 1 and 2 are output from HMPdistance to ref.stool and ref.nasal for 72 nasal samples, and columns 3 and 4 = output from HMPdistance to ref.stool and ref.nasal for 83 saliva samples.

In this test set we have the following groups:

- Group 1 (rows 1-6) are subjects contributing distances from nasal samples only to HMP stool and HMP nasal $(A_{stool}, A_{nasal}, NA, NA)$

- Group 3 (rows 7-72) are subjects contributing distances from nasal and saliva samples to HMP stool and HMP nasal $(A_{stool}, A_{nasal}, B_{stool}, B_{nasal})$

- Group 2 (rows 73-89) are subjects contributing distances from saliva samples only to HMP stool and HMP nasal reference sets $(NA, NA, B_{stool}, B_{nasal})$

To load the test file: 


```r
source_data("https://github.com/NCI-biostats/microbiome-fixed-reference/blob/master/test-matrix.rdata?raw=True")
```

```
## [1] "test.mx"
```

```r
dim(test.mx)
```

```
## [1] 89  4
```
And here are 3 rows from groups 1, 3, and 2:

```r
test.mx[c(1:3, 20:22, 87:89),]
```

```
##            [,1]      [,2]      [,3]      [,4]
##  [1,] 0.7434133 0.2430232        NA        NA
##  [2,] 0.7540581 0.2651087        NA        NA
##  [3,] 0.8127149 0.2575063        NA        NA
##  [4,] 0.8621425 0.2987728 0.6940892 0.7799654
##  [5,] 0.7205882 0.2913344 0.4186379 0.6471613
##  [6,] 0.8590883 0.3073629 0.4786494 0.6340253
##  [7,]        NA        NA 0.4169170 0.6724120
##  [8,]        NA        NA 0.3365672 0.6291558
##  [9,]        NA        NA 0.6546568 0.6203346
```
Here, we have subjects who contribute information to groups 1, 2 and 3:

```r
GHotelling(test.mx, nBoot = 10000, print.details = F, seed = 1)
```

```
## $t2
## [1] 702.3427
## 
## $pval
## [1] 0
```

We also handle scenarios where individuals from one of the groups are missing, for example, where we only have data for individuals in group 3 and 2, 1 and 2, or 1 and 3. Here we show an example where we do not have anyone from group 3, i.e. each individual contributes information from one body site only:

```r
GHotelling(test.mx[-(7:72),], nBoot = 10000, print.details = F, seed = 1)
```

```
## $t2
## [1] 76.25793
## 
## $pval
## [1] 0.077
```



