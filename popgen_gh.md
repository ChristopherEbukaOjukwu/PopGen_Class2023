Problem 2
================
Christopher Ebuka Ojukwu

``` r
library(genio)
```

\#Question 1

``` r
#a. Calculate the allele frequencies. (1pt):
AA = 8  
Aa = 38
aa = 121
n = sum(c(AA, Aa, aa))  # == 167 
p <- (2*AA)/(2*n) + (Aa)/(2*n) # == 0.1616
q <- (2*aa)/(2*n) + (Aa)/(2*n) # == 0.8383
p
```

    ## [1] 0.1616766

``` r
q
```

    ## [1] 0.8383234

``` r
# b. Perform a chi-square test to determine if the genotype frequencies conform to Hardy-Weinberg expectations. (1pt)
FAA <- p^2
FAa <- 2*p*q
Faa <- q^2
obs.genos<-c(AA, Aa, aa)
exp.genos<- c(n*p^2, n*2*p*q, n*q^2) 
X2<- sum((obs.genos-exp.genos)^2 / exp.genos) #chisquare == 4.3063
df <- 1
px2 <- pchisq(X2,df,lower.tail=F) #== 0.0379
X2
```

    ## [1] 4.306359

``` r
px2
```

    ## [1] 0.03797014

\#Question 2

``` r
#Transferred the table 1 into a text file and reading it in.
my_data <- read.delim("~/popgengh/data.txt") 
N<-nrow(my_data) # == 40
my_data
```

    ##    Indiv..Idh1 Idh2 Mdh
    ## 1           Aa   AA  AA
    ## 2           Aa   AA  AA
    ## 3           AA   AA  AA
    ## 4           AA   AA  AA
    ## 5           Aa   AA  Aa
    ## 6           AA   AA  AA
    ## 7           AA   AA  AA
    ## 8           Aa   AA  Aa
    ## 9           AA   Aa  AA
    ## 10          AA   AA  AA
    ## 11          Aa   AA  AA
    ## 12          Aa   AA  Aa
    ## 13          Aa   AA  AA
    ## 14          AA   AA  AA
    ## 15          Aa   AA  Aa
    ## 16          Aa   AA  AA
    ## 17          Aa   AA  Aa
    ## 18          Aa   AA  AA
    ## 19          Aa   AA  AA
    ## 20          Aa   AA  AA
    ## 21          AA   Aa  AA
    ## 22          Aa   AA  Aa
    ## 23          Aa   AA  AA
    ## 24          AA   AA  AA
    ## 25          AA   AA  AA
    ## 26          AA   AA  Aa
    ## 27          AA   AA  AA
    ## 28          AA   AA  Aa
    ## 29          AA   AA  Aa
    ## 30          AA   AA  AA
    ## 31          AA   AA  AA
    ## 32          AA   AA  AA
    ## 33          AA   AA  AA
    ## 34          AA   AA  AA
    ## 35          AA   AA  Aa
    ## 36          AA   AA  AA
    ## 37          AA   AA  AA
    ## 38          AA   AA  Aa
    ## 39          AA   AA  AA
    ## 40          AA   AA  AA

``` r
N
```

    ## [1] 40

``` r
#a. What is the estimated heterozygosity for each of the loci? (2pts)
Aa_in_Idh1 <- sum(my_data[,1] == "Aa") #calculates number of heterozygosity in Idh1 == 15
estimated_Aa_in_Idh1 <- Aa_in_Idh1/N #gets the estimated heterozygosity == 0.375
Aa_in_Idh2 <- sum(my_data[,2] == "Aa") #calculates number of heterozygosity in Idh2 == 2
estimated_Aa_in_Idh2 <- Aa_in_Idh2/N #gets the estimated heterozygosity = 0.05
Aa_in_mdh <- sum(my_data[,3] == "Aa") #calculates number of heterozygosity in mdh == 11
estimated_Aa_in_mdh <- Aa_in_mdh/N #gets the estimated heterozygosity == 0.275
estimated_Aa_in_mdh 
```

    ## [1] 0.275

``` r
#2b. What is the overall estimate of heterozygosity across all loci? (2pts)
sum_overall_estimate <- estimated_Aa_in_Idh1 + estimated_Aa_in_Idh2 + estimated_Aa_in_mdh # == 0.7
average_overall_estimate <- mean(c(estimated_Aa_in_Idh1, estimated_Aa_in_Idh2, estimated_Aa_in_mdh)) # == 0.2333
average_overall_estimate
```

    ## [1] 0.2333333

\#Question 3

``` r
#Reading in data for Table2
table2_data <- read.delim("~/popgengh/data_table2.txt")
table2_data
```

    ##        Locus   FST
    ## 1        pgm 0.028
    ## 2        pgi 0.052
    ## 3         hk 0.291
    ## 4        got 0.017
    ## 5         ak 0.062
    ## 6        bdh 0.034
    ## 7 alpha-gpdh 0.027
    ## 8         to 0.035

``` r
#3a. Estimate the total number of migrants per generation for each locus separately. (2pt)
#Nm = (1 - fst)/(4*fst)

pgm_fst <- table2_data[1,2] #=0.028
Nm_pgm <- (1 - pgm_fst)/(4*pgm_fst) #== 8.6785

pgi_fst <- table2_data[2,2] #=0.052
Nm_pgi <- (1 - pgi_fst)/(4*pgi_fst)  #==4.5576

hk_fst <- table2_data[3,2] #=0.291
Nm_hk <- (1 - hk_fst)/(4*hk_fst) #== 0.609

got_fst <- table2_data[4,2] #== 0.017
Nm_got <- (1 - got_fst)/(4*got_fst) #== 14.4558

ak_fst <- table2_data[5,2] #== 0.062
Nm_ak <- (1 - ak_fst)/(4*ak_fst)  #== 3.7822

bdh_fst <- table2_data[6,2] #== 0.034
Nm_bdh <- (1 - bdh_fst)/(4*bdh_fst) #== 7.1029

alphagpdh_fst <- table2_data[7,2] #== 0.027
Nm_alphagpdh <- (1 - alphagpdh_fst)/(4*alphagpdh_fst) #== 9.0092

to_fst <- table2_data[8,2] #== 0.035
Nm_to <- (1 - to_fst)/(4*to_fst) #== 6.8928
```

``` r
#3b. Calculate an overall estimate of the number of migrants per generation (combining all the data from all loci). (1pt)

overall_Nm <- sum(Nm_pgm, Nm_pgi, Nm_ak, Nm_alphagpdh, Nm_bdh, Nm_got, Nm_hk, Nm_to)/8#==6.8860
overall_Nm
```

    ## [1] 6.886071

``` r
#Question 4
#Reading in data for Table2
table3_data <- read.delim("~/popgengh/data_table3.txt")
N <- 43

#4a. Calculate FIS, FST, and FIT for these data. 

#allele frequency in each subpopulation
#population 1 - 40
pi1 <- table3_data[1,2]
#population 41
pi2 <- table3_data[2,2]
#population 42
pi3 <- table3_data[3,2]
#population 43
pi4 <- table3_data[4,2]

#population 1 - 40
qi1 <- 1 - pi1
#population 41
qi2 <- 1 - pi2
#population 42
qi3 <- 1 - pi3
#population 43
qi4 <- 1 - pi4

#observed heterozygosity in each subpopulation
#population 1 - 40
Hi1 <- table3_data[1,3]
#population 41
Hi2 <- table3_data[2,3]
#population 42
Hi3 <- table3_data[3,3]
#population 43
Hi4 <- table3_data[4,3]

Hi <- (Hi1 + Hi2 + Hi3 + Hi4)/N

#expected heterozygosity in each population
#population 1 - 40
He1 <- 2 * pi1 * qi1
#population 41
He2 <- 2 * pi2 * qi2
#population 42
He3 <- 2 * pi3 * qi3
#population 43
He4 <- 2 * pi4 * qi4

#Get the allele frequency across whole population, and expected total heterozygosity
Hs<- (He1 + He2 + He3 + He4)/N
ptotal <- ((pi1*40) + pi2 + pi3 + pi4)/N
Ht <- 2 * ptotal * (1 - ptotal)

Fst <- (Ht - Hs)/Ht #== 0.3746
Fis <- ((Hs) - (Hi))/(Hs) #== 0.6933
Fit <- Fis + Fst-(Fis * Fst) #== 0.8082

table3_data
```

    ##   Subpopulation   pi   HI
    ## 1         #1-40 1.00 0.00
    ## 2           #41 0.49 0.17
    ## 3           #42 0.83 0.06
    ## 4           #43 0.91 0.06

``` r
Fst
```

    ## [1] 0.374646

``` r
Fis
```

    ## [1] 0.6933813

``` r
Fit
```

    ## [1] 0.8082547
