Chapter One masters thesis
================
Mariam Marand
2024-01-19

## Chapter One Thesis

This document outlines the metagenomic studies utilized in chapter one
of my thesis to analyze microfungi population in biological soil crust
samples at the SRER

``` r
library(tidyverse)
```

    ## Warning: package 'tidyverse' was built under R version 4.1.3

    ## Warning: package 'tibble' was built under R version 4.1.3

    ## Warning: package 'tidyr' was built under R version 4.1.3

    ## Warning: package 'readr' was built under R version 4.1.3

    ## Warning: package 'purrr' was built under R version 4.1.3

    ## Warning: package 'dplyr' was built under R version 4.1.3

    ## Warning: package 'stringr' was built under R version 4.1.3

    ## Warning: package 'forcats' was built under R version 4.1.3

    ## Warning: package 'lubridate' was built under R version 4.1.3

    ## -- Attaching core tidyverse packages ------------------------ tidyverse 2.0.0 --
    ## v dplyr     1.1.2     v readr     2.1.4
    ## v forcats   1.0.0     v stringr   1.5.0
    ## v ggplot2   3.4.3     v tibble    3.2.1
    ## v lubridate 1.9.2     v tidyr     1.3.0
    ## v purrr     1.0.1     
    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()
    ## i Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(vegan)
```

    ## Warning: package 'vegan' was built under R version 4.1.3

    ## Loading required package: permute

    ## Warning: package 'permute' was built under R version 4.1.3

    ## Loading required package: lattice
    ## This is vegan 2.6-4

``` r
library(indicspecies)
library(devtools)
```

    ## Warning: package 'devtools' was built under R version 4.1.3

    ## Loading required package: usethis
    ## 
    ## Attaching package: 'devtools'
    ## 
    ## The following object is masked from 'package:permute':
    ## 
    ##     check

``` r
library(lme4)
```

    ## Warning: package 'lme4' was built under R version 4.1.3

    ## Loading required package: Matrix

    ## Warning: package 'Matrix' was built under R version 4.1.3

    ## 
    ## Attaching package: 'Matrix'
    ## 
    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

``` r
library(patchwork)
library(lmerTest)
```

    ## Warning: package 'lmerTest' was built under R version 4.1.3

    ## 
    ## Attaching package: 'lmerTest'
    ## 
    ## The following object is masked from 'package:lme4':
    ## 
    ##     lmer
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     step

``` r
library(tibble)
##Import OTU counts table and Metadata
FullSRERDataset<-read.csv("Chapter2/FinalSRERdataOTUsandMetaData-ForAnalysis-27June2023.csv", header=TRUE, fileEncoding="UTF-8-BOM")
##Remove samples where richness=0 
filtered_FullData <- subset(FullSRERDataset, Richness != 0)
##Seperate data to generate an OTU table for distance matrix
##MetaDataTable
MetaSRER<-filtered_FullData[,1:12]
MetaSRER <- filtered_FullData[, 1:12]
rownames(MetaSRER) <- MetaSRER$Sample
#MetaSRER <- MetaSRER[, -which(names(MetaSRER) == "Sample")]

#MetaSRER <- as.data.frame(MetaSRER)
#MetaSRER <- MetaSRER %>%
  #column_to_rownames(var = "Sample")
#MetaSRER <- MetaSRER %>%
  #column_to_rownames(var = "Sample")
  #column_to_rownames(MetaSRER, var = "Sample")
 # rownames(MetaSRER) <- MetaSRER$Sample
###Convert distance to a factor since it is being used as a metric to measure disturbance levels
MetaSRER$Distance <- as.factor(MetaSRER$Distance)
###OTU Table
OTU2<-filtered_FullData[,c(1,13:350)]
rownames(OTU2) <- OTU2$Sample
###Filtering samples with low read counts  
sumtable<-function(x){
  if(is.numeric(x)){
    sum(x) > 10
  } else {
    TRUE
  }
}
CH2OTUfilter<-OTU2[, sapply(OTU2,  sumtable)]
##nothing filtered

OTU2$Sample <- NULL
shannon<-diversity(OTU2)
  #column_to_rownames(var="Sample")
###Generate Distance matrix and NMDSplot
dm.method <- 'bray'
dm <- vegdist(OTU2, method=dm.method)
otu.nmds <- metaMDS(dm,
                    k = 2,
                    maxit = 999,
                    trymax = 500,
                    wascores = TRUE)
```

    ## Run 0 stress 0.2311118 
    ## Run 1 stress 0.2303735 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01741825  max resid 0.1535975 
    ## Run 2 stress 0.233275 
    ## Run 3 stress 0.2351785 
    ## Run 4 stress 0.2412195 
    ## Run 5 stress 0.2415532 
    ## Run 6 stress 0.2376479 
    ## Run 7 stress 0.2417427 
    ## Run 8 stress 0.23099 
    ## Run 9 stress 0.2355012 
    ## Run 10 stress 0.2352679 
    ## Run 11 stress 0.2435521 
    ## Run 12 stress 0.2395707 
    ## Run 13 stress 0.2348379 
    ## Run 14 stress 0.2419242 
    ## Run 15 stress 0.2391666 
    ## Run 16 stress 0.2381306 
    ## Run 17 stress 0.232376 
    ## Run 18 stress 0.236268 
    ## Run 19 stress 0.240198 
    ## Run 20 stress 0.2304999 
    ## ... Procrustes: rmse 0.02579717  max resid 0.2468429 
    ## Run 21 stress 0.2430955 
    ## Run 22 stress 0.2420539 
    ## Run 23 stress 0.2404577 
    ## Run 24 stress 0.2354269 
    ## Run 25 stress 0.2339497 
    ## Run 26 stress 0.24448 
    ## Run 27 stress 0.2306885 
    ## ... Procrustes: rmse 0.03573483  max resid 0.2421323 
    ## Run 28 stress 0.2353041 
    ## Run 29 stress 0.2358847 
    ## Run 30 stress 0.2305817 
    ## ... Procrustes: rmse 0.03268409  max resid 0.1865969 
    ## Run 31 stress 0.2336015 
    ## Run 32 stress 0.2397618 
    ## Run 33 stress 0.2400882 
    ## Run 34 stress 0.2571544 
    ## Run 35 stress 0.2412461 
    ## Run 36 stress 0.2423363 
    ## Run 37 stress 0.2452602 
    ## Run 38 stress 0.2456492 
    ## Run 39 stress 0.2408462 
    ## Run 40 stress 0.2337724 
    ## Run 41 stress 0.2375687 
    ## Run 42 stress 0.2336519 
    ## Run 43 stress 0.2309467 
    ## Run 44 stress 0.233449 
    ## Run 45 stress 0.2339771 
    ## Run 46 stress 0.236934 
    ## Run 47 stress 0.2386717 
    ## Run 48 stress 0.2372748 
    ## Run 49 stress 0.2304597 
    ## ... Procrustes: rmse 0.006335218  max resid 0.05450306 
    ## Run 50 stress 0.2465516 
    ## Run 51 stress 0.2403738 
    ## Run 52 stress 0.2377163 
    ## Run 53 stress 0.2407909 
    ## Run 54 stress 0.235338 
    ## Run 55 stress 0.2307792 
    ## ... Procrustes: rmse 0.01696908  max resid 0.1587317 
    ## Run 56 stress 0.2414623 
    ## Run 57 stress 0.2318153 
    ## Run 58 stress 0.2340288 
    ## Run 59 stress 0.2361127 
    ## Run 60 stress 0.2356847 
    ## Run 61 stress 0.2443095 
    ## Run 62 stress 0.2302048 
    ## ... New best solution
    ## ... Procrustes: rmse 0.02196567  max resid 0.2194427 
    ## Run 63 stress 0.2449305 
    ## Run 64 stress 0.2406177 
    ## Run 65 stress 0.230817 
    ## Run 66 stress 0.2392328 
    ## Run 67 stress 0.2353067 
    ## Run 68 stress 0.2374719 
    ## Run 69 stress 0.2412288 
    ## Run 70 stress 0.2400235 
    ## Run 71 stress 0.2364655 
    ## Run 72 stress 0.2375241 
    ## Run 73 stress 0.2357548 
    ## Run 74 stress 0.2515349 
    ## Run 75 stress 0.2354987 
    ## Run 76 stress 0.2455084 
    ## Run 77 stress 0.2453695 
    ## Run 78 stress 0.2448846 
    ## Run 79 stress 0.2375694 
    ## Run 80 stress 0.244391 
    ## Run 81 stress 0.2428907 
    ## Run 82 stress 0.2407469 
    ## Run 83 stress 0.2330731 
    ## Run 84 stress 0.2370084 
    ## Run 85 stress 0.2439545 
    ## Run 86 stress 0.2426 
    ## Run 87 stress 0.2422084 
    ## Run 88 stress 0.232063 
    ## Run 89 stress 0.2375697 
    ## Run 90 stress 0.2447586 
    ## Run 91 stress 0.2309072 
    ## Run 92 stress 0.2439477 
    ## Run 93 stress 0.2399655 
    ## Run 94 stress 0.2445021 
    ## Run 95 stress 0.2420977 
    ## Run 96 stress 0.2413344 
    ## Run 97 stress 0.236744 
    ## Run 98 stress 0.2472835 
    ## Run 99 stress 0.2349245 
    ## Run 100 stress 0.2411864 
    ## Run 101 stress 0.2422868 
    ## Run 102 stress 0.231383 
    ## Run 103 stress 0.2307507 
    ## Run 104 stress 0.2411728 
    ## Run 105 stress 0.2370886 
    ## Run 106 stress 0.237503 
    ## Run 107 stress 0.236158 
    ## Run 108 stress 0.2426006 
    ## Run 109 stress 0.2368048 
    ## Run 110 stress 0.2398384 
    ## Run 111 stress 0.236792 
    ## Run 112 stress 0.2326089 
    ## Run 113 stress 0.2338237 
    ## Run 114 stress 0.2432639 
    ## Run 115 stress 0.2309159 
    ## Run 116 stress 0.2462955 
    ## Run 117 stress 0.2339289 
    ## Run 118 stress 0.2360037 
    ## Run 119 stress 0.2388 
    ## Run 120 stress 0.2302893 
    ## ... Procrustes: rmse 0.01197957  max resid 0.09225557 
    ## Run 121 stress 0.2435346 
    ## Run 122 stress 0.2352067 
    ## Run 123 stress 0.2393994 
    ## Run 124 stress 0.2373174 
    ## Run 125 stress 0.2393507 
    ## Run 126 stress 0.2405885 
    ## Run 127 stress 0.2430396 
    ## Run 128 stress 0.2341067 
    ## Run 129 stress 0.2307847 
    ## Run 130 stress 0.2304376 
    ## ... Procrustes: rmse 0.009567707  max resid 0.05130183 
    ## Run 131 stress 0.2312523 
    ## Run 132 stress 0.2309003 
    ## Run 133 stress 0.2447169 
    ## Run 134 stress 0.4144406 
    ## Run 135 stress 0.2302554 
    ## ... Procrustes: rmse 0.028714  max resid 0.2143388 
    ## Run 136 stress 0.2315082 
    ## Run 137 stress 0.2431286 
    ## Run 138 stress 0.242564 
    ## Run 139 stress 0.2339299 
    ## Run 140 stress 0.2314435 
    ## Run 141 stress 0.2361479 
    ## Run 142 stress 0.2403937 
    ## Run 143 stress 0.236026 
    ## Run 144 stress 0.2447021 
    ## Run 145 stress 0.2345461 
    ## Run 146 stress 0.2426337 
    ## Run 147 stress 0.2462669 
    ## Run 148 stress 0.2441521 
    ## Run 149 stress 0.2347469 
    ## Run 150 stress 0.4137693 
    ## Run 151 stress 0.2499002 
    ## Run 152 stress 0.2397533 
    ## Run 153 stress 0.2426637 
    ## Run 154 stress 0.2441355 
    ## Run 155 stress 0.2378219 
    ## Run 156 stress 0.2420207 
    ## Run 157 stress 0.234478 
    ## Run 158 stress 0.2350123 
    ## Run 159 stress 0.2310774 
    ## Run 160 stress 0.2330286 
    ## Run 161 stress 0.2307071 
    ## Run 162 stress 0.2333247 
    ## Run 163 stress 0.229874 
    ## ... New best solution
    ## ... Procrustes: rmse 0.00836633  max resid 0.05221768 
    ## Run 164 stress 0.2424376 
    ## Run 165 stress 0.2298225 
    ## ... New best solution
    ## ... Procrustes: rmse 0.005704613  max resid 0.05322272 
    ## Run 166 stress 0.2373191 
    ## Run 167 stress 0.235781 
    ## Run 168 stress 0.2313153 
    ## Run 169 stress 0.2421999 
    ## Run 170 stress 0.2344454 
    ## Run 171 stress 0.2412531 
    ## Run 172 stress 0.2329979 
    ## Run 173 stress 0.2311483 
    ## Run 174 stress 0.2323867 
    ## Run 175 stress 0.2313224 
    ## Run 176 stress 0.2410216 
    ## Run 177 stress 0.2441827 
    ## Run 178 stress 0.2429244 
    ## Run 179 stress 0.2343008 
    ## Run 180 stress 0.2365482 
    ## Run 181 stress 0.2474452 
    ## Run 182 stress 0.414438 
    ## Run 183 stress 0.2374456 
    ## Run 184 stress 0.2303739 
    ## Run 185 stress 0.2414447 
    ## Run 186 stress 0.2491588 
    ## Run 187 stress 0.2304066 
    ## Run 188 stress 0.2380872 
    ## Run 189 stress 0.2342438 
    ## Run 190 stress 0.2398165 
    ## Run 191 stress 0.2432829 
    ## Run 192 stress 0.2431003 
    ## Run 193 stress 0.2469993 
    ## Run 194 stress 0.2345741 
    ## Run 195 stress 0.2369605 
    ## Run 196 stress 0.2371722 
    ## Run 197 stress 0.2349055 
    ## Run 198 stress 0.2382493 
    ## Run 199 stress 0.2304925 
    ## Run 200 stress 0.2513144 
    ## Run 201 stress 0.2336272 
    ## Run 202 stress 0.2490437 
    ## Run 203 stress 0.2309564 
    ## Run 204 stress 0.248808 
    ## Run 205 stress 0.2449955 
    ## Run 206 stress 0.2381309 
    ## Run 207 stress 0.237104 
    ## Run 208 stress 0.2427599 
    ## Run 209 stress 0.2443565 
    ## Run 210 stress 0.2380698 
    ## Run 211 stress 0.2331888 
    ## Run 212 stress 0.2312506 
    ## Run 213 stress 0.2431838 
    ## Run 214 stress 0.246595 
    ## Run 215 stress 0.2355752 
    ## Run 216 stress 0.2467255 
    ## Run 217 stress 0.2427038 
    ## Run 218 stress 0.2377013 
    ## Run 219 stress 0.2411676 
    ## Run 220 stress 0.2393654 
    ## Run 221 stress 0.2302665 
    ## ... Procrustes: rmse 0.02535246  max resid 0.1930035 
    ## Run 222 stress 0.2356377 
    ## Run 223 stress 0.2406428 
    ## Run 224 stress 0.2398686 
    ## Run 225 stress 0.2424047 
    ## Run 226 stress 0.2411698 
    ## Run 227 stress 0.2383751 
    ## Run 228 stress 0.2367336 
    ## Run 229 stress 0.2417608 
    ## Run 230 stress 0.2309356 
    ## Run 231 stress 0.2320572 
    ## Run 232 stress 0.2357767 
    ## Run 233 stress 0.2303779 
    ## Run 234 stress 0.2437757 
    ## Run 235 stress 0.2312278 
    ## Run 236 stress 0.233693 
    ## Run 237 stress 0.2443046 
    ## Run 238 stress 0.2304097 
    ## Run 239 stress 0.2388379 
    ## Run 240 stress 0.2416351 
    ## Run 241 stress 0.2378239 
    ## Run 242 stress 0.2430284 
    ## Run 243 stress 0.2338863 
    ## Run 244 stress 0.2380556 
    ## Run 245 stress 0.2356052 
    ## Run 246 stress 0.2304666 
    ## Run 247 stress 0.2386415 
    ## Run 248 stress 0.2303675 
    ## Run 249 stress 0.2423746 
    ## Run 250 stress 0.2354375 
    ## Run 251 stress 0.2426519 
    ## Run 252 stress 0.2485831 
    ## Run 253 stress 0.2382367 
    ## Run 254 stress 0.2384344 
    ## Run 255 stress 0.2412545 
    ## Run 256 stress 0.2360046 
    ## Run 257 stress 0.2460906 
    ## Run 258 stress 0.2308494 
    ## Run 259 stress 0.2315841 
    ## Run 260 stress 0.2446978 
    ## Run 261 stress 0.2375311 
    ## Run 262 stress 0.2312741 
    ## Run 263 stress 0.2456552 
    ## Run 264 stress 0.2303886 
    ## Run 265 stress 0.2406095 
    ## Run 266 stress 0.239732 
    ## Run 267 stress 0.2437258 
    ## Run 268 stress 0.2545817 
    ## Run 269 stress 0.2448984 
    ## Run 270 stress 0.2360743 
    ## Run 271 stress 0.2314184 
    ## Run 272 stress 0.2312752 
    ## Run 273 stress 0.2391451 
    ## Run 274 stress 0.2412956 
    ## Run 275 stress 0.2432769 
    ## Run 276 stress 0.2592456 
    ## Run 277 stress 0.2336602 
    ## Run 278 stress 0.2344857 
    ## Run 279 stress 0.2309209 
    ## Run 280 stress 0.2424674 
    ## Run 281 stress 0.2408062 
    ## Run 282 stress 0.2402001 
    ## Run 283 stress 0.2353977 
    ## Run 284 stress 0.2468045 
    ## Run 285 stress 0.2449821 
    ## Run 286 stress 0.2305814 
    ## Run 287 stress 0.2413396 
    ## Run 288 stress 0.2324939 
    ## Run 289 stress 0.2303625 
    ## Run 290 stress 0.2411159 
    ## Run 291 stress 0.2304896 
    ## Run 292 stress 0.2380648 
    ## Run 293 stress 0.2341904 
    ## Run 294 stress 0.2301731 
    ## ... Procrustes: rmse 0.0153495  max resid 0.1597493 
    ## Run 295 stress 0.2392405 
    ## Run 296 stress 0.230496 
    ## Run 297 stress 0.2373258 
    ## Run 298 stress 0.2305121 
    ## Run 299 stress 0.2368994 
    ## Run 300 stress 0.2453037 
    ## Run 301 stress 0.2392915 
    ## Run 302 stress 0.2425217 
    ## Run 303 stress 0.2311655 
    ## Run 304 stress 0.2300983 
    ## ... Procrustes: rmse 0.006172385  max resid 0.04616822 
    ## Run 305 stress 0.245697 
    ## Run 306 stress 0.2377123 
    ## Run 307 stress 0.2413099 
    ## Run 308 stress 0.4143406 
    ## Run 309 stress 0.2392571 
    ## Run 310 stress 0.2426047 
    ## Run 311 stress 0.235137 
    ## Run 312 stress 0.239088 
    ## Run 313 stress 0.2406957 
    ## Run 314 stress 0.2341704 
    ## Run 315 stress 0.243061 
    ## Run 316 stress 0.2417634 
    ## Run 317 stress 0.2433268 
    ## Run 318 stress 0.2408812 
    ## Run 319 stress 0.230501 
    ## Run 320 stress 0.243184 
    ## Run 321 stress 0.2404494 
    ## Run 322 stress 0.2346455 
    ## Run 323 stress 0.2432994 
    ## Run 324 stress 0.2405594 
    ## Run 325 stress 0.2557522 
    ## Run 326 stress 0.2437991 
    ## Run 327 stress 0.2301508 
    ## ... Procrustes: rmse 0.01050743  max resid 0.08021531 
    ## Run 328 stress 0.2425075 
    ## Run 329 stress 0.2418526 
    ## Run 330 stress 0.2343 
    ## Run 331 stress 0.2402399 
    ## Run 332 stress 0.2489805 
    ## Run 333 stress 0.2311835 
    ## Run 334 stress 0.2511995 
    ## Run 335 stress 0.2416491 
    ## Run 336 stress 0.245379 
    ## Run 337 stress 0.2383627 
    ## Run 338 stress 0.2372063 
    ## Run 339 stress 0.2315515 
    ## Run 340 stress 0.2305251 
    ## Run 341 stress 0.2310062 
    ## Run 342 stress 0.2339642 
    ## Run 343 stress 0.2341972 
    ## Run 344 stress 0.2304921 
    ## Run 345 stress 0.2397604 
    ## Run 346 stress 0.2410775 
    ## Run 347 stress 0.232438 
    ## Run 348 stress 0.2381678 
    ## Run 349 stress 0.2481982 
    ## Run 350 stress 0.2418741 
    ## Run 351 stress 0.2372385 
    ## Run 352 stress 0.2368524 
    ## Run 353 stress 0.2354098 
    ## Run 354 stress 0.2354719 
    ## Run 355 stress 0.2412031 
    ## Run 356 stress 0.2412868 
    ## Run 357 stress 0.2385273 
    ## Run 358 stress 0.2306144 
    ## Run 359 stress 0.2379184 
    ## Run 360 stress 0.232632 
    ## Run 361 stress 0.2308701 
    ## Run 362 stress 0.242285 
    ## Run 363 stress 0.2306755 
    ## Run 364 stress 0.2413381 
    ## Run 365 stress 0.234927 
    ## Run 366 stress 0.2310655 
    ## Run 367 stress 0.2401603 
    ## Run 368 stress 0.2408759 
    ## Run 369 stress 0.2398531 
    ## Run 370 stress 0.2432096 
    ## Run 371 stress 0.2449726 
    ## Run 372 stress 0.245891 
    ## Run 373 stress 0.2542336 
    ## Run 374 stress 0.2301203 
    ## ... Procrustes: rmse 0.009052592  max resid 0.05291007 
    ## Run 375 stress 0.2385174 
    ## Run 376 stress 0.2394191 
    ## Run 377 stress 0.2400223 
    ## Run 378 stress 0.2429666 
    ## Run 379 stress 0.2393658 
    ## Run 380 stress 0.2363065 
    ## Run 381 stress 0.2303683 
    ## Run 382 stress 0.2464508 
    ## Run 383 stress 0.2334409 
    ## Run 384 stress 0.240551 
    ## Run 385 stress 0.2351237 
    ## Run 386 stress 0.2301679 
    ## ... Procrustes: rmse 0.01038614  max resid 0.0800383 
    ## Run 387 stress 0.2320109 
    ## Run 388 stress 0.242174 
    ## Run 389 stress 0.2428362 
    ## Run 390 stress 0.2444211 
    ## Run 391 stress 0.2300404 
    ## ... Procrustes: rmse 0.006715455  max resid 0.04872334 
    ## Run 392 stress 0.2316173 
    ## Run 393 stress 0.2433663 
    ## Run 394 stress 0.2434798 
    ## Run 395 stress 0.236757 
    ## Run 396 stress 0.2304115 
    ## Run 397 stress 0.2413785 
    ## Run 398 stress 0.2445988 
    ## Run 399 stress 0.2350826 
    ## Run 400 stress 0.2393314 
    ## Run 401 stress 0.2368894 
    ## Run 402 stress 0.2365663 
    ## Run 403 stress 0.2304188 
    ## Run 404 stress 0.2332653 
    ## Run 405 stress 0.2346239 
    ## Run 406 stress 0.2429366 
    ## Run 407 stress 0.2367258 
    ## Run 408 stress 0.2312693 
    ## Run 409 stress 0.2392892 
    ## Run 410 stress 0.236325 
    ## Run 411 stress 0.2433748 
    ## Run 412 stress 0.2350968 
    ## Run 413 stress 0.2436169 
    ## Run 414 stress 0.2395468 
    ## Run 415 stress 0.2362774 
    ## Run 416 stress 0.2328546 
    ## Run 417 stress 0.2431526 
    ## Run 418 stress 0.2507738 
    ## Run 419 stress 0.2301508 
    ## ... Procrustes: rmse 0.0105057  max resid 0.0802165 
    ## Run 420 stress 0.2343923 
    ## Run 421 stress 0.2419666 
    ## Run 422 stress 0.2300411 
    ## ... Procrustes: rmse 0.006776056  max resid 0.05226509 
    ## Run 423 stress 0.2458865 
    ## Run 424 stress 0.2437823 
    ## Run 425 stress 0.2407967 
    ## Run 426 stress 0.2460695 
    ## Run 427 stress 0.2402652 
    ## Run 428 stress 0.4143446 
    ## Run 429 stress 0.2447304 
    ## Run 430 stress 0.2461257 
    ## Run 431 stress 0.2313182 
    ## Run 432 stress 0.2382869 
    ## Run 433 stress 0.2329521 
    ## Run 434 stress 0.2451315 
    ## Run 435 stress 0.2396891 
    ## Run 436 stress 0.246772 
    ## Run 437 stress 0.2365313 
    ## Run 438 stress 0.2379107 
    ## Run 439 stress 0.230572 
    ## Run 440 stress 0.2307174 
    ## Run 441 stress 0.2522182 
    ## Run 442 stress 0.2304801 
    ## Run 443 stress 0.2467897 
    ## Run 444 stress 0.239536 
    ## Run 445 stress 0.2407775 
    ## Run 446 stress 0.2375594 
    ## Run 447 stress 0.2392188 
    ## Run 448 stress 0.2301612 
    ## ... Procrustes: rmse 0.009887652  max resid 0.08124835 
    ## Run 449 stress 0.2368338 
    ## Run 450 stress 0.237783 
    ## Run 451 stress 0.2359125 
    ## Run 452 stress 0.2326325 
    ## Run 453 stress 0.2315729 
    ## Run 454 stress 0.2421263 
    ## Run 455 stress 0.2300732 
    ## ... Procrustes: rmse 0.007339121  max resid 0.05017588 
    ## Run 456 stress 0.2350106 
    ## Run 457 stress 0.2429585 
    ## Run 458 stress 0.2436109 
    ## Run 459 stress 0.2428232 
    ## Run 460 stress 0.2392491 
    ## Run 461 stress 0.2380074 
    ## Run 462 stress 0.2302093 
    ## ... Procrustes: rmse 0.00861403  max resid 0.0498273 
    ## Run 463 stress 0.2309412 
    ## Run 464 stress 0.2444173 
    ## Run 465 stress 0.239286 
    ## Run 466 stress 0.2384645 
    ## Run 467 stress 0.2400536 
    ## Run 468 stress 0.2429844 
    ## Run 469 stress 0.2338958 
    ## Run 470 stress 0.2408074 
    ## Run 471 stress 0.2396013 
    ## Run 472 stress 0.242461 
    ## Run 473 stress 0.2310891 
    ## Run 474 stress 0.2403923 
    ## Run 475 stress 0.2399462 
    ## Run 476 stress 0.2330202 
    ## Run 477 stress 0.2409946 
    ## Run 478 stress 0.244358 
    ## Run 479 stress 0.2351705 
    ## Run 480 stress 0.2319743 
    ## Run 481 stress 0.2361379 
    ## Run 482 stress 0.2363276 
    ## Run 483 stress 0.2387586 
    ## Run 484 stress 0.2367824 
    ## Run 485 stress 0.2361304 
    ## Run 486 stress 0.2444793 
    ## Run 487 stress 0.234287 
    ## Run 488 stress 0.2407801 
    ## Run 489 stress 0.2342997 
    ## Run 490 stress 0.2311953 
    ## Run 491 stress 0.2350437 
    ## Run 492 stress 0.2351098 
    ## Run 493 stress 0.2325239 
    ## Run 494 stress 0.2328711 
    ## Run 495 stress 0.2303603 
    ## Run 496 stress 0.2398865 
    ## Run 497 stress 0.2303344 
    ## Run 498 stress 0.2425575 
    ## Run 499 stress 0.2315089 
    ## Run 500 stress 0.2308606 
    ## *** Best solution was not repeated -- monoMDS stopping criteria:
    ##    500: stress ratio > sratmax

``` r
##Check model
stressplot(otu.nmds)
```

![](Chapter_One_Thesis_files/figure-gfm/cars-1.png)<!-- -->

``` r
###NMDS all Pastures
```

## Including Plots

You can also embed plots, for example:

    ## species scores not available

![](Chapter_One_Thesis_files/figure-gfm/plot-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
