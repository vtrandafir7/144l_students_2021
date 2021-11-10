Phyloseq
================
Victor Trandafir
11/10/2021

# Intro

We explore the processed ACIDD 16S sequences using
[phyloseq](https://joey711.github.io/phyloseq/)

# Install phyloseq

``` r
#BiocManager::install("phyloseq")
```

``` r
library(tidyverse) 
library(phyloseq)
library(RColorBrewer)      
library(readxl)
library(lubridate)
```

# Import Data

``` r
count.tab <- read_rds("~/Documents/College/Fourth Year/EEMB 144L/Github/144l_students_2021/Input_Data/week7/seqtab-nochimtaxa.rds") #table of counts for each sequence in each sample
tax.tab <- read_rds("~/Documents/College/Fourth Year/EEMB 144L/Github/144l_students_2021/Input_Data/week7/taxa.rds") #table that matches ASV to sequence

#you will need to download the ACIDD_Exp_Processed_DOC_BGE.rds and the ACIDD_Exp_BactAbund.xlsx files from Justine's Github repository (located in Input_Data/week7/) and put them into your own local Input Data folder
metadata_OG <- read_excel("~/Documents/College/Fourth Year/EEMB 144L/Github/144l_students_2021/Input_Data/week7/ACIDD_Exp_BactAbund.xlsx", sheet = "Metadata")
glimpse(metadata_OG)
```

    ## Rows: 84
    ## Columns: 18
    ## $ Experiment              <chr> "ASH171", "ASH171", "ASH171", "ASH171", "ASH17…
    ## $ Location                <chr> "San Diego", "San Diego", "San Diego", "San Di…
    ## $ Temperature_C           <dbl> 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15…
    ## $ Depth                   <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5…
    ## $ Bottle                  <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", "…
    ## $ Timepoint               <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, …
    ## $ Treatment               <chr> "Control", "Control", "Control", "Control", "C…
    ## $ Target_DOC_Amendment_uM <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ Inoculum_L              <dbl> 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1…
    ## $ Media_L                 <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4…
    ## $ Datetime                <chr> "2017-12-16T21:30", "2017-12-17T10:00", "2017-…
    ## $ TOC_Sample              <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE…
    ## $ DOC_Sample              <lgl> TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, …
    ## $ Parallel_Sample         <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE…
    ## $ Cell_Sample             <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALS…
    ## $ DNA_Sample              <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,…
    ## $ Nutrient_Sample         <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,…
    ## $ DNA_SampleID            <chr> "ASH171-A0_S293", NA, NA, NA, NA, NA, "ASH171-…

``` r
metadata <- metadata_OG %>%
  mutate(Datetime = ymd_hm(Datetime))
glimpse(metadata)
```

    ## Rows: 84
    ## Columns: 18
    ## $ Experiment              <chr> "ASH171", "ASH171", "ASH171", "ASH171", "ASH17…
    ## $ Location                <chr> "San Diego", "San Diego", "San Diego", "San Di…
    ## $ Temperature_C           <dbl> 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15…
    ## $ Depth                   <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5…
    ## $ Bottle                  <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", "…
    ## $ Timepoint               <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, …
    ## $ Treatment               <chr> "Control", "Control", "Control", "Control", "C…
    ## $ Target_DOC_Amendment_uM <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ Inoculum_L              <dbl> 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1…
    ## $ Media_L                 <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4…
    ## $ Datetime                <dttm> 2017-12-16 21:30:00, 2017-12-17 10:00:00, 201…
    ## $ TOC_Sample              <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE…
    ## $ DOC_Sample              <lgl> TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, …
    ## $ Parallel_Sample         <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE…
    ## $ Cell_Sample             <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALS…
    ## $ DNA_Sample              <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,…
    ## $ Nutrient_Sample         <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,…
    ## $ DNA_SampleID            <chr> "ASH171-A0_S293", NA, NA, NA, NA, NA, "ASH171-…

``` r
metadata_doc <- read_rds("~/Documents/College/Fourth Year/EEMB 144L/Github/144l_students_2021/Input_Data/week7/ACIDD_Exp_Processed_DOC_BGE.rds")
glimpse(metadata_doc)
```

    ## Rows: 84
    ## Columns: 64
    ## $ Experiment              <chr> "ASH171", "ASH171", "ASH171", "ASH171", "ASH17…
    ## $ Location                <chr> "San Diego", "San Diego", "San Diego", "San Di…
    ## $ Temperature_C           <dbl> 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15…
    ## $ Depth                   <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5…
    ## $ Bottle                  <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", "…
    ## $ Timepoint               <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, …
    ## $ Treatment               <chr> "Control", "Control", "Control", "Control", "C…
    ## $ Target_DOC_Amendment_uM <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ Inoculum_L              <dbl> 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1…
    ## $ Media_L                 <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4…
    ## $ Datetime                <dttm> 2017-12-16 21:30:00, 2017-12-17 10:00:00, 201…
    ## $ hours                   <dbl> 0.0, 12.5, 25.5, 48.0, 71.5, 95.0, 118.5, 122.…
    ## $ days                    <dbl> 0.0000000, 0.5208333, 1.0625000, 2.0000000, 2.…
    ## $ TOC                     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ sd_TOC                  <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ PTOC                    <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ sd_PTOC                 <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ cells                   <dbl> 1.30e+08, 1.34e+08, 1.28e+08, 1.55e+08, 1.55e+…
    ## $ sd_cells                <dbl> 20900000, 27600000, 22200000, 25200000, 319000…
    ## $ ln_cells                <dbl> 18.68305, 18.71335, 18.66754, 18.85894, 18.858…
    ## $ diff_ln_cells           <dbl> 0.000000000, 0.030305349, -0.045809536, 0.1913…
    ## $ bc                      <dbl> 0.32500, 0.33500, 0.32000, 0.38750, 0.38750, 0…
    ## $ ave_bc                  <dbl> 0.308750, 0.330000, 0.317500, 0.358750, 0.4100…
    ## $ sd_bc                   <dbl> 0.022980970, 0.007071068, 0.003535534, 0.04065…
    ## $ exp_start               <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5…
    ## $ exp_end                 <dbl> 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6…
    ## $ ln_cells_exp_start      <dbl> 18.85894, 18.85894, 18.85894, 18.85894, 18.858…
    ## $ ln_cells_exp_end        <dbl> 19.74776, 19.74776, 19.74776, 19.74776, 19.747…
    ## $ cells_exp_start         <dbl> 1.55e+08, 1.55e+08, 1.55e+08, 1.55e+08, 1.55e+…
    ## $ cells_exp_end           <dbl> 3.77e+08, 3.77e+08, 3.77e+08, 3.77e+08, 3.77e+…
    ## $ days_exp_start          <dbl> 2.979167, 2.979167, 2.979167, 2.979167, 2.9791…
    ## $ days_exp_end            <dbl> 4.937500, 4.937500, 4.937500, 4.937500, 4.9375…
    ## $ mew                     <dbl> 0.4538656, 0.4538656, 0.4538656, 0.4538656, 0.…
    ## $ doubling                <dbl> 1.5272081, 1.5272081, 1.5272081, 1.5272081, 1.…
    ## $ delta_cells             <dbl> 247000000, 247000000, 247000000, 247000000, 24…
    ## $ delta_bc                <dbl> 0.61750, 0.61750, 0.61750, 0.61750, 0.61750, 0…
    ## $ ave_mew                 <dbl> 0.5441279, 0.5441279, 0.5441279, 0.5441279, 0.…
    ## $ sd_mew                  <dbl> 0.09366960, 0.09366960, 0.09366960, 0.09366960…
    ## $ ave_doubling            <dbl> 1.3099139, 1.3099139, 1.3099139, 1.3099139, 1.…
    ## $ sd_doubling             <dbl> 0.22549685, 0.22549685, 0.22549685, 0.22549685…
    ## $ ave_delta_cells         <dbl> 232500000, 232500000, 232500000, 232500000, 23…
    ## $ sd_delta_cells          <dbl> 15047361, 15047361, 15047361, 15047361, 150473…
    ## $ ave_delta_bc            <dbl> 0.581250, 0.581250, 0.581250, 0.581250, 0.5812…
    ## $ sd_delta_bc             <dbl> 0.037618403, 0.037618403, 0.037618403, 0.03761…
    ## $ ave_lag                 <dbl> 3.46875, 3.46875, 3.46875, 3.46875, 3.46875, 3…
    ## $ sd_lag                  <dbl> 0.5080646, 0.5080646, 0.5080646, 0.5080646, 0.…
    ## $ interp_toc              <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ interp_bc               <dbl> 0.3, 0.3, 0.3, 0.4, 0.4, 0.5, 0.9, NA, NA, NA,…
    ## $ doc                     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ bioav_doc               <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ doc_exp_end             <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ delta_doc               <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ tdelta_doc              <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ bge                     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ ave_toc                 <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ sd_toc                  <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ ave_bioav_doc           <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ sd_bioav_doc            <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ ave_delta_doc           <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ sd_delta_doc            <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ ave_tdelta_doc          <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ sd_tdelta_doc           <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ ave_bge                 <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
    ## $ sd_bge                  <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…

``` r
join <- left_join(metadata_doc, metadata)
sample.tab <- join %>%
  drop_na(DNA_SampleID) %>% 
  column_to_rownames(var = "DNA_SampleID") 
```

# Phyloseq Object

We need to create a phyloseq object that merges all three datasets.
Sometimes this doesn’t work beacuse of the format of the data files.
Make sure all the sample names between the sampleinfo.txt and
seqtab-nochimtaxa.txt are the same

``` r
OTU = otu_table(count.tab, taxa_are_rows = TRUE) 
TAX = tax_table(tax.tab)
SAM = sample_data(sample.tab)
ps = phyloseq(OTU,TAX,SAM) 
```

# Filter sequences

We will filter out chloroplasts and mitochondria, because we only
intended to amplify bacterial sequences. It’s good to check you don’t
have anything lurking in the taxonomy table.

``` r
sub_ps <- ps %>%
  # subset_samples(Experiment == "ASH172") %>%  #use this function if you want to only include some subset of your sample set in the subsequent analysis
  subset_taxa(Family  != "mitochondria" & Order  != "Chloroplast")
```

Q1: What did we do in the code chunk above? (What is the difference
between phyloseq object ps and phyloseq object sub_ps?) Why did we do
this?

A1: We removed mitochondria and chloroplasts by creating a new phyloseq
object ‘sub_ps’ to focus our analysis on only bacterial sequences.

# Sample Summary

As a first analysis, we will look at the distribution of read counts
from our samples

<img src="phyloseq_vtrandafir_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />
Q2: Describe what “sequencing depth” means in your own words.

A2: The number of times a nucleotide was detected.

``` r
# mean, max and min of sample read counts
summary(sample_sum_df)
```

    ##       sum      
    ##  Min.   :6047  
    ##  1st Qu.:7153  
    ##  Median :7696  
    ##  Mean   :7696  
    ##  3rd Qu.:8453  
    ##  Max.   :9250

# Beta Diversity

Beta diversity involves calculating metrics such as distances or
dissimilarities based on pairwise comparisons of samples – they don’t
exist for a single sample, but rather only as metrics that relate
samples to each other. i.e. beta diversity = patterns in community
structure between samples

Since differences in sampling depths between samples can influence
distance/dissimilarity metrics, we first need to somehow normalize the
read depth across our samples.

## Subsample

We will rarefy (random subsample with replacement) the read depth of the
samples first (scale to the smallest library size).

[Case for not
subsampling](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531)

[Response blog for
subsampling](https://www.polarmicrobes.org/how-i-learned-to-stop-worrying-and-love-subsampling-rarifying/)

Read depth is an artefact of a machine made by a company in San Diego,
not anything about your samples or their biology. It is totally
artifactual, and controlling for artifacts is critical in science.
Subsampling randomly is the simplest way to control for this, and the
question is whether this is the “best” way of controlling for it. See
links above for alternative arguments about what the best way of
controlling for this artefact is.

A strong reason to subsample is to standardize effort. The bottom line
is that in all experimental design you should not be comparing things to
which you devote different effort in resolution. For instance, you don’t
sample one site once a week and another once a month if you want to
compare the dynamics between the sites. You standardize effort.

With that said, the bigger your differential in mean (or median) read
depth (reads/sample) between pre- and post-subsampling, the greater the
“effect” on beta diversity.

Examples:

-   means reads before = 40k, mean reads after = 1k, big effect.
-   mean reads before = 40k, mean reads after = 20k, small effect.
-   mean reads before = 2k, mean reads after = 1k, small effect.

We will subsample to the minimum read depth of all samples and not
subsample. We’ll then compare the mean reads pre- and post-subsampling
and also compare beta diversity patterns

``` r
ps_min <-  rarefy_even_depth(sub_ps, sample.size = min(sample_sums(sub_ps)))
```

    ## You set `rngseed` to FALSE. Make sure you've set & recorded
    ##  the random seed of your session for reproducibility.
    ## See `?set.seed`

    ## ...

    ## 1OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

``` r
mean(sample_sums(sub_ps)) #7686
```

    ## [1] 7695.625

``` r
mean(sample_sums(ps_min)) #6048 this is also the same as min(sample_sums(sub)ps) 
```

    ## [1] 6047

Q3: Do you think that the subsampling we did here will have a large
effect on our beta diveristy analyses? why or why not?

A3: There will be a small effect because there is a small difference
(1649) in the mean between before and after sub-sampling.

## NMDS

One of the best exploratory analyses for amplicon data is unconstrained
ordinations. Here we will look at non-metric multidimensional scaling
(NMDS) ordinations of our full community samples. For NMDS plots it’s
important to set a seed since the starting positions of samples in the
algorithm is random.

``` r
set.seed(1)
# Ordinate
nmds <- ordinate(sub_ps, method = "NMDS",  distance = "bray") # stress = 0.04
```

    ## Square root transformation
    ## Wisconsin double standardization
    ## Run 0 stress 0.04012717 
    ## Run 1 stress 0.1688387 
    ## Run 2 stress 0.03976973 
    ## ... New best solution
    ## ... Procrustes: rmse 0.02529697  max resid 0.05228924 
    ## Run 3 stress 0.03976973 
    ## ... Procrustes: rmse 1.166456e-06  max resid 2.680971e-06 
    ## ... Similar to previous best
    ## Run 4 stress 0.3118493 
    ## Run 5 stress 0.04012714 
    ## ... Procrustes: rmse 0.02534126  max resid 0.05272595 
    ## Run 6 stress 0.04012714 
    ## ... Procrustes: rmse 0.02536051  max resid 0.05275186 
    ## Run 7 stress 0.03976973 
    ## ... New best solution
    ## ... Procrustes: rmse 2.282144e-06  max resid 6.792273e-06 
    ## ... Similar to previous best
    ## Run 8 stress 0.03976973 
    ## ... New best solution
    ## ... Procrustes: rmse 2.181449e-06  max resid 4.257597e-06 
    ## ... Similar to previous best
    ## Run 9 stress 0.03976973 
    ## ... Procrustes: rmse 2.195192e-06  max resid 6.101203e-06 
    ## ... Similar to previous best
    ## Run 10 stress 0.04012714 
    ## ... Procrustes: rmse 0.02534183  max resid 0.05272568 
    ## Run 11 stress 0.04012713 
    ## ... Procrustes: rmse 0.02535484  max resid 0.05274365 
    ## Run 12 stress 0.1745143 
    ## Run 13 stress 0.04012713 
    ## ... Procrustes: rmse 0.02535307  max resid 0.05274128 
    ## Run 14 stress 0.04012714 
    ## ... Procrustes: rmse 0.02536418  max resid 0.05275666 
    ## Run 15 stress 0.04012714 
    ## ... Procrustes: rmse 0.02532494  max resid 0.05270215 
    ## Run 16 stress 0.182315 
    ## Run 17 stress 0.04012713 
    ## ... Procrustes: rmse 0.02535024  max resid 0.05273699 
    ## Run 18 stress 0.04012714 
    ## ... Procrustes: rmse 0.02532929  max resid 0.05270813 
    ## Run 19 stress 0.04012716 
    ## ... Procrustes: rmse 0.02530767  max resid 0.05267768 
    ## Run 20 stress 0.04012714 
    ## ... Procrustes: rmse 0.02533424  max resid 0.05271508 
    ## *** Solution reached

``` r
set.seed(1)
# Ordinate
nmds_min <- ordinate(ps_min, method = "NMDS",  distance = "bray") # stress = 0.04
```

    ## Square root transformation
    ## Wisconsin double standardization
    ## Run 0 stress 0.04245971 
    ## Run 1 stress 0.1692239 
    ## Run 2 stress 0.04245971 
    ## ... New best solution
    ## ... Procrustes: rmse 7.479737e-06  max resid 1.654479e-05 
    ## ... Similar to previous best
    ## Run 3 stress 0.04245971 
    ## ... Procrustes: rmse 7.202057e-07  max resid 2.03642e-06 
    ## ... Similar to previous best
    ## Run 4 stress 0.04245971 
    ## ... Procrustes: rmse 1.49017e-05  max resid 3.362063e-05 
    ## ... Similar to previous best
    ## Run 5 stress 0.04245971 
    ## ... New best solution
    ## ... Procrustes: rmse 5.658427e-06  max resid 1.243771e-05 
    ## ... Similar to previous best
    ## Run 6 stress 0.04245971 
    ## ... New best solution
    ## ... Procrustes: rmse 2.850551e-06  max resid 5.73355e-06 
    ## ... Similar to previous best
    ## Run 7 stress 0.04245971 
    ## ... Procrustes: rmse 1.115057e-06  max resid 2.274338e-06 
    ## ... Similar to previous best
    ## Run 8 stress 0.04245971 
    ## ... Procrustes: rmse 6.155251e-06  max resid 1.42328e-05 
    ## ... Similar to previous best
    ## Run 9 stress 0.04245971 
    ## ... Procrustes: rmse 1.098891e-05  max resid 2.662573e-05 
    ## ... Similar to previous best
    ## Run 10 stress 0.04245971 
    ## ... Procrustes: rmse 3.169751e-06  max resid 7.147339e-06 
    ## ... Similar to previous best
    ## Run 11 stress 0.04245971 
    ## ... Procrustes: rmse 8.642031e-07  max resid 1.755594e-06 
    ## ... Similar to previous best
    ## Run 12 stress 0.1742776 
    ## Run 13 stress 0.04245971 
    ## ... Procrustes: rmse 5.187519e-06  max resid 1.18356e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.04245971 
    ## ... Procrustes: rmse 2.609999e-06  max resid 4.589371e-06 
    ## ... Similar to previous best
    ## Run 15 stress 0.04245971 
    ## ... Procrustes: rmse 6.83622e-06  max resid 1.464468e-05 
    ## ... Similar to previous best
    ## Run 16 stress 0.04245971 
    ## ... Procrustes: rmse 2.799156e-06  max resid 5.994519e-06 
    ## ... Similar to previous best
    ## Run 17 stress 0.04245971 
    ## ... Procrustes: rmse 4.212666e-06  max resid 9.546163e-06 
    ## ... Similar to previous best
    ## Run 18 stress 0.04245971 
    ## ... Procrustes: rmse 4.573636e-06  max resid 9.663307e-06 
    ## ... Similar to previous best
    ## Run 19 stress 0.04245971 
    ## ... Procrustes: rmse 2.314194e-06  max resid 4.827317e-06 
    ## ... Similar to previous best
    ## Run 20 stress 0.04245971 
    ## ... Procrustes: rmse 1.555441e-06  max resid 2.831626e-06 
    ## ... Similar to previous best
    ## *** Solution reached

<img src="phyloseq_vtrandafir_files/figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

<img src="phyloseq_vtrandafir_files/figure-gfm/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

NMDS plots attempt to show ordinal distances between samples as
accurately as possible in two dimensions. It is important to report the
stress of these plots, because a high stress value means that the
algorithm had a hard time representing the distances between samples in
2 dimensions. The stress of this plot was good - it was .04 (generally
anything below .2 is considered acceptable).

Q4: Which datasets were used to make the two NMDS plots above? Based on
how these two plots look, which dataset are we going to move forward
with and why?

A4: The subsetted data (sub_ps) and read depth sub-sampled data (ps_min)
were used to make these plots. We will use ps_min because there is not a
significant change in the diversity estimates.

Q5: Describe in a couple of sentences what patterns you see in the beta
diversity of the communities in the control and Ash Leachate treatments.

A5: There is generally more diversity in the Ash Leachate treatment than
in the control treatment. Diversity also appears to increase with time
for both treatments. Finally, the treatments from San Diego have higher
diversity than those in Santa Barbara.

# Alpha Diversity

Estimating alpha diversity of microbial communities is
[problematic](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC93182/) no
matter what you do.

We are going to calculate the Chao1 index for richness and the Shannon
diversity index.

**it is important to note that the alpha diversity values are not
interpretable as “real” numbers of anything (due to the nature of
amplicon data), but they can still be useful as relative metrics of
comparison. If Chao1 richness goes up, but Shannon diversity goes down,
it indicates that the sample may have more ASVs but is dominated by a
few of them.**

We will use the subsampled library, which retains estimates of the
species abundance of the real population while standardizing sampling
effort.

[subsampling and alpha diversity
paper](https://www.frontiersin.org/articles/10.3389/fmicb.2019.02407/full)

[Chao1: nonparametric estimation of minimum community
richness](https://www.jstor.org/stable/4615964?seq=1#metadata_info_tab_contents)

``` r
richness <- estimate_richness(ps_min, measures = c("Chao1", "Shannon")) %>% 
  rownames_to_column(., var = "DNA_ID") %>% 
  mutate_at(vars(DNA_ID), str_replace_all, pattern = "171.", "171-") %>% 
   mutate_at(vars(DNA_ID), str_replace_all, pattern = "172.", "172-")
```

Let’s add the sample metadata into this dataframe

``` r
alphadiv <- left_join(richness, sample.tab %>% rownames_to_column(., var = "DNA_ID")) 
```

    ## Joining, by = "DNA_ID"

<img src="phyloseq_vtrandafir_files/figure-gfm/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

Boxes represent the 1.5 interquartile range, with the internal solid
line representing the median. Circles represent data points. p-values
are reported the non-parametric two sample Wilcoxon test, which tests
whether the means between two groups are equal (ns: p \> 0.05, \* : p≤
0.05, \*\* : p ≤ 0.01).

Difference in the alpha diversity indexes among conditions were tested
using pairwise Wilcoxon tests; p \< 0.05 was considered the threshold
significance for a difference between conditions.

From this plot we can see within the treatments that the richness (via
Chao index) of our samples significantly changed, while overall
diversity (via Shannon index) did not change. This suggest that while
richness decreased in both the control and ash leachate treatments, the
eveness was similar between the initial and final conditions.

<img src="phyloseq_vtrandafir_files/figure-gfm/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

From this plot we can see between the treatments that the richness of
the control samples were higher at the initial condition than the ash
leachate, suggesting that there may have been some quality control
issues as we would expect the initial samples to all have the same
richness. By timepoint 6, it looks like the richness was about the same
between the control and the ash leachate. Overall diversity was similar
between the treatments at the initial condition, but not by the end of
the experiment. The ash leachate samples at timepoint 6 may have been
less even.

Q6: Summarize the major takeaways from the alpha diversity plots you
generated.

A6: The richness of the control was approximately the same or slightly
higher than the Ash Leachate treatment at the beginning of the
experiment. By the end of the experiment, the Ash Leachate treatment
appears to have decreased in richness.

# Who??

Which taxa were important? Which taxa were contributing to the change in
community compositon?

**Note: Recovered 16S rRNA gene copy numbers do not equal organism
abundance.**

That said, we can generate a heat map of our samples showing us how the
relative abundance of different taxonomic groups change…potentially
giving us a visual of which taxa are most important to the alpha and
beta diversity patterns we observed. First, we’re going to generate a
custom table that will be easier to work with than a phyloseq object.

## Generate relative abundances

Our data currently shows number gene copies recovered, so we’ll convert
to percentages (relative abundances)

``` r
ps_std <- transform_sample_counts(ps_min, function(x) x/sum(x))
#extract the relative abundance table and coerce into dataframe
ps_std.tab <- as(otu_table(ps_std), "matrix")
ps_std.df = as.data.frame(ps_std.tab) 
```

## Make table

``` r
#first coerce the taxa table into a data frame
tax.df <-  as.data.frame(tax.tab) 
#then combine the data frames
custom.tab <- tax.df %>% 
  rownames_to_column(., var = "asv") %>% 
  left_join(., ps_std.df %>% rownames_to_column(., var = "asv")) %>% 
  #create a new index of that combines the  class, order, family, and genus values, you can play around here!!
  mutate(#pcofg = paste(Phylum, "_", Class, "_", Order,"_", Family, "_", Genus),
         # pcof = paste(Phylum, "_", Class, "_", Order,"_", Family,),
         pco = paste(Phylum, "_", Class, "_", Order)) %>% 
  select(-c(asv:Genus)) %>% 
  # select(pcof,everything()) %>% 
  # group_by(pcof) %>% 
  select(pco,everything()) %>% 
  group_by(pco) %>% 
  #here we are combining the relative abundances based on our grouping
  summarise_at(vars(contains(c("ASH171", "ASH172"))), sum, na.rm = T) %>% 
  ungroup()
```

    ## Joining, by = "asv"

``` r
#save the row names and then make them into the column names
colnames <- custom.tab[,1] 

#transpose the dataframe so we can merge with the sample info table
t_custom.tab <-  as.data.frame(t(custom.tab[,-1]))
# colnames(t_custom.tab) <- colnames$pcof
colnames(t_custom.tab) <- colnames$pco

#merge
sweet.tab <- t_custom.tab %>% 
  rownames_to_column(., var = "sample") %>% 
  left_join(., sample.tab %>% rownames_to_column(., var = "sample") %>% select(sample, Experiment, Location, Bottle, Treatment, Timepoint, days, cells)) %>% 
  select(sample, Experiment:cells, everything())
```

    ## Joining, by = "sample"

``` r
relabund <- sweet.tab %>% 
  select(-c(sample:cells)) %>% 
  #remove groups that are completely absent
  .[ , colSums(.) > 0] %>% 
  #arrange by biggest contributors
  .[, order(colSums(-.))] %>% 
  bind_cols(sweet.tab %>% select(sample:cells), .)
```

## Heatmap

<img src="phyloseq_vtrandafir_files/figure-gfm/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />

Q7: Who uniquely increased in relative abundance in the ash leachate
treatment that did not in the control? what about decrease?

A7: Alteromondales had a much larger increase in relative abundance in
the Ash Leachate treatment than in the control. SAR11 had a much larger
decrease in relative abundance in the Ash Leachate treatment than in the
control.

Everything shown here is just a snapshot of what you can look at with
your community composition data. There are many other resources you can
use to get ideas of how to look at different aspects of your data,
including the [phyloseq tutorial](https://joey711.github.io/phyloseq/)
and [happy belly bioinformatics](https://astrobiomike.github.io). It’s
up to you and your questions!!

# Save and knit

``` r
saveRDS(sweet.tab, "~/Documents/College/Fourth Year/EEMB 144L/Github/144l_students_2021/Output_Data/week7/Custom_ASV_Table.rds")
saveRDS(sub_ps, "~/Documents/College/Fourth Year/EEMB 144L/Github/144l_students_2021/Output_Data/week7/phyloseq_obj.rds")
saveRDS(ps_min, "~/Documents/College/Fourth Year/EEMB 144L/Github/144l_students_2021/Output_Data/week7/subsampled_phyloseq_obj.rds")
saveRDS(alphadiv, "~/Documents/College/Fourth Year/EEMB 144L/Github/144l_students_2021/Output_Data/week7/alphadiv.rds")
```

# Stacked Barplots

<img src="phyloseq_vtrandafir_files/figure-gfm/unnamed-chunk-21-1.png" style="display: block; margin: auto;" />
