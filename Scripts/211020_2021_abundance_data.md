---
output:
  pdf_document: default
  html_document: default
---
211020_2021_Abundance
================
Lynette Lee
10/20/2021

# Goal

Use the data collected of the inoculum from the flow cytometry to
determine the growth curves of the bacterial communities over different
treatments. The treatments include the control, which are filtered sea
water; kelp exudate, which is the dissolved organic carbon (DOC)
released from the macroalgae; kelp exudate supplemented with addition
nutrients, which are nitrogen and phosphorous; and nutrient enriched sea
water, where glucose, nitrogen, and phosphorous was added to the media.

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.5     ✓ dplyr   1.0.7
    ## ✓ tidyr   1.1.4     ✓ stringr 1.4.0
    ## ✓ readr   2.0.2     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(readxl)
library(lubridate)
```

    ## 
    ## Attaching package: 'lubridate'

    ## The following objects are masked from 'package:base':
    ## 
    ##     date, intersect, setdiff, union

# Import Data

``` r
excel_sheets("~/Desktop/GitHub/EEMB 144L/211019_2021_Abundance/144l_students_2021/Input_Data/week4/144L_2021_BactAbund.xlsx")
```

    ## [1] "Metadata" "Data"

``` r
metadata <- read_excel("~/Desktop/GitHub/EEMB 144L/211019_2021_Abundance/144l_students_2021/Input_Data/week4/144L_2021_BactAbund.xlsx", sheet = "Metadata")
#glimpse(metadata)

unique(metadata$Bottle)
```

    ## [1] "A" "B" "C" "D" "E" "F" "G" "H"

``` r
unique(metadata$Treatment)
```

    ## [1] "Control"                        "Kelp Exudate"                  
    ## [3] "Kelp Exudate_Nitrate_Phosphate" "Glucose_Nitrate_Phosphate"

``` r
data <- read_excel("~/Desktop/GitHub/EEMB 144L/211019_2021_Abundance/144l_students_2021/Input_Data/week4/144L_2021_BactAbund.xlsx", sheet = "Data")
#glimpse(data)

joined <- left_join(metadata, data)
```

    ## Joining, by = c("Bottle", "Timepoint")

``` r
glimpse(joined)
```

    ## Rows: 80
    ## Columns: 19
    ## $ Experiment           <chr> "144L_2021", "144L_2021", "144L_2021", "144L_2021…
    ## $ Location             <chr> "Goleta Pier", "Goleta Pier", "Goleta Pier", "Gol…
    ## $ Temperature          <dbl> 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 1…
    ## $ Depth                <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    ## $ Bottle               <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",…
    ## $ Timepoint            <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6…
    ## $ Treatment            <chr> "Control", "Control", "Control", "Control", "Cont…
    ## $ Target_DOC_Amendment <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ Inoculum_L           <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2…
    ## $ Media_L              <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5…
    ## $ Datetime             <chr> "2021-10-04T16:00", "2021-10-05T08:00", "2021-10-…
    ## $ TOC_Sample           <lgl> TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FA…
    ## $ Cell_Sample          <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, T…
    ## $ DAPI_Sample          <lgl> TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FA…
    ## $ DNA_Sample           <lgl> TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FA…
    ## $ Nutrient_Sample      <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, F…
    ## $ all_cells_uL         <chr> "901.48904752420799", "4302.9300548457404", "3944…
    ## $ LNA_cells_uL         <chr> "653.047184033284", "768.27893058466896", "937.24…
    ## $ HNA_cells_uL         <chr> "248.441863490923", "3534.65112426107", "3007.701…

# Prepare Data

The date and time data were separated.

The cell concentrations units were converted from total cells per
microliters (uL) into total cells per liters (L).

``` r
cells <- joined %>%
  mutate(Datetime = ymd_hm(Datetime),
         all_cells_L = as.numeric(all_cells_uL) * 1000 * 1000,
         LNA_cells_L = as.numeric(LNA_cells_uL) * 1000 * 1000,
         HNA_cells_L = as.numeric(HNA_cells_uL) * 1000 * 1000) %>%
  group_by(Treatment, Bottle) %>%
  mutate(interv = interval(first(Datetime), Datetime),
         s = as.numeric(interv),
         hours = s/3600,
         days = hours/24) %>%
  ungroup() %>%
  select(Experiment:DNA_Sample, all_cells_L, hours, days) %>%
  drop_na(all_cells_L)
glimpse(cells)
```

    ## Rows: 60
    ## Columns: 18
    ## $ Experiment           <chr> "144L_2021", "144L_2021", "144L_2021", "144L_2021…
    ## $ Location             <chr> "Goleta Pier", "Goleta Pier", "Goleta Pier", "Gol…
    ## $ Temperature          <dbl> 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 1…
    ## $ Depth                <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    ## $ Bottle               <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", "B",…
    ## $ Timepoint            <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6, 7…
    ## $ Treatment            <chr> "Control", "Control", "Control", "Control", "Cont…
    ## $ Target_DOC_Amendment <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    ## $ Inoculum_L           <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2…
    ## $ Media_L              <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5…
    ## $ Datetime             <dttm> 2021-10-04 16:00:00, 2021-10-05 08:00:00, 2021-1…
    ## $ TOC_Sample           <lgl> TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FA…
    ## $ Cell_Sample          <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, T…
    ## $ DAPI_Sample          <lgl> TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FA…
    ## $ DNA_Sample           <lgl> TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FA…
    ## $ all_cells_L          <dbl> 901489048, 4302930055, 3944945700, 3467546844, 23…
    ## $ hours                <dbl> 0, 16, 28, 40, 52, 64, 76, 88, 100, 0, 16, 28, 40…
    ## $ days                 <dbl> 0.0000000, 0.6666667, 1.1666667, 1.6666667, 2.166…

# Plot Growth Curve

The growth curves of the microbial communities under different
treatments over the span of 5 days were graphed together. The microbial
growth peaked before the end of the 1st day for the control and kelp
exudate treatments, while microbial growth peaked during the 2nd day for
the supplemented kelp exudate and glucose-nitrate-phosphate treatments.
The control also had similar growth curve to the kelp exudate
treatmehnt.

There were some missing time points since some of the labels on the FCM
tubes fell off when storing them in liquid nitrogen. As a result, the
data contained some NA values which were hidden earlier.

``` r
custom.colors <- c("Control" = "#FAD7A0", "Kelp Exudate" = "#DAF7A6", "Kelp Exudate_Nitrate_Phosphate" = "#A3E4D7", "Glucose_Nitrate_Phosphate" = "#D7BDE2")

levels <- c("Control", "Kelp Exudate", "Kelp Exudate_Nitrate_Phosphate", "Glucose_Nitrate_Phosphate")

cells %>%
  mutate(dna = ifelse(DNA_Sample == T, "*", NA)) %>%
  ggplot(aes(x=days, y=all_cells_L, group = interaction(Treatment, Bottle))) +
  geom_line(aes(color = factor(Treatment, levels = levels)), size = 1) + 
  geom_point(aes(fill = factor(Treatment, levels = levels)), size = 3, color = "black", shape = 21) +
  geom_text(aes(label = dna), size = 10, color = "#F5B7B1") + 
  labs(x = "Days", y = expression(paste("Cells, L"^-1)), fill = " ") +
  guides(color = "none") +
  scale_color_manual(values = custom.colors) +
  scale_fill_manual(values = custom.colors) +
  theme_bw()
```

    ## Warning: Removed 40 rows containing missing values (geom_text).

![](211020_2021_abundance_data_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

# Calculations

To graph the natural log of the growth curve, the data was converted
into the natural log values of the microbial cell concentrations. The
incremental increase in cell concentration was found by subtracting the
cellular concentration of Tn by the cellular concentration of Tn-1,
which is the earlier time point.

``` r
ln_cells <- cells %>%
  group_by(Treatment, Bottle) %>%
  mutate(ln_cells = log(all_cells_L),
         diff_ln_cells = ln_cells - lag(all_cells_L, default = first(ln_cells)))
```

# Plotting Logarithmic Data

The logarithmic data was graphed and the data was separated by the
different treatments and the bottles. The data shows the change in
cellular concentration ovaer time. The exponential growth of the marine
microbial communities appear to peak between the first treatment and day
1 for most treatment bottles. However, the kelp exudate with nitrate and
phosphate treatment appears to have later exponential growth in between
day 1 and 2.

``` r
ln_cells %>%
  mutate(dna = ifelse(DNA_Sample == T, "*", NA)) %>%
  ggplot(aes(x=days, y=diff_ln_cells, group = interaction(Treatment, Bottle))) +
  geom_line(aes(color = factor(Treatment, levels = levels)), size = 1) +
  geom_point(aes(fill = factor(Treatment, levels = levels)), size = 3, color = "black", shape = 21) +
  geom_text(aes(label = dna), size = 10, color = "#F5B7B1") +
  labs(x = "Days", y = expression(paste("∆ln cells, L"^-1)), fill = " ") +
  guides(color = "none") +
  scale_color_manual(values = custom.colors) +
  scale_fill_manual(values = custom.colors) +
  facet_wrap("Bottle", ncol = 2) +
  theme_bw()
```

    ## Warning: Removed 40 rows containing missing values (geom_text).

![](211020_2021_abundance_data_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

# Identifying Exponential Growth

To graph the logarithmic graph of the cellular concentrations over time,
the natural log values of the cellular concentrations were plotted. This
graph better shows when the exponential growth occurred.

Based on this data, the control and the kelp exudate treatment
experienced exponential growth between days 0-1. The kelp exudate
supplemented with nutrient treatment had one of its bottle showing
exponential growth in between days 0 and 1, however one of its bottle
has inconclusive data due to the NA value at T0. In addition, the
glucose-nitrate-phosphate treatment also appears to have one of its data
missing a T0 value while the other bottle displays exponential growth
between T0 and T1.

``` r
ln_cells %>%
  mutate(dna = ifelse(DNA_Sample == T, "*", NA)) %>%
  ggplot(aes(x=days, y=ln_cells, group = interaction(Treatment, Bottle))) +
  geom_line(aes(color = factor(Treatment, levels = levels)), size = 1) +
  geom_point(aes(fill = factor(Treatment, levels = levels)), size = 3, color = "black", shape = 21) +
  geom_text(aes(label = dna), size = 10, color = "#F5B7B1") + 
  labs(x = "Days", y = expression(paste("ln cells, L"^-1)), fill = " ") +
  guides(color = "none") +
  scale_color_manual(values = custom.colors) +
  scale_fill_manual(values = custom.colors) +
  facet_wrap("Bottle", ncol = 2) +
  theme_bw()
```

    ## Warning: Removed 40 rows containing missing values (geom_text).

![](211020_2021_abundance_data_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

# Conclusion

The RStudio program efficiently gathers data from multiple sources and
combines them, which was observed when the “metadata” was joined to
“data.” Then the newly joined data could be rearranged so that it can be
later used to calculate new values and plot graphs. For example, when
preparing the data for calculations, the datetime data points were
revised so that the actual dates became seaparate from the time values.
Once this new data was collected and labeled as “cells,” this new data
set was used to calculate natural log values to create logarithmic
graphs. These logarithmic graphs better display the exponential growth
of the marine microbial communities than simply observing change in
cellular concentration.
