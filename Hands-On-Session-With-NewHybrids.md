Hands On Session With NewHybrids
================

## Make You have the data

Marina has posted a folder of data. That folder (and its contents)
should be saved to a directory called `data` in your current working
directory in R.

## Let’s Investigate the Data

Johanna provided the data in “Variant Call Format.” This is a standard
format for genotypes derived from sequencing data. If you have never
looked at a VCF file before, it is worth doing. We can do that on the
Unix/Linux terminal like this:

``` sh
gzip -cd  data/HZ_ybch_diff95.vcf.gz  | less -S
```

You can get out of the `less` viewer by hitting the `q` key.

We can talk about this format—it can be daunting at first.

## Let us write some tools to get the data into NewHybrids format

There are a number of wonderful tools for manipulating VCF files on the
Unix command line (like `bcftools`); however, we can also use R, with
the package ‘vcfR’ to manipulate the data into shape for NewHybrids.
People here are probably more familiar with R that with bcftools, so we
will use R.

For a quick heads up on where we are going, we want to get the data into
the format that NewHybrids uses, which is described in [the newhybrids
manual](https://github.com/eriqande/newhybrids/blob/master/new_hybs_doc1_1Beta3.pdf),
and which looks like this:

    NumIndivs 79 
    NumLoci 49
    Digits 1
    Format Lumped

    LocusNames  sAAT1 sAAT2 sAAT3 ADA1 ADA2 ADH mAH1 mAH3 sAH mAAT1 CKA1 CKB FH
                bGLUA bGALA GDA2 G3PDH1 GPIB2 GPIA GR mIDHP1 mIDHP2 sIDHP2 LDHA1 LDHA2 LDHB1
                LDHB2 LDHC PEPB1 PEPLT sMDHA1 sMDHA2 sMDHB1 sMDHB2 sMEP2 mMEP1 MPI NTP PGDH
                PGK1 PGK2 PGM2 PEPD1 IDDH1 sSOD1 TPI1 TPI2 TPI4 TPI3

    1   11  11  11  0  11  11  11  11  11  11  11  11  11  32  11  11  21  11  11  11  11 
                 11  21  11  11  13  11  11  11  21  11  11  11  11  11  11  11  11  11  11  11 
                 11  11  11  11  11  11  11  11 
    2   21  11  21  11  11  12  11  11  12  11  11  11  0  11  11  11  12  11  12  11  11 
                 11  11  11  11  11  12  11  11  22  11  11  11  11  11  11  11  11  11  11  11 
                 11  11  11  11  11  11  11  11 
    3   11  11  21  12  11  21  11  11  11  11  11  11  11  22  11  11  31  11  11  11  11 
                 11  11  11  11  33  22  11  11  21  11  11  11  11  11  11  11  11  11  11  11 

So, we will get the data in an R data frame (a tidy tibble!) which we
can manipulate how we want and then we will have a function to punch out
the data set in the above format.

Notice that missing data is denoted by a 0.

### Get the data into a tidy tibble

So, we will get all the data a format that has columns: `chrom`, `pos`,
`indiv`, `geno`, where the `geno` column is 11 for `0/0` (i.e. the
reference homozygote), 22 for `1/1`, and either 12 or 21 for `0/1`, or
`1/0`, respectively. And 0 for a missing genotype.

Here is a function that does that

``` r
library(tidyverse)
library(vcfR)

#' Read a VCF file into a long-format tibble that will be easy to make into a newhybrids input
#' @param vcf_path the path to the VCF file
vcf2nh_tidy <- function(vcf_path) {
  v <- read.vcfR(vcf_path, verbose = FALSE)
  vt <- vcfR2tidy(v)
  ck <- vt$fix %>%
    distinct(ChromKey, CHROM)
  vt$gt %>%
    left_join(ck, by = join_by(ChromKey)) %>%
    rename(chrom = CHROM, pos = POS, indiv = Indiv, gt = gt_GT) %>%
    select(chrom, pos, indiv, gt) %>%
    mutate(
      geno = case_when(
        is.na(gt) ~ "0",
        gt == "0/0" ~ "11",
        gt == "1/1" ~ "22",
        gt == "0/1" | gt == "1/0" ~ "12",
        TRUE ~ "PROBLEM!!"
      )
    ) %>%
    select(-gt)
}
```

And we can use it to read in the two data sets that Johanna has: the HZ
one, which is focused on the hybrid zone birds, and the “all” one which
has them all.

``` r
all_long <- vcf2nh_tidy("data/all_ybch_diff85.vcf.gz")
hz_long <- vcf2nh_tidy("data/HZ_ybch_diff95.vcf.gz")
```

### Look over these samples and order them appropriately

We have individual sample names in these that look like:

``` r
unique(all_long$indiv) %>% head()
```

    ## [1] "Y_96N0042_CKDL240031457-1A_22GH2JLT4_L6"
    ## [2] "Y_96N0238_CKDL240031457-1A_22GH2JLT4_L6"
    ## [3] "Y_96N0349_CKDL240031457-1A_22GH2JLT4_L6"
    ## [4] "Y_97N5522_CKDL240031457-1A_22GH2JLT4_L6"
    ## [5] "Y_98N3098_CKDL240031457-1A_22GH2JLT4_L6"
    ## [6] "Y_98N3156_CKDL240031457-1A_22GH2JLT4_L6"

and

``` r
unique(hz_long$indiv) %>% head()
```

    ## [1] "YBCH1_23E29J01_S1" "YBCH3_23F01J01_S2" "YBCH5_23F01A02_S3"
    ## [4] "YBCH6_23F02A01_S4" "YBCH7_23F02J01_S5" "YBCH8_23F02J02_S6"

In order to view the results in an intuitive way, it will be nice to
arrange these birds in an order where all the putatively pure birds are
together and the putative hybrids are also near one another in the data
set. We can learn which birds are which from the “popmap” files for each
of them. These files look like:

###### popmap_allybch.csv

``` r
popmap_all <- read_csv("data/popmap_allybch.csv")
popmap_all
```

    ## # A tibble: 235 × 10
    ##    BGP_ID    Lat   Long State color   Population sample    pop      pop2    pop1
    ##    <chr>   <dbl>  <dbl> <chr> <chr>   <chr>      <chr>     <chr>   <dbl>   <dbl>
    ##  1 96N0042  41.1 -122.  CA    #004488 western    Y_96N004… p2    1.00e+0 1.01e-4
    ##  2 96N0238  36.3 -122.  CA    #004488 western    Y_96N023… p2    1.00e+0 1.02e-4
    ##  3 96N0349  40.9 -110.  UT    #004488 western    Y_96N034… p2    9.96e-1 4.41e-3
    ##  4 97N5522  40.8 -124.  CA    #004488 western    Y_97N552… p2    1.00e+0 1.01e-4
    ##  5 98N3098  42.3  -75.8 NY    #5aae61 eastern    Y_98N309… p1    1.01e-4 1.00e+0
    ##  6 98N3156  32.8  -80.0 SC    #5aae61 eastern    Y_98N315… p1    1.01e-4 1.00e+0
    ##  7 98N3248  42.6 -109.  WY    #004488 western    Y_98N324… p2    1.00e+0 1.02e-4
    ##  8 98N5070  39.0  -96.8 KS    #5aae61 eastern    Y_98N507… p1    1.00e-4 1.00e+0
    ##  9 98N5244  39.3  -94.9 KS    #5aae61 eastern    Y_98N524… p1    4.55e-2 9.54e-1
    ## 10 99N0449  42.5 -121.  OR    #004488 western    Y_99N044… p2    1.00e+0 1.02e-4
    ## # ℹ 225 more rows

###### popmap_hz.csv

``` r
popmap_hz <- read_csv("data/popmap_hz.csv")
popmap_hz
```

    ## # A tibble: 93 × 10
    ##    BGP_ID   Lat  Long State Color   Missingness_16milSNPs  pop1     pop2 sample 
    ##    <chr>  <dbl> <dbl> <chr> <chr>                   <dbl> <dbl>    <dbl> <chr>  
    ##  1 YBCH1   30.3 -97.6 TX    #5aae61               0.0179  1.00  0.0001   YBCH1_…
    ##  2 YBCH3   30.5 -95.6 TX    #5aae61               0.0283  0.966 0.0337   YBCH3_…
    ##  3 YBCH5   30.5 -95.6 TX    #5aae61               0.0845  0.971 0.0290   YBCH5_…
    ##  4 YBCH6   30.6 -95.6 TX    #5aae61               0.0625  1.00  0.000100 YBCH6_…
    ##  5 YBCH7   30.5 -95.7 TX    #5aae61               0.0823  0.895 0.105    YBCH7_…
    ##  6 YBCH8   30.6 -95.6 TX    #5aae61               0.0578  0.940 0.0597   YBCH8_…
    ##  7 YBCH11  34.1 -94.6 OK    #5aae61               0.0496  0.971 0.0292   YBCH11…
    ##  8 YBCH12  34.1 -94.6 OK    #5aae61               0.0202  0.965 0.0354   YBCH12…
    ##  9 YBCH14  34.5 -95.3 OK    #5aae61               0.0280  0.973 0.0273   YBCH14…
    ## 10 YBCH15  34.5 -95.4 OK    #5aae61               0.00641 1.00  0.000100 YBCH15…
    ## # ℹ 83 more rows
    ## # ℹ 1 more variable: pop <chr>

It appears that the `sample` column in both cases gives the sample names
that are used in the VCF files.

### A function to write the NewHybrids file

We now have the ingredients we need to write samples to a NewHybrids
file, sorted in a way that makes sense. One of the silly things about
NewHybrids is that the input file takes a sample number (in order from
1, 2, 3,…), rather than a sample name, so we will also want to write out
a file that gives the index/number that corresponds to each sample.

Let’s write a function to do this, so that we can use it over again, as
needed.

``` r
#' Write out a newhybrids file
#' 
#' @param ordered_samples  a tibble with a column `sample` that are in the order
#' you want individuals to appear in the newhybrids file. Note that if you want to
#' remove individuals, you just make sure they do not appear in `ordered_samples`
#' @param nh_tidy a tibble that is like the output of `vcf2nh_tidy()`.  All samples
#' requested in `ordered_samples` must appear in the `indiv` column.  If there is a
#' column `ZS` in nh_tidy, then its contents will for the z and s options string for
#' each row.
#' @param outpath  the path to the newhybrids output file to write.  New directories will
#' be created as needed and the file giving index/number corresonding to each
#' sample will be written out as `outpath + _sample_index`.
#' @output This is called for its side effect of writing two files that can be
#' used in NewHybrids.
write_nh <- function(ordered_samples, nh_tidy, outpath) {
  levs <- ordered_samples$sample
  
  if("ZS" %in% names(ordered_samples)) {
    os <- ordered_samples %>%
      select(sample, ZS)
  } else {
    os = ordered_samples %>%
      select(sample)
  }
  wide <- os %>%
    left_join(nh_tidy, by = join_by(sample == indiv)) %>%
    mutate(
      chrompos = str_c(chrom, ":", pos)
    ) %>%
    select(-chrom, -pos) %>%
    pivot_wider(
      names_from = chrompos,
      values_from = geno
    ) %>%
    mutate(idx = 1:n(), .before = sample)
  
  # write that out
  dn <- dirname(outpath)
  if(dn != ".") {
    dir.create(dn, showWarnings = FALSE, recursive = TRUE)
  }
  outNH <- outpath
  outIDX <- str_c(outpath, "_sample_index")
  
  Locus_vec <- names(wide) %>%
    setdiff(c("idx", "ZS", "sample"))
  
  # put the preamble/header stuff in the file
  cat(
    "NumIndivs ", nrow(wide),
    "\nNumLoci ", length(Locus_vec),
    "\nDigits 1\nFormat Lumped",
    "\n\n",
    "LocusNames ", paste(Locus_vec, sep = " ", collapse = " "),
    "\n\n",
    file = outNH,
    sep = ""
  )
  
  # write the data themselves
  wide %>%
    select(-sample) %>%
    write_tsv(outNH, col_names = FALSE, append = TRUE)
  
  # then write the indices
  wide %>%
    select(idx, sample) %>%
    write_tsv(file = outIDX)
}
```

## Running NewHybrids

Now we can make input files for NewHybrids and run it.

### No “prior” information

First, we throw all the individuals into the analysis without giving
newhybrids any hints about who is “pure” and who is from from the hybrid
zone, etc.

We will make each data set in its own directory. These are “vanilla”
runs because we are not using any information about pure and likely
hybrid birds.

``` r
popmap_all %>%
  arrange(desc(pop1)) %>%
  write_nh(all_long, "results/all/vanilla/dat.txt")

popmap_hz %>%
  arrange(desc(pop1)) %>%
  write_nh(all_long, "results/hz/vanilla/dat.txt")
```

And now we can run those.

With many loci and not many samples, NewHybrids does not converge to the
right part of the space and nearly everyone is inferred to be an F2!
Looking at the data, there is clear structure, it’s just that the MCMC
doesn’t get there.

### Pre-assign birds to different pure categories based on structure results

By giving NewHybrids a few hints about who is an eastern and who is a
western bird, the Markov chain converges quite quickly to a reasonable
looking part of the space. This is done by setting the Z option for each
individual.

We can set the Z for individuals with an admixture fraction \>0.999 or
less than \<0.0002, and see what happens then.

``` r
popmap_all %>%
  arrange(desc(pop1)) %>%
  mutate(
    ZS = case_when(
      pop1 > 0.999 ~ "z0",
      pop1 < 0.0002 ~ "z1",
      TRUE ~ ""
    )
  ) %>%
  write_nh(all_long, "results/all/a_few_Zs/dat.txt")
```

### Running newhybrids

This has to be done in the terminal.

After running that and looking at the results, we have what appears to
be some contamination, maybe. Suspect individuals are:

``` r
contam_suspects <- c(20, 23, 24, 39)
```

We will want to kick them out.

Also let us apply the s option to only the ones in the hybrid zone.
