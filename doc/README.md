CUPLR: Cancer of Unknown Primary Location Resolver
================

CUPLR is a classifier for identifying the primary tumor location of
metastatic tumors based on features derived from whole genome sequencing
data. The model itself is an ensemble of binary random forest
classifiers which each predict the probability of one cancer type. The
final predicted cancer type is the one with the highest probability.

The features used by CUPLR are extracted from the output of the [HMF
pipeline](https://github.com/hartwigmedical/pipeline5), specifically
from: (i) the somatic VCF file containing SBS, DBS and indel mutations,
(ii)
[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple)
output, and (iii)
[LINX](https://github.com/hartwigmedical/hmftools/tree/master/sv-linx)
output. Below is a summary of the features.

| Feature type | Source | Description                                                                               |
|--------------|--------|-------------------------------------------------------------------------------------------|
| sigs         | VCF    | SBS, DBS and indel [mutational signatures](https://cancer.sanger.ac.uk/signatures/)       |
| mut\_load    | VCF    | Total number of SBSs, DBSs and indels                                                     |
| rmd          | VCF    | Regional mutational density signatures                                                    |
| kataegis     | VCF    | Counts of kataegis foci                                                                   |
| chrom\_arm   | PURPLE | Chromosome arm gains/losses compared to the overall genome ploidy                         |
| gender       | PURPLE | Gender as derived from copy number data                                                   |
| gene         | LINX   | Deep deletions, amplifications, biallelic losses and mutations of cancer associated genes |
| sv           | LINX   | Simple and complex structural variants                                                    |
| fusion       | LINX   | Presence of gene fusions                                                                  |
| viral\_ins   | LINX   | Presence of viral insertions                                                              |

# Installation

## Dependencies

CUPLR depends on several R packages for basic functionality. These can
be installed as follows:

``` r
install.packages(c('randomForest','reshape2'))
install.packages('rfFC', repos='http://R-Forge.R-project.org')

install.packages('devtools')
devtools::install_github('https://github.com/UMCUGenetics/mutSigExtractor/')
```

Certain functions require other R packages (e.g. ggplot2 for plotting).
For details on optional dependencies, see the DESCRIPTION files.

## Loading CUPLR

CUPLR is composed of the following R packages:

-   `cuplr`: random forest training and prediction
-   `featureExtractor`: extract all features used by CUPLR
-   `statsExtra`: statistics used for univariate feature selection for
    training CUPLR
-   `nmf`: non-negative matrix factorization for generating regional
    mutational density signatures

To download CUPLR, run the following commands in the terminal,:

    cd /working/dir/
    git clone https://github.com/UMCUGenetics/cuplr/

Then load the `featureExtractor`, `cuplr`, and `mutSigExtractor`
packages. The remaining CUPLR packages (`statsExtra` and `nmf`) are only
needed for training CUPLR.

``` r
cuplr_dir <- '/path/to/cuplr/'
devtools::load_all(paste0(cuplr_dir,'/featureExtractor/')
devtools::load_all(paste0(cuplr_dir,'/cuplr/')
```

# Using CUPLR

This tutorial will demonstrate how cancer type prediction can be
performed using example input files located at `doc/data/`. These are 3
primary tumor samples from the PCAWG consortium, whose BAM files were
analyzed with the HMF pipeline.

## Extracting features

Extraction of features per sample is performed using
`extractFeaturesCuplr()`.

### Method 1

Paths to the required input files can be provided to `input.paths` as a
named character vector. Make sure that this vector has the names shown
below. The below code will return a 1-row dataframe with the values for
each feature.

``` r
sample_dir <- paste0(cuplr_dir,'/doc/data/DO48977/')
features <- extractFeaturesCuplr(
   input.paths=c(
      purple.smnv=        paste0(sample_dir,'/DO48977T.purple.somatic.vcf.gz'),   ## *.purple.somatic.vcf.gz
      purple.cnv=         paste0(sample_dir,'/DO48977T.purple.cnv.somatic.tsv'),  ## *.purple.cnv.somatic.tsv
      purple.purity=      paste0(sample_dir,'/DO48977T.purple.purity.tsv'),       ## *.purple.purity.tsv
      linx.drivers=       paste0(sample_dir,'/DO48977T.linx.driver.catalog.tsv'), ## *.linx.driver.catalog.tsv
      linx.fusions=       paste0(sample_dir,'/DO48977T.linx.fusion.tsv'),         ## *.linx.fusion.tsv
      linx.viral.inserts= paste0(sample_dir,'/DO48977T.linx.viral_inserts.tsv'),  ## *.linx.viral_inserts.tsv
      linx.vis.sv.data=   paste0(sample_dir,'/DO48977T.linx.vis_sv_data.tsv')     ## *.linx.vis_sv_data.tsv
   )
)
features[,1:5]
```

    ##   sigs.SBS1 sigs.SBS2 sigs.SBS3 sigs.SBS4 sigs.SBS5
    ## 1 0.1765092  0.047903         0         0 0.0634117

### Method 2

Alternatively, the path to a folder with all the required input files
can be provided to `in.dir`. The filenames should end with those shown
in the comments and the above example (e.g. \*.purple.somatic.vcf.gz).

An output folder can optionally be specified to `out.dir`. With this
approach, intermediate files are also written and allows resuming the
feature extraction if a crash occurs (instead of starting the extraction
from the beginning).

For this tutorial, we will use `in.dir` and `out.dir` to load the
example input data and write the outputs. The below code extracts the
features for all 3 example samples.

``` r
sample_dirs <- list.dirs(paste0(cuplr_dir,'doc/data/'), recursive=F, full.names=T)

for(in_dir in sample_dirs){
   message('\nExtracting features at: ', in_dir)

   out_dir <- paste0(in_dir,'/output/')
   dir.create(out_dir, showWarnings=F)

   extractFeaturesCuplr(
      in.dir=in_dir, ## Provide `input` to in.dir instead of `input.paths`
      out.dir=out_dir,
      verbose=1
   )
}
```

We can then read the `all_features.txt.gz` files which are the 1-row
dataframes (as mentioned above). It is important to specify
`check.names=F` to `read.delim()` as the column names have some illegal
characters which would otherwise be modified, and we don’t want this to
happen.

``` r
out_paths <- paste0(sample_dirs,'/output/all_features.txt.gz')

## Read txt files
features <- lapply(out_paths, function(i){ 
   read.delim(i, check.names=F) 
})

## Merge 1-row dataframes into a single dataframe
features <- do.call(rbind, features)

## Assign sample names to rownames
sample_names <- basename(sample_dirs)
rownames(features) <- sample_names
features[,1:5]
```

    ##           sigs.SBS1  sigs.SBS2 sigs.SBS3 sigs.SBS4  sigs.SBS5
    ## DO220848 0.06598833 0.04049707         0         0 0.02095571
    ## DO48977  0.17650919 0.04790300         0         0 0.06341170
    ## DO51095  0.00000000 0.00000000         0         0 0.00000000

### Summary of feature names

There are 4000+ features; too many to print on screen. However, with
below code we can see the names of the first couple of features in each
group.

``` r
feature_groups <- groupFeaturesByTag(colnames(features), rm.tags=F)
lapply(feature_groups, head)
```

    ## $sigs
    ## [1] "sigs.SBS1" "sigs.SBS2" "sigs.SBS3" "sigs.SBS4" "sigs.SBS5" "sigs.SBS6"
    ## 
    ## $mut_load
    ## [1] "mut_load.snv"   "mut_load.indel" "mut_load.dbs"  
    ## 
    ## $rmd
    ## [1] "rmd.1p_1" "rmd.1p_2" "rmd.1p_3" "rmd.1p_4" "rmd.1p_5" "rmd.1p_6"
    ## 
    ## $sv
    ## [1] "sv.n_events"                "sv.complex.n_events"       
    ## [3] "sv.complex.largest_cluster" "sv.double_minutes"         
    ## [5] "sv.foldbacks"               "sv.LINEs"                  
    ## 
    ## $chrom_arm
    ## [1] "chrom_arm.gain.1p" "chrom_arm.gain.1q" "chrom_arm.gain.2p"
    ## [4] "chrom_arm.gain.2q" "chrom_arm.gain.3p" "chrom_arm.gain.3q"
    ## 
    ## $gene
    ## [1] "gene.ACVR1.amp"       "gene.ACVR1.deep_del"  "gene.ACVR1.biall"    
    ## [4] "gene.ACVR1.monoall"   "gene.ACVR2A.amp"      "gene.ACVR2A.deep_del"
    ## 
    ## $fusion
    ## [1] "fusion.@IGH_*"    "fusion.@IGL_BCL6" "fusion.@IGL_MYC"  "fusion.*_ABL1"   
    ## [5] "fusion.*_ALK"     "fusion.*_BRAF"   
    ## 
    ## $viral_ins
    ## [1] "viral_ins.AAV" "viral_ins.EBV" "viral_ins.HBV" "viral_ins.HCV"
    ## [5] "viral_ins.HIV" "viral_ins.HPV"
    ## 
    ## $gender
    ## [1] "gender.gender"
    ## 
    ## $kataegis
    ## [1] "kataegis.foci"

## Predicting cancer type

Now that we have the features, we can predict the cancer type. We first
need to load CUPLR itself, as well as the probability calibration
curves.

The scores outputted by a random forest need to be calibrated to yield
true probabilities. A well calibrated classifier should classify the
samples such that among the samples to which it gave a score close to
e.g 0.8, approximately 80% actually belong to the positive class. For
more info, see the [scikit-learn
documentation](https://scikit-learn.org/stable/modules/calibration.html)

``` r
model <- readRDS(MODEL_PATH)
prob_cal_curves <- read.delim(PROB_CAL_CURVES_PATH)
```

``` r
model <- cacheAndReadData(MODEL_PATH)
prob_cal_curves <- read.delim(PROB_CAL_CURVES_PATH)
```

We can then use `predict()` to get the probabilities of each cancer type
for each sample. The output of `predict()` is a list containing several
objects. Here we have specified `calc.feat.contrib=T` which will
calculate the contribution of each feature to each cancer type
prediction. This is required to generate the patient report (see next
section)

``` r
pred_report <- predict(
   object=model,
   newdata=features,
   prob.cal.curves=prob_cal_curves,
   calc.feat.contrib=T
)
pred_report
```

    ## Objects in list:
    ## $prob $prob_scaled $class_pred $feat_contrib
    ## 
    ## $prob_scaled
    ##              Biliary       Breast       Cervix CNS_Glioma CNS_Medullo
    ## DO220848 0.000000000 6.206318e-04 0.0000278392          0           0
    ## DO48977  0.001187102 7.016440e-04 0.8207116788          0           0
    ## DO51095  0.004693888 1.914623e-05 0.0000000000          0           0
    ##          CNS_PiloAstro   Colorectum Colorectum_NET      Gastric HeadAndNeck_ACC
    ## DO220848  0.000000e+00 0.0012737281   0.0000000000 0.0000000000               0
    ## DO48977   0.000000e+00 0.0026945400   0.0000000000 0.0001003129               0
    ## DO51095   6.998627e-05 0.0005451065   0.0007053233 0.0098180250               0
    ##          HeadAndNeck_Other Kidney       Liver     Lung_NSC Lung_SC     Lymphoid
    ## DO220848      0.0001201604      0 0.000000000 0.0057080125       0 3.880966e-05
    ## DO48977       0.0077212784      0 0.000000000 0.0063026494       0 1.052110e-04
    ## DO51095       0.0000000000      0 0.000255996 0.0001332586       0 4.366619e-02
    ##           Mesothelium Myeloid        Ovary Pancreas Pancreas_NET  Prostate
    ## DO220848 0.0000000000       0 5.013721e-04        0  0.001203551 0.0000000
    ## DO48977  0.0001084994       0 5.993926e-05        0  0.000000000 0.0000000
    ## DO51095  0.0179715590       0 0.000000e+00        0  0.000000000 0.8006185
    ##          Sarcoma_GIST Sarcoma_Leiomyo Sarcoma_Lipo Sarcoma_Other Skin_Melanoma
    ## DO220848            0    0.0006902753            0  0.0010987696             1
    ## DO48977             0    0.0000000000            0  0.0008740941             0
    ## DO51095             0    0.0000000000            0  0.0006535543             0
    ##           Skin_Other     Thyroid   Urothelial       Uterus
    ## DO220848 0.004528138 0.003746576 3.852320e-04 0.000000e+00
    ## DO48977  0.000000000 0.010817591 1.757926e-02 2.099301e-05
    ## DO51095  0.003908766 0.000000000 2.984463e-05 0.000000e+00

However, the above raw prediction output is not informative at a glance.
We can use `summary()` to show the prediction in a neat table.

``` r
summary(
   pred_report, 
   top.n.classes=3, ## Optional: Shows the top 3 classes and their probabilities
   top.n.feat=3     ## Optional: Shows the top 3 features contributing to the top class prediction
)
```

    ##     sample  pred_class.1 pred_class.2 pred_class.3 pred_prob.1 pred_prob.2
    ## 1 DO220848 Skin_Melanoma     Lung_NSC   Skin_Other       1.000       0.006
    ## 2  DO48977        Cervix   Urothelial      Thyroid       0.821       0.018
    ## 3  DO51095      Prostate     Lymphoid  Mesothelium       0.801       0.044
    ##   pred_prob.3                   feat.1             feat.2            feat.3
    ## 1       0.005 rmd.Skin_Melanoma.1=0.17    sigs.SBS7=0.135  sigs.SBS38=0.098
    ## 2       0.011      viral_ins.HPV=0.305 rmd.Cervix.1=0.279  sigs.SBS13=0.021
    ## 3       0.018 fusion.TMPRSS2_ERG=0.271  mut_load.snv=0.18 mut_load.dbs=0.07

## Graphical patient report

Using `patientReport()` we can plot the class probabilities and feature
contributions per class. This is primarily useful for clinical
reporting.

``` r
patient_report <- patientReport(
   prob=pred_report$prob_scaled, 
   feat.contrib=pred_report$feat_contrib, 
   sample.name='DO51095'
)
```

``` r
grid::grid.draw(patient_report)
```

![](/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/doc/supp/patient_report_example.png)

We can see that sample DO51095 is a prostate cancer sample, with and
this prediction is supported by this sample having a TMPRSS2-ERG fusion,
a well-known prostate cancer event.

In cases where the probabilities are more uncertain, more feature
contribution panels corresponding to the top predicted classes will be
shown, which will further aid in determining the cancer type.
