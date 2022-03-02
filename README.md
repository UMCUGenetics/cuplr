CUPLR: Cancer of Unknown Primary Location Resolver
================

CUPLR is a classifier for identifying the primary tumor location of
metastatic tumors based on features derived from whole genome sequencing
data. The model itself is an ensemble of binary random forest
classifiers which each predict the probability of one cancer type. The
final predicted cancer type is the one with the highest probability.

The features used by CUPLR are extracted from the output of the [HMF
pipeline 5](https://github.com/hartwigmedical/pipeline5), specifically
from: (i) the somatic VCF file containing SBS, DBS and indel mutations,
(ii)
[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple)
(v2.53) output, and (iii)
[LINX](https://github.com/hartwigmedical/hmftools/tree/master/linx)
(v1.17) output. Below is a summary of the features.

| Feature type | Source | Description                                                                                                   |
|--------------|--------|---------------------------------------------------------------------------------------------------------------|
| sigs         | VCF    | Relative contribution of SBS, DBS and indel [mutational signatures](https://cancer.sanger.ac.uk/signatures/)  |
| mut_load     | VCF    | Total number of SBSs, DBSs and indels                                                                         |
| rmd          | VCF    | [Regional mutational density (RMD)](https://www.nature.com/articles/nature14221) profiles as extracted by NMF |
| chrom_arm    | PURPLE | Chromosome arm gains/losses compared to the overall genome ploidy                                             |
| gender       | PURPLE | Gender as derived from copy number data                                                                       |
| genome       | PURPLE | Genome properties, including overall ploidy, diploid proportion, and presence of whole genome duplication     |
| gene         | LINX   | Deep deletions, amplifications, biallelic losses and mutations of cancer associated genes                     |
| sv           | LINX   | Simple and complex structural variants                                                                        |
| fusion       | LINX   | Presence of gene fusions                                                                                      |
| viral_ins    | LINX   | Presence of viral sequence insertions                                                                         |

CUPLR was trained with tumor samples from 6756 patients from the Hartwig
Medical Foundation (HMF) and the Pan-Cancer Analysis of Whole Genomes
(PCAWG) consortium. The model can predict the primary tumor location
with an accuracy of 0.90 based on cross-validation and 0.89 based on
predictions on a held out test set.

For details on performance, the top features used by CUPLR, and other
details, please see the plots at `doc/perf/`.

# Reference

**Machine learning-based tissue of origin classification for cancer of
unknown primary diagnostics using genome-wide mutation features**  
*Luan Nguyen, Arne Van Hoeck, Edwin Cuppen*  
<https://www.biorxiv.org/content/10.1101/2021.10.05.463244v1.full>

For now, please refer to the paper on bioRxiv. Once the paper is
published, this link will be updated.

# Installation

## Dependencies

CUPLR depends on several R packages for basic functionality. These can
be installed as follows:

``` r
install.packages(c('randomForest','reshape2','ggplot2','ggrepel'))

## For calculating random forest feature contributions:
install.packages('rfFC', repos='http://R-Forge.R-project.org')

## For reading VCFs and extracting mutational signatures:
install.packages('devtools')
devtools::install_github('https://github.com/UMCUGenetics/mutSigExtractor/')
```

Certain functions require other R packages. For details on optional
dependencies, please see the DESCRIPTION files.

## Loading CUPLR

CUPLR is composed of the following R packages:

-   `cuplr`: random forest training, prediction, and performance
-   `featureExtractor`: extract the features used by CUPLR
-   `statsExtra`: statistics used for univariate feature selection for
    training CUPLR
-   `nmf`: wrapper around the NNLM package for non-negative matrix
    factorization for generating RMD profiles.

To download CUPLR, run the following commands in the terminal:

    cd /working/dir/
    git clone https://github.com/UMCUGenetics/cuplr/

Then load the `featureExtractor`, `cuplr`, and `mutSigExtractor`
packages. The remaining CUPLR packages (`statsExtra` and `nmf`) are only
needed for training CUPLR.

``` r
cuplr_dir <- '/path/to/cuplr/'
devtools::load_all(paste0(cuplr_dir,'/featureExtractor/'))
devtools::load_all(paste0(cuplr_dir,'/cuplr/'))
```

# Using CUPLR

This tutorial will demonstrate how cancer type prediction can be
performed using example input files located at `doc/data/`. These are 3
primary tumor samples from the PCAWG consortium, whose BAM files were
analyzed with the HMF pipeline.

## Extracting features

Extraction of features per sample is performed using
`extractFeaturesCuplr()`. This typically takes about 15-30 seconds per
sample.

### Method 1

Paths to the required input files can be provided to `input.paths` as a
named character vector. Make sure that this vector has the names shown
below. The below code will return a 1-row dataframe with the values for
each feature.

``` r
sample_dir <- paste0(cuplr_dir,'/doc/data/DO35886/')
features <- extractFeaturesCuplr(
   input.paths=c(
      purple.smnv=        paste0(sample_dir,'/DO35886T.purple.somatic.vcf.gz'),   ## *.purple.somatic.vcf.gz
      purple.cnv=         paste0(sample_dir,'/DO35886T.purple.cnv.somatic.tsv'),  ## *.purple.cnv.somatic.tsv
      purple.purity=      paste0(sample_dir,'/DO35886T.purple.purity.tsv'),       ## *.purple.purity.tsv
      linx.drivers=       paste0(sample_dir,'/DO35886T.linx.driver.catalog.tsv'), ## *.linx.driver.catalog.tsv
      linx.fusions=       paste0(sample_dir,'/DO35886T.linx.fusion.tsv'),         ## *.linx.fusion.tsv
      linx.viral.inserts= paste0(sample_dir,'/DO35886T.linx.viral_inserts.tsv'),  ## *.linx.viral_inserts.tsv
      linx.vis.sv.data=   paste0(sample_dir,'/DO35886T.linx.vis_sv_data.tsv')     ## *.linx.vis_sv_data.tsv
   )
)
features[,1:5]
```

    ##    sigs.SBS1  sigs.SBS2   sigs.SBS3 sigs.SBS4 sigs.SBS5
    ## 1 0.04792343 0.01030294 0.009803203         0         0

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

    ##             sigs.SBS1  sigs.SBS2   sigs.SBS3 sigs.SBS4  sigs.SBS5
    ## DO220844   0.07132016 0.02757172 0.000000000         0 0.06116417
    ## DO35886    0.04792343 0.01030294 0.009803203         0 0.00000000
    ## HMF000185A 0.11327875 0.00000000 0.000000000         0 0.27081490

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
    ## [1] "sv.sv_load"           "sv.LINEs"             "sv.DEL_[0,300)"      
    ## [4] "sv.DEL_[300,1e+03)"   "sv.DEL_[1e+03,1e+04)" "sv.DEL_[1e+04,1e+05)"
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
    ## [1] "fusion.@IGH_*"     "fusion.@IGL_MYC"   "fusion.*_ETV1/4/5"
    ## [4] "fusion.*_RET"      "fusion.*_TFE3"     "fusion.*_TFEB"    
    ## 
    ## $viral_ins
    ## [1] "viral_ins.AAV" "viral_ins.EBV" "viral_ins.HBV" "viral_ins.HCV"
    ## [5] "viral_ins.HIV" "viral_ins.HPV"
    ## 
    ## $genome
    ## [1] "genome.ploidy"             "genome.diploid_proportion"
    ## [3] "genome.has_wgd"           
    ## 
    ## $gender
    ## [1] "gender.gender"

## Predicting cancer type

Now that we have the features, we can predict the cancer type. We first
need to load CUPLR itself as well as the probability calibration curves.

The scores outputted by a random forest need to be calibrated to yield
true probabilities. This means that a probability of 0.9 should mean
that for example: amongst 100 hypothetical samples predicted as breast
cancer, 90 of them are actually breast cancer. For more info, see the
[scikit-learn
documentation](https://scikit-learn.org/stable/modules/calibration.html)
about this topic.

``` r
model <- readRDS(MODEL_PATH)
prob_cal_curves <- read.delim(PROB_CAL_CURVES_PATH)
```

We can then use `predict()` to get the probabilities of each cancer type
for each sample. The output of `predict()` is a list containing several
objects. Here we have specified `calc.feat.contrib=T` which will
calculate the contribution of each feature to each cancer type
prediction. This is required to generate the patient report (see next
section).

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
    ## Calibrated probabilities:
    ## $prob_scaled
    ##                 Biliary   Breast Cervix  CNS_Glioma CNS_Medullo CNS_PiloAstro
    ## DO220844   0.0005773561 0.213022      0 0.000000000 0.000000000     0.0000000
    ## DO35886    0.0000000000 0.000000      0 0.001242909 0.005507125     0.9135719
    ## HMF000185A 0.0000000000 0.000000      0 0.000000000 0.000000000     0.0000000
    ##              Colorectal Gastroesophageal HeadAndNeck_Other
    ## DO220844   0.0003462293     0.0002115484       0.001378098
    ## DO35886    0.0004806873     0.0016663170       0.000000000
    ## HMF000185A 0.0000000000     0.0000000000       0.000000000
    ##            HeadAndNeck_SalivaryGland Kidney_Chromophobe Kidney_ClearCell
    ## DO220844                0.0000000000       0.0000000000     0.0000000000
    ## DO35886                 0.0000000000       0.0007175474     0.0001646324
    ## HMF000185A              0.0004121185       0.0000000000     0.0001097549
    ##            Kidney_Papillary       Liver Lung_NonSmallCell Lung_SmallCell
    ## DO220844       0.000000e+00 0.000000000       0.003274019              0
    ## DO35886        0.000000e+00 0.001093323       0.004961213              0
    ## HMF000185A     9.299103e-05 0.000000000       0.001241159              0
    ##            Lymphoid Mesothelium      Myeloid NET_Gastrointestinal     NET_Lung
    ## DO220844          0           0 0.0000000000         1.507647e-04 0.000000e+00
    ## DO35886           0           0 0.0004481103         4.186282e-05 0.000000e+00
    ## HMF000185A        0           0 0.0000000000         1.116342e-04 6.116457e-05
    ##            NET_Pancreas      Ovarian     Pancreas Prostate Sarcoma_GIST
    ## DO220844   0.0004615711 0.0015382461 0.0011009668        0            0
    ## DO35886    0.0383353453 0.0001209926 0.0002253971        0            0
    ## HMF000185A 0.0392612195 0.0000000000 0.0004507942        0            0
    ##            Sarcoma_Leiomyo Sarcoma_Lipo Sarcoma_Osteo Sarcoma_Other
    ## DO220844      0.0003592354    0.0000000  0.0022210274  0.0004005747
    ## DO35886       0.0007545717    0.0000000  0.0031881031  0.0010086661
    ## HMF000185A    0.0023550295    0.5667544  0.0006195843  0.4688011203
    ##            Skin_Carcinoma Skin_Melanoma      Thyroid   Urothelial Uterus
    ## DO220844      0.000000000     0.8560547 0.0000000000 0.0003200161      0
    ## DO35886       0.000000000     0.0000000 0.0037040469 0.0004442942      0
    ## HMF000185A    0.003430203     0.0000000 0.0003721351 0.0007926048      0

However, the above raw prediction output is not informative at a glance.
We can use `summary()` to show the prediction in a neat table.

``` r
summary(
   pred_report, 
   top.n.classes=3, ## Optional: Shows the top 3 classes and their probabilities
   top.n.feat=3     ## Optional: Shows the top 3 features contributing to the top class prediction
)
```

    ##       sample  pred_class.1  pred_class.2      pred_class.3 pred_prob.1
    ## 1   DO220844 Skin_Melanoma        Breast Lung_NonSmallCell       0.856
    ## 2    DO35886 CNS_PiloAstro  NET_Pancreas       CNS_Medullo       0.914
    ## 3 HMF000185A  Sarcoma_Lipo Sarcoma_Other      NET_Pancreas       0.567
    ##   pred_prob.2 pred_prob.3           feat_contrib.1.1         feat_contrib.1.2
    ## 1       0.213       0.003  rmd.Skin_Melanoma.1=0.512         sigs.SBS38=0.267
    ## 2       0.038       0.006 fusion.KIAA1549_BRAF=0.214     mut_load.indel=0.165
    ## 3       0.469       0.039     fusion.FUS_DDIT3=0.325 rmd.Sarcoma_Lipo.2=0.319
    ##                   feat_contrib.1.3
    ## 1 sv.COMPLEX.largest_cluster=0.021
    ## 2                  sigs.ID12=0.082
    ## 3                  sigs.SBS8=0.017

## Graphical patient report

Using `patientReport()` we show plot the output of CUPLR graphically for
one patient. The left panel of the patient report shows the cancer type
probabilities. The right panels show the values of the most important
features contributing to each of the top predicted cancer types.

In the right panels, the feature value averages in patients with the
respective predicted cancer type and patients with other cancer types
are also shown to provide context to the feature values of the patient.
For numeric features, (e.g. SBS mutational load) the average refers to
the interquartile mean. For boolean features, (e.g. presence of a gene
fusion) the average refers to the proportion of patients with the
feature, with the patient feature values of 0% and 100% indicating
absence/presence respectively.

``` r
patient_report <- patientReport(
   probs=pred_report$prob_scaled, 
   feat.contrib=pred_report$feat_contrib, 
   sample.name='HMF000185A',
   
   ## Shows all 3 panels. Set to 2 to show only the first 2, or 1 to show only the first
   level.of.detail=3 
)
```

``` r
plot(patient_report)
```

![](doc/supp/patient_report_example.png) In the first panel we can see
that this sample was predicted as a liposarcoma (this sample was indeed
labelled as a liposarcoma).

In the middle panel, we can see that the prediction is primarily
supported by this sample having a FUS-DDIT3 fusion (fusion.FUS_DDIT3),
as well as having the RMD profile of liposarcoma (rmd.Sarcoma_Lipo.2).

The third panel shows the feature values compared to the average feature
values in liposarcoma (pink) and in other cancer types (red). We can see
that the liposarcoma prediction for this sample was supported by this
sample having a higher proportion of SBSs contributing (33.5% of SBSs)
to the RMD profile of liposarcoma compared to in other cancer types
(0%). We can also see that FUS-DDIT3 fusions occur in 18.4% of
liposarcomas, but in 0% of other cancer types. The sample value is 1
(100% occurrence) indicating that this feature is present in this
sample.

We can also see in the first panel that Sarcoma_Other also has a high
probability. In cases where the probabilities are more uncertain (like
as shown here), more feature contribution panels corresponding to the
top predicted classes will be shown. This should aid in determining the
cancer type, since it’s possible that e.g. the 2nd predicted cancer type
better matches with evidence from pathology.
