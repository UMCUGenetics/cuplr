CUPLR: Cancer of Unknown Primary Location Resolver
================

CUPLR is a suite of R packages for detecting the primary tumor location
of metastatic tumors based on features derived from whole genome
sequencing data.

This includes the following R packages:

  - `cuplr`: random forest training and prediction
  - `featureExtractor`: extract all features used by CUPLR
  - `geneDriverAnnotator`: annotate gene amplifications and biallelic
    loss
  - `commonUtils`: code used by all CUPLR R packages

The features used by CUPLR fall under the following types:

  - `sigs`: SNV signatures, indel contexts, MNV counts
  - `kataegis`: counts of kataegis foci
  - `gene_amp`: focal gene amplifications
  - `gene_def`: biallelic gene loss
  - `aneuploidy`: chromosome arm gains/losses compared to the overall
    genome ploidy
  - `purple`: output from PURPLE; includes gender and whole genome
    duplication status
  - `sv_types`: \[LINX\] simple and complex SV events
  - `viral_ins`: \[LINX\] presence of viral insertions
  - `fusion`: \[LINX\] presence of gene fusions
  - `rep_elem`: counts of repetitive elements
  - `rmd`: regional mutational density

## Dependencies

CUPLR depends on several R packages. These can be installed as follows:

``` r
## cuplr dependencies
install.packages('rfFC')
install.packages('randomForest')

## featureExtractor dependencies
install.packages("seqminer")

## commonUtils dependencies
install.packages("BiocManager")
BiocManager::install("GenomeInfoDb")

## Devtools required to load CUPLR R packages directly
install.packages("devtools")
```

`geneDriverAnnotator` also requires SnpEff/SnpSift and java. However,
these have been embedded at `geneDriverAnnotator/dep/` and do not need
to be installed separately

In bash, run the following commands to download CUPLR and
mutSigExtractor:

    cd /working/dir/
    git clone https://github.com/UMCUGenetics/cuplr/
    git clone https://github.com/UMCUGenetics/mutSigExtractor/

## Loading CUPLR

In R, load CUPLR and its dependencies:

``` r
setwd('/working/dir/')

## These need to be installed (done in previous step) and loaded with library()
## This is because these packages contain C code that must be first compiled
library('randomForest')
library('rfFC')
library("seqminer")

## Load these packages directly with devtools
devtools::load_all('mutSigExtractor/') 

devtools::load_all('cuplr/commonUtils/')
## Load commonUtils first! The other CUPLR R packages rely on code in commonUtils

devtools::load_all('cuplr/cuplr/')
devtools::load_all('cuplr/featureExtractor/')
devtools::load_all('cuplr/geneDriverAnnotator/')
```

## Extract features

Extract the features from the HMF pipeline outputs. Note that the gene
driver annotation takes some time (\~2-5min). The other features take in
general a few seconds each. However, extracting signatures, kataegis and
RMD can take some time (\>5min each) in the case of hypermutators
(i.e. large somatic SNV/indel vcf).

``` r
extractFeaturesCuplr(
   ## vcfs
   germ.vcf.path='/path/to/*.annotated.vcf.gz',
   som.vcf.path='/path/to/*.purple.somatic.vcf.gz',
   sv.vcf.path='/path/to/*.purple.sv.ann.vcf.gz',

   ## PURPLE output
   purple.cnv.path='/path/to/*.purple.cnv.somatic.tsv',
   purple.purity.path='/path/to/*.purple.purity.tsv',

   ## LINX output
   linx.fusion.path='/path/to/*.linx.fusion.tsv',
   linx.viral.inserts.path='/path/to/*.linx.viral_inserts.tsv',
   linx.vis.sv.data.path='/path/to/*.linx.vis_sv_data.tsv',

   out.dir='/path/to/output/dir/', 
   sample.name='sample_name', 
   verbose=2
   )
```

This will generate `features.txt.gz` within `out.dir`, a 1-row dataframe
with feature name header.

Use `readFeaturesCuplr()` to read `features.txt.gz` into R. This
function will load the txt file as a dataframe and assign factor levels
to categorical features (e.g. `gene_def`). Below is what the first few
columns look like. Feature names in the header follow the regex format
`^feature_type.feature_name`.

``` r
devtools::load_all('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/commonUtils/')
devtools::load_all('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/featureExtractor/')
devtools::load_all('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/cuplr/')
features <- readFeaturesCuplr('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/featureExtractor/test/CPCT02020731T/features.txt.gz')
features[,1:5]
```

    ##   sigs.snv.SBS1 sigs.snv.SBS2 sigs.snv.SBS3 sigs.snv.SBS4 sigs.snv.SBS5
    ## 1      613.0238      135.3868             0      70.66962      813.3634

## Predict cancer type

First load CUPLR with `readRDS()`. Then use
`predict.randomForestEnsemble()` to obtain the probabilities for each
cancer type.

``` r
cuplr <- readRDS(CUPLR) ## `CUPLR` refers to the default path to the random forest model
predict.randomForestEnsemble(cuplr, features, type='prob')
```

    ##   Biliary Bone_SoftTissue Breast Colon_Rectum Gastric Head_and_neck Kidney
    ## 1       0               0      0            0       0             0  0.002
    ##   Liver Lung Lymphoid Mesothelioma Nervous_system Neuroendocrine Ovary Pancreas
    ## 1     0    0        0            0              0              0     0        0
    ##   Prostate Skin Urinary_tract Uterus
    ## 1    0.996    0         0.002      0
    ## attr(,"class")
    ## [1] "matrix" "votes"
