options(stringsAsFactors=F)
devtools::load_all('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/geneDriverAnnotator/')

## Load genes list ----------------------------------------------------------------
genes <- read.delim('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/geneDriverAnnotator/scripts/genes_lists/genes_list.txt')
genes <- genes[!duplicated(genes$gene_name),'gene_name']
genes <- sort(unique(genes))

df <- data.frame(
   gene_name=genes, 
   ensembl_gene_id=geneNamesToEnsg(genes)
)

write.tsv(
   df,
   '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/geneDriverAnnotator/scripts/genes_lists/genes_list_ann.txt'
)
df <- read.delim('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/geneDriverAnnotator/scripts/genes_lists/genes_list_ann.txt')

## Query Biomart ----------------------------------------------------------------
mart <- biomaRt::useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl', host='grch37.ensembl.org')
#biomaRt::listAttributes(mart)
genes_bed <- biomaRt::getBM(
   mart=mart, verbose=F,
   attributes=c(
      'chromosome_name', 'start_position', 'end_position',
      'hgnc_id','hgnc_symbol','ensembl_gene_id'
   ),
   filters='ensembl_gene_id',
   values=df$ensembl_gene_id
)

## Clean up table
genes_bed <- genes_bed[
   genes_bed$ensembl_gene_id %in% df$ensembl_gene_id & ## Only keep genes with ENSG
   !apply(genes_bed,1,anyNA) & ## Rm rows with NA
   !grepl('PATCH',genes_bed$chromosome_name) ## Rm genome patch genes
,]

genes_bed <- genes_bed[order(genes_bed$hgnc_symbol),]

## Manually add missing ensg ids and coords ----------------------------------------------------------------
missing_genes <- df$ensembl_gene_id[ !(df$ensembl_gene_id %in% genes_bed$ensembl_gene_id) ]
df[df$ensembl_gene_id %in% missing_genes,]

write.tsv(
   genes_bed,
   '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/geneDriverAnnotator/scripts/gene_list/genes_pre.txt'
)

## Verify that ENSG and HGNC ids can be converted back and forth
# anyNA(ensgToHgncSymbol(genes_bed$ensembl_gene_id))
# anyNA(geneNamesToEnsg(genes_bed$hgnc_symbol))

## Exons ----------------------------------------------------------------
genes_bed <- read.delim('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/geneDriverAnnotator/scripts/gene_list/genes.txt')

exons_bed_pre <- biomaRt::getBM(
   mart=mart, verbose=F,
   attributes=c(
      'chromosome_name', 'start_position', 'end_position',
      'hgnc_id','hgnc_symbol','ensembl_gene_id','ensembl_exon_id'
   ),
   filters='ensembl_gene_id',
   values=genes_bed$ensembl_gene_id
)

exons <- biomaRt::getBM(
   mart=mart, verbose=F,
   attributes=c(
      'ensembl_exon_id','exon_chrom_start','exon_chrom_end'
   ),
   filters='ensembl_exon_id',
   values=exons_bed_pre$ensembl_exon_id
)
colnames(exons) <- c('ensembl_exon_id','exon_start','exon_end')

exons_bed <- merge(exons_bed_pre, exons, by='ensembl_exon_id', all=T)
sel_cols <- c(
   chrom='chromosome_name',
   exon_start='exon_start',
   exon_end='exon_end',
   ensembl_exon_id='ensembl_exon_id',
   hgnc_symbol='hgnc_symbol',
   ensembl_gene_id='ensembl_gene_id',
   gene_start='start_position',
   gene_end='end_position'
)
exons_bed$hgnc_symbol <- as.character(genes_bed$hgnc_symbol[ match(exons_bed$ensembl_gene_id, genes_bed$ensembl_gene_id) ])
exons_bed <- exons_bed[order(exons_bed$hgnc_symbol),]
exons_bed <- exons_bed[,sel_cols]
colnames(exons_bed) <- names(sel_cols)
exons_bed <- exons_bed[exons_bed$ensembl_gene_id %in% genes_bed$ensembl_gene_id,]

#table(exons_bed$chrom)

## Fix missing, and PATCH ----------------------------------------------------------------
genes_bed$ensembl_gene_id[ !(genes_bed$ensembl_gene_id %in% exons_bed$ensembl_gene_id) ]
write.tsv(
   exons_bed,
   '/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/geneDriverAnnotator/scripts/gene_list/exons_pre.txt'
)
