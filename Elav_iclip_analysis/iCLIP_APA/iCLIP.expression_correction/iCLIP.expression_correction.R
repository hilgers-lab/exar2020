library(optparse)

option_list <- list( 
  
  make_option(c("-c", "--iclip.tpm_table"),type="character",
              help="iCLIP xl-count tpm table"),
  make_option(c("-x", "--expression.tpm_table"), type='character',
              help="tpm corrected expression table (tsv)"),
  make_option("--gtf", type="character",
              help="genome annotation file (gtf)"),
  make_option(c('-o',"--outprefix"), type="character",
              help="outprefix (.tsv)"),
  
  make_option('--expression.cutoff', default = 1, type = 'double',
              help = "expression cutoff (%default)"),
  make_option('--verbose', default = FALSE, action="store_true",
              help = "expression cutoff (%default)")
  
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

verbose = opt$verbose

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

# Goal is to correct iCLIP TPMs by gene tpm-corrected expression, gene-level

# TPM corrected iCLIP counts
cat(">>> Importing data\n")

f0 = opt$iclip.tpm_table
names(f0) <- gsub('(.*)\\.tsv$','\\1',basename(f0))
iclip.tab <- read.table(f0,header=TRUE)

dm6.annotation = opt$gtf
dm6.genes = import.gff(dm6.annotation, feature.type = 'gene')
gid2biotype = deframe(mcols(dm6.genes)[c('gene_id','gene_biotype')])
gid2width = deframe(data.frame(mcols(dm6.genes)$gene_id, width(dm6.genes), stringsAsFactors = FALSE))


# expression table, TPM corrected
expr.tab = read_tsv(opt$expression.tpm_table)
expr.tab = expr.tab %>% filter(gid2biotype[gene_id] == 'protein_coding') 

cat(">>> Expression correction\n")

iclip.tab = iclip.tab %>% 
  group_by(group_id) %>% 
  summarize(tpm_iclip = mean(tpm))

# Set of tables, populated with iCount peaks --scores output at different cutoffs
tab0 = iclip.tab %>% 
  left_join(expr.tab, by = c('group_id'= 'gene_id')) %>% 
  mutate(width = gid2width[group_id]) %>% 
  filter(tpm_expression > opt$expression.cutoff)

## Correct for gene expression through linear regression

lm0 = lm(formula = log2(tpm_iclip) ~ log2(tpm_expression), data = tab0)
if(verbose)
  summary(lm0)
iclip.expected = predict(lm0, tab0['tpm_expression'])

iclip_quantification = data.frame(gene_id = tab0$group_id, 
                                       width.gene = tab0$width,
                                       tpm_expr = tab0$tpm_expression,
                                       tpm_iclip = tab0$tpm_iclip,
                                       iclip.expression_fit = iclip.expected) %>% 
  mutate(iclip.expression_ratio = tpm_iclip/tpm_expr,
         iclip.expression_regression = tpm_iclip / 2**iclip_fit)

cat(">>> Write output\n")
# output all tables
outnameSet <- paste0(opt$outprefix,'.expression_correction.tsv')
write_tsv(x = iclip_quantification, path = outnameSet)
