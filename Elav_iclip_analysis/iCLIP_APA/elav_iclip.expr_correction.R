library(readr)
library(rtracklayer)
library(tibble)
library(dplyr)
library(reshape2)
library(ggplot2)

# Goal is to correct iCLIP TPMs by gene expression

# TPM corrected iCLIP counts
f0 = c("<set of tsv>")
names(f0) <- gsub('(.*)\\.tsv$','\\1',basename(f0))
iclip.tabSet <- lapply(f0, read.table,header=TRUE)

dm6.annotation = "<dm6_ensembl>.gtf"
dm6.genes = import.gff(dm6.annotation, feature.type = 'gene')
gid2biotype = deframe(mcols(dm6.genes)[c('gene_id','gene_biotype')])
gid2width = deframe(data.frame(mcols(dm6.genes)$gene_id, width(dm6.genes), stringsAsFactors = FALSE))


# expression table, TPM corrected
expr.tab = read_tsv("<tpm values of gene expression>")
expr.tab = expr.tab %>% filter(gid2biotype[gene_id] == 'protein_coding') 
colnames(expr.tab)[6] = 'tpm_expression'


iclip.tabSet <- lapply(iclip.tabSet, function(tab) { tab %>% group_by(group_id) %>% 
  summarize(tpm_iclip = mean(tpm))})

# Set of tables, populated with iCount peaks --scores output at different cutoffs
tab0.set <- lapply(iclip.tabSet, function(tab, tab.expr) {
  tab %>% left_join(tab.expr, by = c('group_id'= 'gene_id'))
}, expr.tab[c('gene_id','tpm_expression')])

tab0.set <- lapply(tab0.set, function(tab, g2w){
  tab %>% mutate(width = g2w[group_id])
}, gid2width)

tab1.set = lapply(tab0.set, function(tab) tab %>% filter(tpm_expression > 1))
## Correct for gene expression through linear regression

# per iCLIP peaks FDR cutoff
lm0.set <- lapply(tab1.set, function(tab){
  lm(formula = log2(tpm_iclip) ~ log2(tpm_expression), 
     data = tab)
})

# lapply(lm0.set, summary)

iclip.expectSet <- mapply(predict, lm0.set, lapply(tab0.set,'[','tpm_expression'))


elav.iclip_quantificationSet <- mapply(function(tab, lm0) {
  iclip.expect <- predict(lm0, tab['tpm_expression'])
  
  elav.iclip_quantification = data.frame(gene_id = tab$group_id, 
                                         width.gene = tab$width,
                                         tpm_expr = tab$tpm_expression,
                                         tpm_iclip = tab$tpm_iclip,
                                         iclip_fit = iclip.expect) %>% 
    mutate(iclip.regression = tpm_iclip / 2**iclip_fit,
           iclip.ratio = tpm_iclip/tpm_expr)
  return(elav.iclip_quantification)
},tab0.set, lm0.set, SIMPLIFY = FALSE)

  
# output all tables
outnameSet <- paste0('./elav_iclip.expr_corrected/elav_iclip.expression_corrected.peak_',
                     gsub('elav_iclip.*\\.(fdr.*)$','\\1',
                          names(elav.iclip_quantificationSet)),'.tsv')

void <- mapply(function(tab, outfile){
    write_tsv(x = tab, path = outfile)
  }, elav.iclip_quantificationSet, outnameSet, SIMPLIFY = FALSE)
elementNROWS(void)

