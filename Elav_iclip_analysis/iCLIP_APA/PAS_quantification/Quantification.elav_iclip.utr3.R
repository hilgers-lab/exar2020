suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(tibble))

tpm4 <- function(counts, counts.total, len) {
  x <- counts/len
  return((x*1e6)/counts.total)
}

dm6.genes = import.gff('../dm6_ensembl96.gtf', feature.type = 'gene')
gid2symbol = mcols(dm6.genes)[c('gene_id','gene_name')] %>% deframe

workdir="" # snakemake workdir

ctsMatSet_files <- list.files(paste0('../',workdir,'/coverage/'),pattern = '\\.bed$', full.names = TRUE)
names(ctsMatSet_files) = gsub('elav_iclip\\.(.*_cDNA)-elav_wt.*sorted\\.(.*)\\.bed','\\1.\\2', basename(ctsMatSet_files))
ctsMatSet = lapply(ctsMatSet_files, read_tsv)

zero.anno <- lapply(ctsMatSet, function(x) x %>% mutate(is.zero = as.logical(Counts == 0)) %>% as.data.frame)

# NOTE: report features with zero counts before adding the pseudocount
for(n in names(zero.anno)){
  fname0 = paste0('./Zero.annotation/',n,'.counts_replicate.tsv')
  write_tsv(as.data.frame(zero.anno[[n]]), fname0)
}

ctsMatSet = lapply(ctsMatSet, function(x) {x$width = x$End - x$Start; x})
#  add psuedocount
ctsMatSet <- lapply(ctsMatSet, function(x) {x$Counts = x$Counts + 1; x})

libsizeSet_files = list.files(paste0('../',workdir,'/libsize/'), pattern = 'libsize$', full.names = TRUE)
names(libsizeSet_files) = gsub('.*(NNN.*NN).*\\.libsize$','\\1', basename(libsizeSet_files))
libsizes <- sapply(lapply(libsizeSet_files, read_tsv, col_names = FALSE),'[[',1)


# precomputed expression fit
expression.fit = "<>" # require columns gene_id, iclip_fit

ge.fit = read_tsv(expression.fit)
ge.fit.map = deframe(ge.fit[c('gene_id','iclip_fit')])

## TODO: Update name of replicates and expand, if necessary
gene_ids0 = ctsMatSet$`elav_iclip.NNNGGTTNN_cDNA-Redundancy.elav_wt.early.sorted.breakpoints_flank100.bed`$Gene
width0 = ctsMatSet$`elav_iclip.NNNGGTTNN_cDNA-Redundancy.elav_wt.early.sorted.breakpoints_flank100.bed`$width
cts.dd_rep <- data.frame(gene_id = gene_ids0, 
                         gene_name = gid2symbol[gene_ids0],
                         width = width0,
                         breakpts_<rep> = ctsMatSet[[replicate-breakpts]]$Counts,
                         tes_<rep> = ctsMatSet[[replicate-tes]]$Counts,
                         stringsAsFactors = FALSE)

# TODO: rename output file for count table for replicates
write_tsv(cts.dd_rep, '<dataset>.counts_replicates.tsv')

# tpm normalization
# TODO: update sample name
tpm.dd_rep = cts.dd_rep
for(n in setdiff(colnames(cts.dd_rep), c('gene_id','gene_name', 'width'))){
  bc0 = toupper(gsub('.*_(.*)$','NNN\\1NN',n))
  tpm.dd_rep[,n] = tpm4(cts.dd_rep[,n], libsizes[bc0] + nrow(cts.dd_rep), cts.dd_rep$width)
}
# TODO: rename output file for tpm table for replicates
write_tsv(tpm.dd_rep, '<dataset>.tpm_replicates.tsv')


# pool replicates
tpm.dd <- data.frame(gene_id = tpm.dd_rep$gene_id, 
                     gene_name = tpm.dd_rep$gene_name, 
                     width = tpm.dd_rep$width, 
                     breakpts = rowMeans(tpm.dd_rep[,grepl('^breakpts',colnames(tpm.dd_rep))]),
                     tes = rowMeans(tpm.dd_rep[,grepl('^tes',colnames(tpm.dd_rep))]),stringsAsFactors = FALSE)
tpm.dd %>% arrange(desc(breakpts)) %>% head(n = 20)
# TODO: update name for average iCLIP TPM 
write_tsv(tpm.dd, './<dataset>.tpm.tsv')
# ge correction
rna.ge <- ge.fit[match(tpm.dd$gene_id, ge.fit$gene_id),]
tpm_ge.dd <- tpm.dd %>% mutate(breakpts = breakpts / 2**rna.ge$iclip_fit,
                               tes = tes / 2**rna.ge$iclip_fit)
# TODO: update name for average iCLIP TPM, gene expression corrected 
write_tsv(tpm_ge.dd, './<dataset>.tpm_GEcorrected.tsv')

