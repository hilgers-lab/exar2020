library(optparse)

option_list <- list( 
  make_option(c("-p", "--peakset"),type="character",
              help="iCount peak scores (tsv), comma-separated list"),
  make_option("--gtf", type="character",
              help="genome annotation file (gtf)"),
  make_option(c("-o","--outprefix"), type="character",
              help="outprefix (.tsv)"),  
  make_option("--fdr.cutoff", default = 0.05, type = 'double', 
               help = "Default chromosomes (%default)"),
  make_option("--chroms", default = '2L,2R,3L,3R,4,X', type = 'character', 
               help = "Default chromosomes (%default)")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))


suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))

# michael's version
# https://support.bioconductor.org/p/91218/
#
# vector version. matrix wont work
tpm3 <- function(counts.vector,len) {
  x <- counts.vector/len
  return((x*1e6)/sum(x))
}

params=list()
params$chroms.valid = strsplit(opt$chroms,',')[[1]]
params$fdr.cutoff = opt$fdr.cutoff

dm6.transcripts = import.gff(opt$gtf, feature.type = 'transcript')
tx2gid = deframe(mcols(dm6.transcripts)[c('transcript_id','gene_id')])

dm6.genes = import.gff(opt$gtf, feature.type = 'gene')
gid2width = deframe(data.frame(mcols(dm6.genes)$gene_id, width(dm6.genes)))
gid2gene_biotype = deframe(mcols(dm6.genes)[c('gene_id','gene_biotype')])


txdb0 = makeTxDbFromGFF(opt$gtf)
dm6.utr3 = threeUTRsByTranscript(txdb0, use.names = TRUE)
utr2width = dm6.utr3 %>% as.data.frame %>%
  mutate(gene_id = tx2gid[group_name]) %>%
  group_by(gene_id) %>%
  summarize(width = sum(width)) %>% deframe

f0 = strsplit(opt$peakset,',')[[1]]
names(f0) = gsub('.tsv$','',basename(f0))

coltypes0 = cols(
  chrom = col_character(),
  position = col_double(),
  strand = col_character(),
  name = col_character(),
  group_id = col_character(),
  score = col_double(),
  score_extended = col_character(),
  FDR = col_double()
)
xltabSet = lapply(f0, read_tsv, col_names = TRUE, col_types = coltypes0)

xltabSet_filtered = lapply(xltabSet, subset, 
                           FDR < params$fdr.cutoff & 
                             chrom %in% params$chroms.valid)
  
xlsitesSet = lapply(xltabSet_filtered, function(tab) {
  makeGRangesFromDataFrame(tab, keep.extra.columns = TRUE,
                           seqnames.field = 'chrom',
                           start.field = 'position', 
                           end.field = 'position')})
  
# discard ambiguous xl-sites
coord = lapply(xlsitesSet, function(gr) {paste0(seqnames(gr),'_',start(gr))})
coord.uniq = lapply(coord, function(x) names(which(table(x) == 1)))
  
coord.select = mapply(function(x, y) x %in% y, coord, coord.uniq)
xlsitesSet.uniq = mapply(function(gr, b) gr[b], xlsitesSet, coord.select, SIMPLIFY = FALSE) %>% 
  GRangesList()
  
tab0 = xlsitesSet.uniq %>% as.data.frame

targets.stats = tab0 %>% group_by(group_id, group_name) %>% 
  summarize(xlsite_count = n(),
            score_sum = sum(score),
            padj.min = min(FDR)) %>% 
  mutate(width = gid2width[group_id],
         gene_biotype = gid2gene_biotype[group_id]) %>% 
  ungroup %>%
  group_by(group_name) %>% 
  mutate(tpm = tpm3(score_sum, width)) %>% 
  ungroup %>% 
  arrange(group_name)

outfname = paste0(opt$outprefix,'.genebody_tpm.fdr',params$fdr.cutoff,'.tsv')
print(outfname)
write_tsv(targets.stats, outfname)
