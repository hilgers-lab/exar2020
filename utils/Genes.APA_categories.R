suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(UpSetR))

# this script takes up the genome annotation 
# and identifies the APA type of a gene
# APA types may coincide
#
# output: is a table per gene, describing the 
# APA status {TRUE, FALSE} per category

# parameter
params = list()
params$chroms.valid = c('2L','2R','3L','3R','4','X')

# 
txdb0 = makeTxDbFromGFF('/data/hilgers/group/Paper_ELAV-dependentTargets/Annotations/dm6_ensembl96.gtf')
dm6.transcripts = import.gff('/data/hilgers/group/Paper_ELAV-dependentTargets/Annotations/dm6_ensembl96.gtf',
                             feature.type = 'transcript')
dm6.genes = import.gff('/data/hilgers/group/Paper_ELAV-dependentTargets/Annotations/dm6_ensembl96.gtf',
                       feature.type = 'gene')

# derivatives
tx2gid = deframe(mcols(dm6.transcripts)[c('transcript_id','gene_id')])
geneBiotype.map = deframe(mcols(dm6.genes)[c('gene_id','gene_biotype')])

# stop codon count
dm6.cds = cdsBy(txdb0,by = 'tx', use.names =TRUE) %>% 
  as.data.frame %>% 
  filter(seqnames %in% params$chroms.valid)
dm6.cds_last = dm6.cds %>% 
  group_by(group_name) %>% 
  mutate(is.last = ifelse(unique(strand) == '-',
                          min(cds_id) ,
                          max(cds_id)) == cds_id,
         gene_id = tx2gid[group_name]) %>% 
  filter(is.last) %>% 
  ungroup  %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)

dm6.stop_codon = unique(resize(dm6.cds_last, width = 1, fix = 'end'))

dm6.stop_codon.map = dm6.stop_codon %>% 
  as.data.frame %>% 
  group_by(gene_id) %>% 
  add_tally() %>% 
  dplyr::select(gene_id, n) %>%
  deframe

dm6.utr3_df = threeUTRsByTranscript(txdb0, use.names = TRUE) %>% as.data.frame %>%
  filter(seqnames %in% params$chroms.valid)
dm6.utr3_df$transcript_id = dm6.utr3_df$group_name
dm6.utr3_df$gene_id = tx2gid[dm6.utr3_df$group_name]
dm6.utr3_df$group_name <- NULL

dm6.utr3_byTx = makeGRangesListFromDataFrame(dm6.utr3_df, split.field = 'transcript_id', keep.extra.columns = TRUE)
dm6.utr3_byGene = makeGRangesListFromDataFrame(dm6.utr3_df, split.field = 'gene_id', keep.extra.columns = TRUE)

###################
# Set definitions #
GID.categories = list()

## gene w/o UTR3, i.e. nonCoding variants
# gid_nonCoding = genes(txdb0) %>% as.data.frame %>% 
#   filter(!gene_id %in% names(dm6.utr3_byGene) & seqnames %in% params$chroms.valid) %>% 
#   dplyr::select(gene_id) %>%
#   unlist
# gid_nonCoding = split(gid_nonCoding, geneBiotype.map[gid_nonCoding])
# GID.categories = c(GID.categories, gid_nonCoding)

## singletons - one transcript per gene
gid.singletons <- c(names(which(elementNROWS(unique(split(dm6.transcripts, mcols(dm6.transcripts)$gene_id))) == 1)),
                    names(which(elementNROWS(dm6.utr3_byGene %>% unique) == 1)))

GID.categories$Singletons = gid.singletons

gid.singleTES <- names(which(elementNROWS(unique(resize(dm6.utr3_byGene, fix = 'end', width = 1))) == 1))
GID.categories$has.SingleTES <- gid.singleTES

## APA genes
gid.apa <- setdiff(names(which(elementNROWS(unique(dm6.utr3_byGene)) > 1)),
                   names(which(elementNROWS(unique(dm6.utr3_byTx)) > 1)))
GID.categories$APA_targets <- gid.apa

intersect(gid.singletons, gid.apa)
## spliced UTR3s - multiple UTR3s per transript
# type 1) simple UTR3 - same stop codon, different TES. Variant: spliced UTR3
# type 2) 

utr3_byGene.apa_splicedUTR3 <- tx2gid[names(which(elementNROWS(unique(dm6.utr3_byTx)) > 1))]
GID.categories$SplicedUTR3 = utr3_byGene.apa_splicedUTR3


## Tandem UTR3
dm6.utr3_byTx.uniq <- unique(dm6.utr3_byGene) %>% as.data.frame %>% 
  makeGRangesListFromDataFrame(split.field = 'transcript_id', keep.extra.columns = TRUE)

dm6.utr3_byGene[gid.apa] %>% reduce %>% as.data.frame 
gid.tandem = names(which((dm6.utr3_byGene[gid.apa] %>% reduce %>% elementNROWS()) < 
  (dm6.utr3_byGene[gid.apa] %>% unique %>% elementNROWS())))

gid.tandemUTR3 <- gid.tandem

### exclude range utr3, for AS singletons, i.e. all TX-tes are the same
GID.categories$TandemUTR3 = gid.tandemUTR3

## AL genes - inclusive. Multiple UTR3 per genes with offset
# dm6.cds_DS = resize(makeGRangesFromDataFrame(dm6.cds, keep.extra.columns = TRUE), width = 1, fix = 'end')

# results in 3'utrs overlapping stop codons

GID.categories$AlternativeLast = names(which((dm6.utr3_byGene[gid.apa] %>% reduce %>% elementNROWS()) > 1))

intersect(GID.categories$AlternativeLast,
        c(GID.categories$Singletons, GID.categories$has.SingleTES))


gene_table = as_tibble(mcols(dm6.genes)[,c('gene_id','gene_name', 'gene_biotype')])
gene_table$stop.codon_count = dm6.stop_codon.map[gene_table$gene_id]
gene_table$is.Singleton = gene_table$gene_id %in% GID.categories$Singletons
gene_table$is.TandemUTR3 = gene_table$gene_id %in% setdiff(GID.categories$TandemUTR3, 
                                                           c(GID.categories$Singletons,GID.categories$has.SingleTES))
gene_table$is.AlternativeLast = gene_table$gene_id %in% setdiff(GID.categories$AlternativeLast,
                                                                c(GID.categories$Singletons, GID.categories$has.SingleTES))
gene_table$has.SplicedUTR3 = gene_table$gene_id %in% GID.categories$SplicedUTR3
gene_table$has.singleTES = gene_table$gene_id %in% GID.categories$has.SingleTES



c.cols=colnames(gene_table)[-c(1:4)]
ddlist <- list()
for(i in seq_along(c.cols))
  ddlist[[c.cols[i]]] = gene_table$gene_id[which(gene_table[[c.cols[i]]])]
upset(fromList(ddlist), nsets = 32, nintersects = 32)

# write_tsv(gene_table, './Genes.APA_classification/dm6_ensembl96.APA_status.tsv')
# 
# write_tsv(x = (gene_table %>% dplyr::filter(has.SplicedUTR3 & !(is.AlternativeLast|is.TandemUTR3|is.Singleton)) %>%
#                  as.data.frame)[,1:4],
#           path = './Genes.APA_classification/dm6_ensembl96.APA_status.unmatched_category.tsv')

