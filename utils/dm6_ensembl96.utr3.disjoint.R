library(GenomicFeatures)
library(rtracklayer)
library(tibble)

dm6_tx <- import.gff('./dm6_ensembl96.gtf', feature.type = 'transcript')
tx2gid <- deframe(mcols(dm6_tx)[,c('transcript_id','gene_id')])
txdb1 <- makeTxDbFromGFF('./dm6_ensembl96.gtf')
dm6.utr3 <- threeUTRsByTranscript(txdb1, use.names = TRUE)

utr3.segments <- disjoin(dm6.utr3)
# add segment order
df <- as.data.frame(utr3.segments)
df.grouped <- split(df, f = df$group_name)
df.grouped_anno <- lapply(df.grouped, function(x) {
  if(unique(x$strand == '+'))
    x$order = order(x$start)
  else if(unique(x$strand == '-'))
    x$order = rev(order(x$start))
  else
    x$order = NA
  return(x)
})
df.grouped_anno <- do.call(rbind,df.grouped_anno)
gr1 <- makeGRangesFromDataFrame(df.grouped_anno, keep.extra.columns = TRUE)
names(gr1) = NULL
mcols(gr1)$gene_id <- tx2gid[mcols(gr1)$group_name]

export(object = unique(gr1), './dm6_ensembl96.utr3.disjoint.gtf')
