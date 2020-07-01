rule segment:
    input:
        gtf=config['annotation'],
        fai=config['fasta_index']
    output:
        gtf="segment/dm6_icount.segment.gtf",
        regions="segment/regions.gtf.gz"
    params:
        logging="-P segment/dm6_icount.segment.log --file_log 1 --results_file segment/dm6_icount.segment.metric"
    shell:
        "iCount segment {input.gtf} {output.gtf} {input.fai} {params}"

# demux needs be done externally
rule mapstar:
    input:
        fastq=os.path.join(config["fastq_directory"], "{sample}.fastq.gz"),
        star_index=config['star_index']
    output:
        outbam='mapstar/{sample}/Aligned.sortedByCoord.out.bam'
    params:
        outdir='mapstar/{sample}',
        logging="--file_logpath mapstar/{sample}.log --file_log 1 --results_file mapstar/{sample}.metrics",
        annotation=config['annotation']
    threads: 8
    shell:
        "iCount mapstar {input.fastq} {input.star_index} {params.outdir} --annotation {params.annotation} --threads {threads} {params.logging}"

rule xlsites:
    input:
        bam=rules.mapstar.output.outbam,
        segments=rules.segment.output.gtf
    output:
        singles=("xlsites/{sample}_{quantification_mode}_singles.bed"),
        multi=("xlsites/{sample}_{quantification_mode}_multi.bed"),
        skipped=("xlsites/{sample}_{quantification_mode}_skipped.bam")
    params:
        quantification="{quantification_mode}",
        mapq_thrs=3,
        logging="--report_progress --file_logpath xlsites/{sample}_{quantification_mode}.log --file_log 1 --results_file xlsites/{sample}_{quantification_mode}.metric",
        segmentation="--segmentation {input.segments}" if config['use_segmentation'] else ""
    benchmark:
        "xlsites/benchmark/{sample}_{quantification_mode}_singles.benchmark.txt"
    shell:
        "iCount xlsites --group_by start --quant {params.quantification} --mapq_th {params.mapq_thrs} \
        {params.segmentation} {params.logging} \
        {input.bam} {output.singles} {output.multi} {output.skipped}"

rule link_bams:
    input:
        bam=rules.mapstar.output.outbam
    output:
        bam=os.path.join("Alignments","{sample}.bam"),
        bai=os.path.join("Alignments","{sample}.bam.bai")
    threads: 4
    shell:
        "ln -s $(realpath {input.bam}) {output.bam} && samtools index -@ {threads} {output.bam}"

rule peaks:
    input:
        annotation=rules.segment.output.gtf,
        xlsites=rules.xlsites.output.singles
    output:
        peaks=os.path.join("Peak_calling","peaks","{sample}_{quantification_mode}_peaks.bed"),
    params:
        half_window="--half_window 3",
        fdr="--fdr " + str(config_peakcalling['fdr']),
        perm="--perms " + str(config_peakcalling['nperm']),
        scores="--scores Peak_calling/peaks/{sample}_{quantification_mode}_peaks.tsv",
        logging="--file_logpath Peak_calling/peaks/{sample}_{quantification_mode}_peaks.log --file_log 1 --results_file Peak_calling/peaks/{sample}_{quantification_mode}_peaks.metrics",
    shell:
        "iCount peaks {params} {input.annotation} {input.xlsites} {output.peaks}"

rule peaks_transcript:
    input:
        annotation=rules.segment.output.gtf,
        xlsites=rules.xlsites.output.singles
    output:
        peaks=os.path.join("Peak_calling","peaks_{feature}","{sample}_{quantification_mode}_peaks.bed"),
    params:
        feature = '--features {feature} -g transcript_id',
        half_window="--half_window 3",
        fdr="--fdr " + str(config_peakcalling['fdr']),
        perm="--perms " + str(config_peakcalling['nperm']),
        scores="--scores Peak_calling/peaks_{feature}/{sample}_{quantification_mode}_peaks.tsv",
        logging="--file_logpath Peak_calling/peaks_{feature}/{sample}_{quantification_mode}_peaks.log \
        --file_log 1 --results_file Peak_calling/peaks_{feature}/{sample}_{quantification_mode}_peaks.metrics",
    shell:
        "iCount peaks {params} {input.annotation} {input.xlsites} {output.peaks}"


rule clusters:
    input:
        sites=rules.xlsites.output.singles,
        peaks=rules.peaks.output.peaks,
    output:
        clusters=os.path.join("Peak_calling", "clusters", "{sample}_{quantification_mode}_clusters.bed")
    params:
        dist="--dist 20",
        slop="--slop 3",
        logging="--file_logpath Peak_calling/clusters/{sample}_{quantification_mode}_clusters.log \
            --file_log 1 --results_file Peak_calling/clusters/{sample}_{quantification_mode}_clusters.metric"
    shell:
        "iCount clusters {params} {input.sites} {input.peaks} {output.clusters}"


rule clusters_transcript:
    input:
        sites=rules.xlsites.output.singles,
        peaks=rules.peaks_transcript.output.peaks,
    output:
        clusters=os.path.join("Peak_calling", "clusters_{feature}", "{sample}_{quantification_mode}_clusters.bed")
    params:
        dist="--dist 20",
        slop="--slop 3",
        logging="--file_logpath Peak_calling/clusters_{feature}/{sample}_{quantification_mode}_clusters.log \
            --file_log 1 --results_file Peak_calling/clusters_{feature}/{sample}_{quantification_mode}_clusters.metric"
    shell:
        "iCount clusters {params} {input.sites} {input.peaks} {output.clusters}"

rule annotate_xlsites:
    input:
        annotation=rules.segment.output.gtf,
        xlsites=rules.xlsites.output.singles
    output:
        annotated=os.path.join("annotate_xlsites", "{sample}_{quantification_mode}_singles_{annotation_type}.bed")
    params:
        subtype="--subtype {annotation_type}",
        logging="--file_logpath annotate_xlsites/{sample}_{quantification_mode}_{annotation_type}.log \
        --file_log 1 --results_file annotate_xlsites/{sample}_{quantification_mode}_{annotation_type}.metrics"
    shell:
        "iCount annotate {params} {input.annotation} {input.xlsites} {output.annotated}"

rule annotate_peaks:
    input:
        annotation=rules.segment.output.gtf,
        peaks=rules.peaks.output.peaks
    output:
        annotated=os.path.join("annotate_peaks", "{sample}_{quantification_mode}_peaks_{annotation_type}.bed")
    params:
        subtype="--subtype {annotation_type}",
        logging="--file_logpath annotate_peaks/{sample}_{quantification_mode}_peaks_{annotation_type}.log \
        --file_log 1 --results_file annotate_peaks/{sample}_{quantification_mode}_peaks_{annotation_type}.metrics"
    shell:
        "iCount annotate {params} {input.annotation} {input.peaks} {output.annotated}"


rule summary_xlsites:
    input:
        annotation=rules.segment.output.regions,
        xlsites=rules.xlsites.output.singles
    output:
        genesummary=os.path.join("summary_xlsites","summary_{sample}_{quantification_mode}","summary_gene.tsv")
    params:
        directory="summary_xlsites/summary_{sample}_{quantification_mode}",
        logging="--file_logpath summary_xlsites/summary_{sample}_{quantification_mode}.log --file_log 1 --results_file summary_xlsites/summary_{sample}_{quantification_mode}.metrics"
    shell:
        "iCount summary {params.logging} {input.annotation} {input.xlsites} {params.directory}"


rule summary_peaks:
    input:
        annotation=rules.segment.output.regions,
        peaks=rules.peaks.output.peaks
    output:
        genesummary=os.path.join("summary_peaks","summary_{sample}_{quantification_mode}","summary_gene.tsv")
    params:
        directory=os.path.join("summary_peaks", "summary_{sample}_{quantification_mode}"),
        logging=' '.join(["--file_logpath ", os.path.join(" summary_peaks","summary_{sample}_{quantification_mode}.log"), \
        "--file_log 1", "--results_file", os.path.join("summary_peaks", "summary_{sample}_{quantification_mode}.metrics")])
    shell:
        "iCount summary {params.logging} {input.annotation} {input.peaks} {params.directory}"


## do quantification on annotation_type = biotype. Select right input, and output
rule summary_peaks_per_feature:
    input:
        annotation=rules.segment.output.regions,
        peaks=rules.annotate_peaks.output.annotated
    output:
        bedfiltered=os.path.join('summary_per_feature_peaks','summary_{sample}_{quantification_mode}_peaks_{annotation_type}.{feature}.bed'),
        genesummary=os.path.join("summary_per_feature_peaks","summary_{sample}_{quantification_mode}_peaks_{annotation_type}.{feature}","summary_gene.tsv")
    params:
        feature="\"{feature}\"",
        directory=os.path.join("summary_per_feature_peaks", "summary_{sample}_{quantification_mode}_peaks_{annotation_type}.{feature}"),
        logging="--file_logpath summary_per_feature_peaks/summary_{sample}_{quantification_mode}_peaks_{annotation_type}.{feature}.log " \
        "--file_log 1 --results_file summary_per_feature_peaks/summary_{sample}_{quantification_mode}_peaks_{annotation_type}.{feature}.metrics"
    shell:
        "grep {params.feature} {input.peaks} > {output.bedfiltered} && "
        "iCount summary {params.logging} {input.annotation} {output.bedfiltered} {params.directory}"


# signal tracks, cannot be converted to bigwig for bedgraph is strand specific
rule bedgraph:
    input:
        xlsites=rules.peaks.output.peaks
    output:
        bedgraph='bedgraph/{sample}_{quantification_mode}_peaks.bedgraph'
    params:
        name="--name {sample}_{quantification_mode}",
        logging="--file_logpath bedgraph/{sample}_{quantification_mode}_peaks.log \
        --file_log 1 --results_file bedgraph/{sample}_{quantification_mode}_peaks.metrics"
    shell:
        "iCount bedgraph {params} {input.xlsites} {output.bedgraph}"

rule bedgraph_transcript:
    input:
        xlsites=rules.peaks_transcript.output.peaks
    output:
        bedgraph='bedgraph_{feature}/{sample}_{quantification_mode}_peaks.bedgraph'
    params:
        name="--name {sample}_{quantification_mode}",
        logging="--file_logpath bedgraph_{feature}/{sample}_{quantification_mode}_peaks.log \
        --file_log 1 --results_file bedgraph_{feature}/{sample}_{quantification_mode}_peaks.metrics"
    shell:
        "iCount bedgraph {params} {input.xlsites} {output.bedgraph}"
