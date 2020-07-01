rule fastqc:
    input:
        fastq=None
    output:
        directory=None
    theads: 8
    shell:
        "fastqc -t {threads} -o {output.directory} {input.fastq}"

rule multiqc:
    input:
        directory=None
    output:
        directory=None
    params:
        pass
    shell:
        "multiqc -o {output.directory} {input.directory}"
