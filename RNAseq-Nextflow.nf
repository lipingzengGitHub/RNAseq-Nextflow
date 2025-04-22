nextflow.enable.dsl=2

params {
    reads       = 'RawData/*_{1,2}.fq.gz'
    genome      = 'Ref/genome.fa'
    gtf         = 'Ref/genes.gtf'
    index       = 'Ref/hisat2_index/genome'   // HISAT2 index basename
    outdir      = 'results'
    rscript     = 'scripts/R_script_DEG.R'
}

process fastqc {
    tag "$sample_id"

    container 'biocontainers/fastqc:v0.11.9_cv8'

    input:
    tuple val(sample_id), path(reads)

    output:
    path("${sample_id}/fastqc")

    script:
    """
    mkdir -p ${sample_id}/fastqc
    fastqc ${reads.join(' ')} -o ${sample_id}/fastqc
    """
}

process hisat2_align {
    tag "$sample_id"

    container 'biocontainers/hisat2:v2.2.1_cv2'

    input:
    tuple val(sample_id), path(reads)

    output:
    path("${sample_id}/aligned.bam")

    script:
    """
    mkdir -p $sample_id
    hisat2 -x ${params.index} -1 ${reads[0]} -2 ${reads[1]} | samtools sort -o ${sample_id}/aligned.bam
    """
}

process count_reads {
    tag "$sample_id"

    container 'biocontainers/htseq:v0.11.2_cv3'

    input:
    path(bam)
    val sample_id

    output:
    path("${sample_id}/counts.txt")

    script:
    """
    htseq-count -f bam -r pos -s no -i gene_id ${bam} ${params.gtf} > ${sample_id}/counts.txt
    """
}

process deg_analysis {
    tag "DESeq2"

    container 'rocker/tidyverse:4.2.2'  // You can customize this to include DESeq2 & clusterProfiler

    input:
    path(counts)

    output:
    path("DEG_results.csv")
    path("PCA_plot.png")
    path("Volcano_plot.png")
    path("Heatmap.png")
    path("GO_plot.png")

    script:
    """
    Rscript ${params.rscript} ${counts}
    """
}

workflow {
    read_pairs_ch = Channel
        .fromFilePairs(params.reads, flat: true)
        .map { sample_id, reads -> tuple(sample_id, reads) }

    fastqc_out = read_pairs_ch | fastqc
    align_out = read_pairs_ch | hisat2_align

    count_out = align_out.map { bam -> tuple(bam.getParent().getName(), bam) } | count_reads

    count_out.collect() | deg_analysis
}

