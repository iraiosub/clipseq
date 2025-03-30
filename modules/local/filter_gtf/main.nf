process FILTER_GTF_BY_TRANSCRIPT {
    tag "$gtf"
    label "process_single"

    conda "bioconda::pyranges=0.1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pyranges:0.1.2--pyhdfd78af_1':
        'quay.io/biocontainers/pyranges:0.1.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(gtf)
    tuple val(meta), path(transcript)

    output:
    tuple val(meta), path("*_representative_transcript_filtered.gtf")       ,emit: filtered_gtf
    path  "*.log"                                                           ,emit: log
    path  "versions.yml"                                                    ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    process_name = task.process
    output       = task.ext.output ?: "${gtf.simpleName}_representative_transcript"
    template 'filter_gtf_by_transcripts.py'
}
