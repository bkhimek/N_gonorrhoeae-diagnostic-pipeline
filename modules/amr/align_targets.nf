process ALIGN_TARGETS {
    tag "$gene_name"
    label 'mafft'
    publishDir "${params.results_dir}/alignments", mode: 'copy'

    input:
    tuple val(gene_name), path(combined_fasta)

    output:
    tuple val(gene_name), path("${gene_name}_aligned.faa")

    script:
    """
    mafft \\
        --auto \\
        --thread ${task.cpus} \\
        --reorder \\
        ${combined_fasta} \\
        > ${gene_name}_aligned.faa
    """
}
