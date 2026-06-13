/*
 * call_mutations.nf
 *
 * Runs parse_mutations.py on a per-gene MSA to extract resistance mutations.
 * FA-1909 (reference) must be the first sequence in the alignment.
 * Outputs a per-strain, per-gene mutation table.
 *
 * Input:  tuple val(gene_name), path(aligned_faa)
 * Output: tuple val(gene_name), path("${gene_name}_mutations.tsv")
 */

process CALL_MUTATIONS {
    tag "$gene_name"
    label 'python'

    publishDir "${params.results_dir}/mutations", mode: 'copy'

    input:
    tuple val(gene_name), path(aligned_faa)

    output:
    tuple val(gene_name), path("${gene_name}_mutations.tsv")

    script:
    def qrdr_args = ""
    if (gene_name == "gyrA") {
        qrdr_args = "--qrdr-start ${params.qrdr_gyra_start} --qrdr-end ${params.qrdr_gyra_end}"
    } else if (gene_name == "parC") {
        qrdr_args = "--qrdr-start ${params.qrdr_parc_start} --qrdr-end ${params.qrdr_parc_end}"
    }
    """
    python3 ${projectDir}/bin/parse_mutations.py \\
        --alignment  ${aligned_faa} \\
        --gene       ${gene_name} \\
        --reference  ${params.reference_strain} \\
        ${qrdr_args} \\
        --output     ${gene_name}_mutations.tsv
    """
}
