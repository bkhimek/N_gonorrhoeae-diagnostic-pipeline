/*
 * amrfinderplus.nf
 *
 * Runs AMRFinderPlus in protein-only mode on each strain's Prokka .faa.
 * Detects resistance genes and protein-level point mutations (POINTP).
 *
 * Note: GFF + nucleotide mode is incompatible with Prokka GFF3 format
 * (attribute mismatch: Prokka uses ID=/locus_tag=, AMRFinder expects protein_id=).
 * Protein-only mode still captures gene-level hits and POINTP mutations,
 * which is sufficient for our cross-check role. QRDR point mutations
 * (GyrA, ParC) are handled by the custom alignment + parse_mutations.py.
 *
 * Input:  tuple val(strain_id), path(faa)
 * Output: tuple val(strain_id), path("${strain_id}_amrfinder.tsv")
 */

process AMRFINDERPLUS {
    tag "$strain_id"
    label 'amrfinder'

    publishDir "${params.results_dir}/amrfinderplus", mode: 'copy'

    input:
    tuple val(strain_id), path(faa)

    output:
    tuple val(strain_id), path("${strain_id}_amrfinder.tsv")

    script:
    """
    amrfinder \\
        --protein   ${faa} \\
        --organism  Neisseria_gonorrhoeae \\
        --name      ${strain_id} \\
        --output    ${strain_id}_amrfinder.tsv \\
        --database  ${params.amrfinder_db} \\
        --threads   ${task.cpus}
    """
}
