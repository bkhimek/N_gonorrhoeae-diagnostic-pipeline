nextflow.enable.dsl = 2

include { AMRFINDERPLUS   } from '../modules/amr/amrfinderplus'
include { EXTRACT_TARGETS } from '../modules/amr/extract_targets'
include { COLLECT_GENES   } from '../modules/amr/collect_genes'
include { ALIGN_TARGETS   } from '../modules/amr/align_targets'
include { CALL_MUTATIONS  } from '../modules/amr/call_mutations'
include { RESISTANCE_PROFILE } from '../modules/amr/resistance_profile'

workflow AMR_ANALYSIS {

    take:
    strains_ch   // tuple val(strain_id), path(faa)

    main:

    // 1. AMRFinderPlus — parallel gene-level detection
    AMRFINDERPLUS(strains_ch)

    // 2. BLAST — extract target genes from each strain
    ref_genes = Channel.value(file("${params.reference_dir}/genes"))
    EXTRACT_TARGETS(strains_ch, ref_genes)

    // 3. Collect all per-strain target files; prepend reference; split by gene
    all_targets = EXTRACT_TARGETS.out.map { sid, f -> f }.collect()
    COLLECT_GENES(all_targets, ref_genes)

    // 4. Fan out — align each gene
    gene_ch = COLLECT_GENES.out
        .flatten()
        .map { f -> tuple(f.name.replace('.combined.faa', ''), f) }
    ALIGN_TARGETS(gene_ch)

    // 5. Call mutations per gene
    CALL_MUTATIONS(ALIGN_TARGETS.out)

    // 6. Merge all outputs into final profiles + HTML report
    mutations_files  = CALL_MUTATIONS.out.map { g, f -> f }.collect()
    amrfinder_files  = AMRFINDERPLUS.out.map { s, f -> f }.collect()
    RESISTANCE_PROFILE(mutations_files, amrfinder_files)

    emit:
    profiles = RESISTANCE_PROFILE.out[0]
    report   = RESISTANCE_PROFILE.out[1]
}
