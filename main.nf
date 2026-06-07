params.gffs           = "/home/himek/s_pyogenes_project/prokka_output/*/*.gff"
params.blast_db       = "/home/himek/s_pyogenes_project/streptococcus_non_pyogenes_db"
params.core_threshold = 0.99
params.identity_thr   = 98.0
params.min_pident     = 85.0
params.min_qcovs      = 80.0
params.panaroo_results = "results/panaroo"

process SPLIT_AND_ALIGN {
    container 'strep-analysis:1.0'
    publishDir "results/split", mode: 'copy'

    input:
    path gpa
    path gene_data
    path cds

    output:
    path "split", emit: split_dir

    script:
    """
    mkdir -p split
    split_core100_from_panaroo.py \
        --gpa ${gpa} \
        --gene_data ${gene_data} \
        --cds ${cds} \
        --outdir split \
        --min_presence ${params.core_threshold} \
        --threads 8
    """
}

process CALC_IDENTITY {
    container 'strep-analysis:1.0'
    publishDir "results/inclusivity", mode: 'copy'

    input:
    path split_dir

    output:
    path "identity.tsv", emit: identity

    script:
    """
    calculate_identity_with_names.py \
        --split_dir ${split_dir} \
        --output identity.tsv
    """
}

process FILTER_CANDIDATES {
    container 'strep-analysis:1.0'
    publishDir "results/inclusivity", mode: 'copy'

    input:
    path identity

    output:
    path "candidates.tsv", emit: candidates

    script:
    """
    filter_high_identity_genes.py \
        --input ${identity} \
        --output candidates_raw.tsv \
        --thr ${params.identity_thr}
    grep -v [.]raw candidates_raw.tsv > candidates.tsv
    """
}

process BUILD_CONSENSUS {
    container 'strep-analysis:1.0'
    publishDir "results/consensus", mode: 'copy'

    input:
    path split_dir
    path candidates

    output:
    path "consensus.fasta", emit: consensus

    script:
    """
    build_consensus_from_split.py \
        --split_dir ${split_dir} \
        --out_fasta consensus.fasta \
        --genes_txt ${candidates}
    """
}

process BLAST_EXCLUSIVITY {
    container 'strep-analysis:1.0'
    publishDir "results/exclusivity", mode: 'copy'

    input:
    path consensus

    output:
    path "blast_results.tsv", emit: blast_tsv

    script:
    """
    blastn \
        -query ${consensus} \
        -db ${params.blast_db} \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs sscinames staxids stitle" \
        -perc_identity 1 \
        -qcov_hsp_perc 50 \
        -num_threads 8 \
        -out blast_results.tsv
    """
}

process SUMMARISE_EXCLUSIVITY {
    container 'strep-analysis:1.0'
    publishDir "results/exclusivity", mode: 'copy'

    input:
    path blast_tsv

    output:
    path "exclusivity.tsv", emit: exclusivity

    script:
    """
    summarize_blast_exclusivity.py \
        --blast_tsv ${blast_tsv} \
        --out_tsv exclusivity.tsv \
        --min_pident ${params.min_pident} \
        --min_qcovs ${params.min_qcovs}
    """
}

process EXTRACT_PASS {
    container 'strep-analysis:1.0'
    publishDir "results/pass_fastas", mode: 'copy'

    input:
    path consensus
    path exclusivity

    output:
    path "pass/*.fasta", emit: pass_fastas

    script:
    """
    mkdir -p pass
    extract_pass_consensus.py \
        --consensus ${consensus} \
        --exclusivity_tsv ${exclusivity} \
        --outdir pass
    """
}

workflow {
    gpa       = Channel.fromPath("${params.panaroo_results}/gene_presence_absence.csv")
    gene_data = Channel.fromPath("${params.panaroo_results}/gene_data.csv")
    cds       = Channel.fromPath("${params.panaroo_results}/combined_DNA_CDS.fasta")

    SPLIT_AND_ALIGN(gpa, gene_data, cds)
    CALC_IDENTITY(SPLIT_AND_ALIGN.out.split_dir)
    FILTER_CANDIDATES(CALC_IDENTITY.out.identity)
    BUILD_CONSENSUS(SPLIT_AND_ALIGN.out.split_dir, FILTER_CANDIDATES.out.candidates)
    BLAST_EXCLUSIVITY(BUILD_CONSENSUS.out.consensus)
    SUMMARISE_EXCLUSIVITY(BLAST_EXCLUSIVITY.out.blast_tsv)
    EXTRACT_PASS(BUILD_CONSENSUS.out.consensus, SUMMARISE_EXCLUSIVITY.out.exclusivity)
}
