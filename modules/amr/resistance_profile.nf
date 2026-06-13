process RESISTANCE_PROFILE {
    label 'python'

    publishDir "${params.results_dir}/profiles", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.results_dir}/reports",  mode: 'copy', pattern: "*.html"

    input:
    path mutation_files
    path amrfinder_files

    output:
    path "resistance_profiles.tsv"
    path "amr_summary.html"

    script:
    """
    mkdir -p mut_dir amr_dir

    for f in *_mutations.tsv; do
        [ -f "\$f" ] && mv "\$f" mut_dir/
    done
    for f in *_amrfinder.tsv; do
        [ -f "\$f" ] && mv "\$f" amr_dir/
    done

    python3 ${projectDir}/bin/summarize_amr.py \
        --mutations-dir mut_dir \
        --amrfinder-dir amr_dir \
        --output-tsv    resistance_profiles.tsv \
        --output-html   amr_summary.html
    """
}
