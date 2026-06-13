process EXTRACT_TARGETS {
    tag "$strain_id"
    label 'blast'
    publishDir "${params.results_dir}/blast_targets", mode: 'copy'

    input:
    tuple val(strain_id), path(faa)
    path reference_genes_dir

    output:
    tuple val(strain_id), path("${strain_id}_targets.faa")

    script:
    """
    makeblastdb -in ${faa} -dbtype prot -out strain_db 2>/dev/null

    > ${strain_id}_targets.faa

    for gene in penA mtrR gyrA parC porB; do
        ref="${reference_genes_dir}/\${gene}.faa"

        blastp \\
            -query \$ref \\
            -db strain_db \\
            -outfmt "6 sseqid pident qcovs evalue bitscore" \\
            -max_target_seqs 5 \\
            -evalue 1e-10 \\
            -num_threads ${task.cpus} \\
            -out blast_\${gene}.txt

        hit=\$(awk -v pid=${params.blast_pident} -v cov=${params.blast_qcov} \\
            '\$2 >= pid && \$3 >= cov {print \$1; exit}' blast_\${gene}.txt)

        if [ -n "\$hit" ]; then
            samtools faidx ${faa} "\$hit" | \\
                sed "1s/.*/>${strain_id}__\${gene}/" >> ${strain_id}_targets.faa
        else
            echo "WARNING: no hit for \${gene} in ${strain_id}" >&2
        fi
    done

    found=\$(grep -c ">" ${strain_id}_targets.faa || true)
    echo "INFO: ${strain_id} — \${found}/5 genes extracted" >&2
    """
}
