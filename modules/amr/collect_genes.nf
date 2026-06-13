process COLLECT_GENES {
    label 'python'
    publishDir "${params.results_dir}/blast_targets", mode: 'copy'

    input:
    path all_targets
    path ref_genes_dir

    output:
    path "*.combined.faa"

    script:
    """
    for gene in penA mtrR gyrA parC porB; do
        # Start with FA-1909 reference sequence
        cp ${ref_genes_dir}/\${gene}.faa \${gene}.combined.faa

        # Append matching gene sequence from each strain targets file
        python3 - \${gene} << 'PYEOF'
import sys, os
gene = sys.argv[1]
for fname in sorted(f for f in os.listdir('.') if f.endswith('_targets.faa')):
    with open(fname) as fh:
        header, seq_lines = None, []
        for line in fh:
            if line.startswith('>'):
                if header and gene in header:
                    with open(f'{gene}.combined.faa', 'a') as out:
                        out.write(header)
                        out.writelines(seq_lines)
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header and gene in header:
            with open(f'{gene}.combined.faa', 'a') as out:
                out.write(header)
                out.writelines(seq_lines)
PYEOF
        count=\$(grep -c ">" \${gene}.combined.faa)
        echo "INFO: \${gene} — \${count} sequences" >&2
    done
    """
}
