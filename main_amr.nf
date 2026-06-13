#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { AMR_ANALYSIS } from './subworkflows/amr_analysis'

workflow {
    strains_ch = Channel
        .fromPath("${params.prokka_dir}/*/", type: 'dir')
        .map { dir ->
            def sid = dir.name
            def faa = file("${dir}/${sid}.faa")
            tuple(sid, faa)
        }
        .filter { sid, faa -> faa.exists() }

    AMR_ANALYSIS(strains_ch)
}
