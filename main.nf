#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { EXTRACT_ADATA} from "./modules/01_convert_extract_adata"
include { CMS_CALLER} from "./modules/02_cmscaller"

workflow {
    adata_ch = Channel.fromPath(params.input_path)

    EXTRACT_ADATA(adata_ch)
    ch_out_count_mat = EXTRACT_ADATA.out.counts_matrix
    ch_out_count_mat.view()
    CMS_CALLER(ch_out_count_mat)    
}
