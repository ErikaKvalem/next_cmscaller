nextflow.enable.dsl=2

out_dir = file(params.outdir)
mode = params.publish_dir_mode

process EXTRACT_ADATA {
    publishDir "${out_dir}", mode: "$mode"

    input:
    path(adata_ch)

    output:
    path("*.mtx"), emit: counts_matrix

    script:
    """
    convert_extract_adata.py --adata=${adata_ch}
    """


}