nextflow.enable.dsl=2

out_dir = file(params.outdir)
mode = params.publish_dir_mode

process CMS_CALLER {
    publishDir "${out_dir}", mode: "$mode"

    input:
    path(ch_out_count_mat)

    output:
    path("*.csv"), emit: cmscaller_results

    script:
    """
    CMScaller.R --count_mat=${ch_out_count_mat}   
    """

   
}