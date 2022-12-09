
rule decoupler_expr:
    input:
        expr="input_data/rna_expr_shift.rds",
        meta="input_data/rna_meta.rds",
        network="input_data/pathway_commons_decoupler.rds"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    output:
        results_out="decoupler_workflow/results/{component}/decoupler_subset_results.rds",
        enr_weights="decoupler_workflow/results/{component}/decoupler_priori_weights.tsv"
    singularity:
        "library://yasharw/priori/decoupler_env_multiply_f-stat_corr"
    shell:
        """
        mkdir -p decoupler_workflow/results/{params.component}
        Rscript ./scripts/run_rna_benchmark.R {input.expr} {input.meta} {input.network} {output.results_out} {params.component} {output.enr_weights}
        """

rule generate_results_expr:
    input:
        results = "decoupler_workflow/results/{component}/decoupler_subset_results.rds"
    output:
        res_file = "decoupler_workflow/out_files/{component}_results.tsv",
        auc_file = "decoupler_workflow/out_files/{component}_summary_results.tsv"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://yasharw/priori/decoupler_env_multiply_f-stat_corr"
    shell:
        """
        mkdir -p decoupler_workflow/out_files
        Rscript ./scripts/generate_results.R {input.results} {output.res_file} {output.auc_file}
        """
