
rule decoupler_w_kd:
    input:
        expr="input_data/Knocktf_median_impute_expr.rds",
        meta="input_data/KnockTF_joined_meta.rds",
        network="input_data/KnockTF_intersected_regulon.rds"
    params:
        kd = lambda wildcards: "{}".format(wildcards.kd),
        cell = lambda wildcards: "{}".format(wildcards.cell)
    output:
        temp_expr=temp("decoupler_workflow/results/{kd}/{cell}/decoupler_subset_expr.rds"),
        temp_meta=temp("decoupler_workflow/results/{kd}/{cell}/decoupler_subset_meta.rds"),
        temp_netw=temp("decoupler_workflow/results/{kd}/{cell}/decoupler_subset_network.rds"),
        results_out="decoupler_workflow/results/{kd}/{cell}/decoupler_subset_results.rds" 
    shell:
        """
        mkdir -p decoupler_workflow/results/{params.kd}/{params.cell}
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/run_decoupler_kd_analysis.R {input.expr} {input.meta} {input.network} {params.kd} {params.cell} {output.results_out} {output.temp_expr} {output.temp_meta} {output.temp_netw}
        """

rule decoupler_wo_kd:
    input:
        expr="input_data/Knocktf_median_impute_expr.rds",
        meta="input_data/KnockTF_joined_meta.rds",
        network="input_data/KnockTF_intersected_regulon.rds"
    params:
        cell = lambda wildcards: "{}".format(wildcards.cell),
    output:
        temp_expr=temp("decoupler_workflow/kd_agnostic_results/{cell}/decoupler_subset_expr.rds"),
        temp_meta=temp("decoupler_workflow/kd_agnostic_results/{cell}/decoupler_subset_meta.rds"),
        temp_netw=temp("decoupler_workflow/kd_agnostic_results/{cell}/decoupler_subset_network.rds"),
        results_out="decoupler_workflow/kd_agnostic_results/{cell}/decoupler_subset_results.rds"
    shell:
        """
        mkdir -p decoupler_workflow/kd_agnostic_results/{params.cell}
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/run_decoupler_wo_kd_analysis.R {input.expr} {input.meta} {input.network} {params.cell} {output.results_out} {output.temp_expr} {output.temp_meta} {output.temp_netw}
        """

rule generate_supple_fig_2_w_kd:
    input:
        results = "decoupler_workflow/results/{kd}/{cell}/decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/w_kd_plots/{kd}_{cell}_supplemental_figure2.pdf"
    shell:
        """
        mkdir -p decoupler_workflow/w_kd_plots
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/generate_heat_and_jaccard.R {input.results} {output.fig} 
        """

rule generate_supple_fig_2_wo_kd:
    input:
        results = "decoupler_workflow/kd_agnostic_results/{cell}/decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/wo_kd_plots/{cell}_supplemental_figure2.pdf"
    shell:
        """
        mkdir -p decoupler_workflow/wo_kd_plots
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/generate_heat_and_jaccard.R {input.results} {output.fig} 
        """

rule generate_supple_fig_3_w_kd:
    input:
        results = "decoupler_workflow/results/{kd}/{cell}/decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/w_kd_plots/{kd}_{cell}_supplemental_figure3.pdf",
        csv_out = "decoupler_workflow/w_kd_plots/{kd}_{cell}_supplemental_table.csv"
    shell:
        """
        mkdir -p decoupler_workflow/w_kd_plots
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/generate_box_and_scatter.R {input.results} {output.fig} {output.csv_out}
        """

rule generate_supple_fig_3_wo_kd:
    input:
        results = "decoupler_workflow/kd_agnostic_results/{cell}/decoupler_subset_results.rds" 
    output:
        fig = "decoupler_workflow/wo_kd_plots/{cell}_supplemental_figure3.pdf",
        csv_out = "decoupler_workflow/wo_kd_plots/{cell}_supplemental_table.csv"
    shell:
        """
        mkdir -p decoupler_workflow/wo_kd_plots
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/generate_box_and_scatter.R {input.results} {output.fig} {output.csv_out}
        """
