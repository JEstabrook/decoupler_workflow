
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

rule kd_decoupler_base_plots:
    input:
        bench_results = "decoupler_workflow/results/{kd}/{cell}/decoupler_subset_results.rds"
    output:
        roc_plot = "decoupler_workflow/w_kd_plots/{kd}_{cell}_base_plot_ROC.pdf",
        pr_plot = "decoupler_workflow/w_kd_plots/{kd}_{cell}_base_plot_PR.pdf",
        auroc_heat = "decoupler_workflow/w_kd_plots/{kd}_{cell}_base_plot_AUROC_HEAT.pdf",
        pr_heat = "decoupler_workflow/w_kd_plots/{kd}_{cell}_base_plot_PR_HEAT.pdf",
    shell:
        """
        mkdir -p decoupler_workflow/w_kd_plots/
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/generate_decoupler_base_plots.R {input.bench_results} {output.roc_plot} {output.pr_plot} {output.auroc_heat} {output.pr_heat}
        """

rule wo_kd_decoupler_base_plots:
    input:
        bench_results = "decoupler_workflow/kd_agnostic_results/{cell}/decoupler_subset_results.rds"
    output:
        roc_plot = "decoupler_workflow/wo_kd_plots/{cell}_base_plot_ROC.pdf",
        pr_plot = "decoupler_workflow/wo_kd_plots/{cell}_base_plot_PR.pdf",
        auroc_heat = "decoupler_workflow/wo_kd_plots/{cell}_base_plot_AUROC_HEAT.pdf",
        pr_heat = "decoupler_workflow/wo_kd_plots/{cell}_base_plot_PR_HEAT.pdf",
    shell:
        """
        mkdir -p decoupler_workflow/wo_kd_plots/
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/generate_decoupler_base_plots.R {input.bench_results} {output.roc_plot} {output.pr_plot} {output.auroc_heat} {output.pr_heat}
        """


rule down_sample_decoupler_base_plots:
    input:
        bench_results = "decoupler_workflow/downsample/{nokdcontrol}_{downsample}_{cell}_decoupler_subset_results.rds"
    output:
        roc_plot = "decoupler_workflow/nokddownsample/baseplots/{nokdcontrol}__{downsample}_{cell}_base_plot_ROC.pdf",
        pr_plot = "decoupler_workflow/nokddownsample/baseplots/{nokdcontrol}__{downsample}_{cell}_base_plot_PR.pdf",
        auroc_heat = "decoupler_workflow/nokddownsample/baseplots/{nokdcontrol}__{downsample}_{cell}_base_plot_AUROC_HEAT.pdf",
        pr_heat = "decoupler_workflow/nokddownsample/baseplots/{nokdcontrol}__{downsample}_{cell}_base_plot_PR_HEAT.pdf",
    shell:
        """
        mkdir -p decoupler_workflow/nokddownsample/baseplots
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/generate_decoupler_base_plots.R {input.bench_results} {output.roc_plot} {output.pr_plot} {output.auroc_heat} {output.pr_heat}
        """

rule down_sample_kd_decoupler_base_plots:
    input:
        bench_results = "decoupler_workflow/downsample/{control}_{downsample}_{kd}_{cell}_decoupler_subset_results.rds"
    output:
        roc_plot = "decoupler_workflow/wkddownsample/baseplots/{control}_{kd}_{downsample}_{cell}_base_plot_ROC.pdf",
        pr_plot = "decoupler_workflow/wkddownsample/baseplots/{control}_{kd}_{downsample}_{cell}_base_plot_PR.pdf",
        auroc_heat = "decoupler_workflow/wkddownsample/baseplots/{control}_{kd}_{downsample}_{cell}_base_plot_AUROC_HEAT.pdf",
        pr_heat = "decoupler_workflow/wkddownsample/baseplots/{control}_{kd}_{downsample}_{cell}_base_plot_PR_HEAT.pdf",
    shell:
        """
        mkdir -p decoupler_workflow/wkddownsample/baseplots
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/generate_decoupler_base_plots.R {input.bench_results} {output.roc_plot} {output.pr_plot} {output.auroc_heat} {output.pr_heat}
        """


rule kd_downsample_all:
    input:
        bench_results = "decoupler_workflow/results/{kd}/{cell}/decoupler_subset_results.rds"
    params:
        down = lambda wildcards: "{}".format(wildcards.downsample),
    output:
        results_out = "decoupler_workflow/downsample/kdallcontrol_{downsample}_{kd}_{cell}_decoupler_subset_results.rds"
    shell:
        """
        mkdir -p decoupler_workflow/downsample/
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/kd_all_control_sampling.R {input.bench_results} {params.down} {output.results_out}
        """

rule kd_downsample_only:
    input:
        bench_results = "decoupler_workflow/results/{kd}/{cell}/decoupler_subset_results.rds"
    params:
        down = lambda wildcards: "{}".format(wildcards.downsample),
    output:
        results_out = "decoupler_workflow/downsample/kdonlycontrol_{downsample}_{kd}_{cell}_decoupler_subset_results.rds"
    shell:
        """
        mkdir -p decoupler_workflow/downsample/
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/kd_specific_sampling.R {input.bench_results} {params.down} {output.results_out}
        """

rule nokd_downsample_all:
    input:
        bench_results = "decoupler_workflow/kd_agnostic_results/{cell}/decoupler_subset_results.rds" 
    params:
        down = lambda wildcards: "{}".format(wildcards.downsample),
    output:
        results_out = "decoupler_workflow/downsample/nokdallcontrol_{downsample}_{cell}_decoupler_subset_results.rds"
    shell:
        """
        mkdir -p decoupler_workflow/downsample/
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/kd_all_control_sampling.R {input.bench_results} {params.down} {output.results_out}
        """

rule nokd_downsample_only:
    input:
        bench_results = "decoupler_workflow/kd_agnostic_results/{cell}/decoupler_subset_results.rds"
    params:
        down = lambda wildcards: "{}".format(wildcards.downsample),
    output:
        results_out = "decoupler_workflow/downsample/nokdonlycontrol_{downsample}_{cell}_decoupler_subset_results.rds"
    shell:
        """
        mkdir -p decoupler_workflow/downsample/
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/kd_specific_sampling.R {input.bench_results} {params.down} {output.results_out}
        """

rule down_sample_generate_supple_fig_2_kdallcontrol:
    input:
        results = "decoupler_workflow/downsample/kdallcontrol_{downsample}_{kd}_{cell}_decoupler_subset_results.rds" 
    output:
        fig = "decoupler_workflow/downsample/w_kd_plots/kdallcontrol_{downsample}_{kd}_{cell}_supplemental_figure2.pdf"
    shell:
        """
        mkdir -p decoupler_workflow/downsample/w_kd_plots
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/generate_heat_and_jaccard.R {input.results} {output.fig}
        """

rule down_sample_generate_supple_fig_2_kdonlycontrol:
    input:
        results = "decoupler_workflow/downsample/kdonlycontrol_{downsample}_{kd}_{cell}_decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/downsample/w_kd_plots/kdonlycontrol_{downsample}_{kd}_{cell}_supplemental_figure2.pdf"
    shell:
        """
        mkdir -p decoupler_workflow/downsample/w_kd_plots
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/generate_heat_and_jaccard.R {input.results} {output.fig}
        """

rule down_sample_generate_supple_fig_2_nokdallcontrol:
    input:
        results = "decoupler_workflow/downsample/nokdallcontrol_{downsample}_{cell}_decoupler_subset_results.rds" 
    output:
        fig = "decoupler_workflow/downsample/wo_kd_plots/nokdallcontrol_{downsample}_{cell}_supplemental_figure2.pdf"
    shell:
        """
        mkdir -p decoupler_workflow/downsample/wo_kd_plots
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/generate_heat_and_jaccard.R {input.results} {output.fig}
        """

rule down_sample_generate_supple_fig_2_nokdonlycontrol:
    input:
        results = "decoupler_workflow/downsample/nokdonlycontrol_{downsample}_{cell}_decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/downsample/wo_kd_plots/nokdonlycontrol_{downsample}_{cell}_supplemental_figure2.pdf"
    shell:
        """
        mkdir -p decoupler_workflow/downsample/wo_kd_plots
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/generate_heat_and_jaccard.R {input.results} {output.fig}
        """
### figure 3

rule down_sample_generate_supple_fig_3_kdallcontrol:
    input:
        results = "decoupler_workflow/downsample/kdallcontrol_TRUE_{kd}_{cell}_decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/downsample/w_kd_plots/kdallcontrol_TRUE_{kd}_{cell}_supplemental_figure3.pdf",
        csv_out = "decoupler_workflow/downsample/w_kd_plots/kdallcontrol_TRUE_{kd}_{cell}_supplemental_table.csv"
    shell:
        """
        mkdir -p decoupler_workflow/downsample/w_kd_plots
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/generate_box_and_scatter.R {input.results} {output.fig} {output.csv_out}
        """

rule down_sample_generate_supple_fig_3_kdonlycontrol:
    input:
        results = "decoupler_workflow/downsample/kdonlycontrol_TRUE_{kd}_{cell}_decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/downsample/w_kd_plots/kdonlycontrol_TRUE_{kd}_{cell}_supplemental_figure3.pdf",
        csv_out = "decoupler_workflow/downsample/w_kd_plots/kdonlycontrol_TRUE_{kd}_{cell}_supplemental_table.csv"
    shell:
        """
        mkdir -p decoupler_workflow/downsample/w_kd_plots
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/generate_box_and_scatter.R {input.results} {output.fig} {output.csv_out}
        """

rule down_sample_generate_supple_fig_3_nokdallcontrol:
    input:
        results = "decoupler_workflow/downsample/nokdallcontrol_TRUE_{cell}_decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/downsample/wo_kd_plots/nokdallcontrol_TRUE_{cell}_supplemental_figure3.pdf",
        csv_out = "decoupler_workflow/downsample/wo_kd_plots/nokdallcontrol_TRUE_{cell}_supplemental_table.csv"
    shell:
        """
        mkdir -p decoupler_workflow/downsample/wo_kd_plots
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/generate_box_and_scatter.R {input.results} {output.fig} {output.csv_out}
        """

rule down_sample_generate_supple_fig_3_nokdonlycontrol:
    input:
        results = "decoupler_workflow/downsample/nokdonlycontrol_TRUE_{cell}_decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/downsample/wo_kd_plots/nokdonlycontrol_TRUE_{cell}_supplemental_figure3.pdf",
        csv_out = "decoupler_workflow/downsample/wo_kd_plots/nokdonlycontrol_TRUE_{cell}_supplemental_table.csv"
    shell:
        """
        mkdir -p decoupler_workflow/downsample/wo_kd_plots
        set +u
        source deactivate
        source activate decoupler_env
        Rscript ./scripts/generate_box_and_scatter.R {input.results} {output.fig} {output.csv_out}
        """

