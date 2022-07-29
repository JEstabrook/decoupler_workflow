
rule decoupler_w_kd:
    input:
        expr="input_data/Knocktf_median_impute_expr.rds",
        meta="input_data/KnockTF_joined_meta.rds",
        network="input_data/KnockTF_intersected_regulon.rds"
    params:
        kd = lambda wildcards: "{}".format(wildcards.kd),
        cell = lambda wildcards: "{}".format(wildcards.cell),
        component = lambda wildcards: "{}".format(wildcards.component)
    output:
        temp_expr=temp("decoupler_workflow/results/{component}/{kd}/{cell}/decoupler_subset_expr.rds"),
        temp_meta=temp("decoupler_workflow/results/{component}/{kd}/{cell}/decoupler_subset_meta.rds"),
        temp_netw=temp("decoupler_workflow/results/{component}/{kd}/{cell}/decoupler_subset_network.rds"),
        results_out="decoupler_workflow/results/{component}/{kd}/{cell}/decoupler_subset_results.rds" 
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/results/{params.component}/{params.kd}/{params.cell}
        Rscript ./scripts/run_decoupler_kd_analysis.R {input.expr} {input.meta} {input.network} {params.kd} {params.cell} {output.results_out} {output.temp_expr} {output.temp_meta} {output.temp_netw} {params.component}
        """

rule decoupler_wo_kd:
    input:
        expr="input_data/Knocktf_median_impute_expr.rds",
        meta="input_data/KnockTF_joined_meta.rds",
        network="input_data/KnockTF_intersected_regulon.rds"
    params:
        cell = lambda wildcards: "{}".format(wildcards.cell),
        component = lambda wildcards: "{}".format(wildcards.component)
    output:
        temp_expr=temp("decoupler_workflow/kd_agnostic_results/{component}/{cell}/decoupler_subset_expr.rds"),
        temp_meta=temp("decoupler_workflow/kd_agnostic_results/{component}/{cell}/decoupler_subset_meta.rds"),
        temp_netw=temp("decoupler_workflow/kd_agnostic_results/{component}/{cell}/decoupler_subset_network.rds"),
        results_out="decoupler_workflow/kd_agnostic_results/{component}/{cell}/decoupler_subset_results.rds"
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/kd_agnostic_results/{params.cell}
        Rscript ./scripts/run_decoupler_wo_kd_analysis.R {input.expr} {input.meta} {input.network} {params.cell} {output.results_out} {output.temp_expr} {output.temp_meta} {output.temp_netw} {params.component}
        """

rule decoupler_w_kd_mod:
    input:
        expr="input_data/Knocktf_median_impute_expr.rds",
        meta="input_data/KnockTF_joined_meta.rds",
        network="input_data/KnockTF_intersected_regulon.rds"
    params:
        kd = lambda wildcards: "{}".format(wildcards.kd),
        cell = lambda wildcards: "{}".format(wildcards.cell),
        component = lambda wildcards: "{}".format(wildcards.component)
    output:
        temp_expr=temp("decoupler_workflow/mod_results/{component}/{kd}/{cell}/decoupler_subset_expr.rds"),
        temp_meta=temp("decoupler_workflow/mod_results/{component}/{kd}/{cell}/decoupler_subset_meta.rds"),
        temp_netw=temp("decoupler_workflow/mod_results/{component}/{kd}/{cell}/decoupler_subset_network.rds"),
        results_out="decoupler_workflow/mod_results/{component}/{kd}/{cell}/decoupler_subset_results.rds"
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/mod_results/{params.component}/{params.kd}/{params.cell}
        Rscript ./scripts/run_decoupler_kd_analysis_curves.R {input.expr} {input.meta} {input.network} {params.kd} {params.cell} {output.results_out} {output.temp_expr} {output.temp_meta} {output.temp_netw} {params.component}
        """

rule decoupler_wo_kd_mod:
    input:
        expr="input_data/Knocktf_median_impute_expr.rds",
        meta="input_data/KnockTF_joined_meta.rds",
        network="input_data/KnockTF_intersected_regulon.rds"
    params:
        cell = lambda wildcards: "{}".format(wildcards.cell),
        component = lambda wildcards: "{}".format(wildcards.component)
    output:
        temp_expr=temp("decoupler_workflow/mod_kd_agnostic_results/{component}/{cell}/decoupler_subset_expr.rds"),
        temp_meta=temp("decoupler_workflow/mod_kd_agnostic_results/{component}/{cell}/decoupler_subset_meta.rds"),
        temp_netw=temp("decoupler_workflow/mod_kd_agnostic_results/{component}/{cell}/decoupler_subset_network.rds"),
        results_out="decoupler_workflow/mod_kd_agnostic_results/{component}/{cell}/decoupler_subset_results.rds"
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/mod_kd_agnostic_results/{params.cell}
        Rscript ./scripts/run_decoupler_wo_kd_analysis_curves.R {input.expr} {input.meta} {input.network} {params.cell} {output.results_out} {output.temp_expr} {output.temp_meta} {output.temp_netw} {params.component}
        """

### Including sign

rule decoupler_w_kd_mod_sign:
    input:
        expr="input_data/Knocktf_median_impute_expr.rds",
        meta="input_data/KnockTF_joined_meta.rds",
        network="input_data/KnockTF_intersected_regulon.rds"
    params:
        kd = lambda wildcards: "{}".format(wildcards.kd),
        cell = lambda wildcards: "{}".format(wildcards.cell),
        component = lambda wildcards: "{}".format(wildcards.component)
    output:
        temp_expr=temp("decoupler_workflow/sign_mod_results/{component}/{kd}/{cell}/decoupler_subset_expr.rds"),
        temp_meta=temp("decoupler_workflow/sign_mod_results/{component}/{kd}/{cell}/decoupler_subset_meta.rds"),
        temp_netw=temp("decoupler_workflow/sign_mod_results/{component}/{kd}/{cell}/decoupler_subset_network.rds"),
        results_out="decoupler_workflow/sign_mod_results/{component}/{kd}/{cell}/decoupler_subset_results.rds",
        enr_weights="decoupler_workflow/sign_mod_results/{component}/{kd}/{cell}/decoupler_priori_weights.tsv"
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/sign_mod_results/{params.component}/{params.kd}/{params.cell}
        Rscript ./scripts/run_decoupler_kd_analysis_sign_curves.R {input.expr} {input.meta} {input.network} {params.kd} {params.cell} {output.results_out} {output.temp_expr} {output.temp_meta} {output.temp_netw} {params.component} {output.enr_weights}
        """

rule decoupler_wo_kd_mod_sign:
    input:
        expr="input_data/Knocktf_median_impute_expr.rds",
        meta="input_data/KnockTF_joined_meta.rds",
        network="input_data/KnockTF_intersected_regulon.rds"
    params:
        cell = lambda wildcards: "{}".format(wildcards.cell),
        component = lambda wildcards: "{}".format(wildcards.component)
    output:
        temp_expr=temp("decoupler_workflow/sign_mod_kd_agnostic_results/{component}/{cell}/decoupler_subset_expr.rds"),
        temp_meta=temp("decoupler_workflow/sign_mod_kd_agnostic_results/{component}/{cell}/decoupler_subset_meta.rds"),
        temp_netw=temp("decoupler_workflow/sign_mod_kd_agnostic_results/{component}/{cell}/decoupler_subset_network.rds"),
        results_out="decoupler_workflow/sign_mod_kd_agnostic_results/{component}/{cell}/decoupler_subset_results.rds",
        enr_weights="decoupler_workflow/sign_mod_kd_agnostic_results/{component}/{cell}/decoupler_priori_weights.tsv"
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/sign_mod_kd_agnostic_results/{params.cell}
        Rscript ./scripts/run_decoupler_wo_kd_analysis_sign_curves.R {input.expr} {input.meta} {input.network} {params.cell} {output.results_out} {output.temp_expr} {output.temp_meta} {output.temp_netw} {params.component} {output.enr_weights}
        """

##

rule generate_supple_fig_2_w_kd:
    input:
        results = "decoupler_workflow/results/{component}/{kd}/{cell}/decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/w_kd_plots/{component}/{kd}_{cell}_supplemental_figure2.pdf"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/w_kd_plots/{params.component}
        Rscript ./scripts/generate_heat_and_jaccard.R {input.results} {output.fig} 
        """

rule generate_supple_fig_2_wo_kd:
    input:
        results = "decoupler_workflow/kd_agnostic_results/{component}/{cell}/decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/wo_kd_plots/{component}/{cell}_supplemental_figure2.pdf"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/wo_kd_plots/{params.component}
        Rscript ./scripts/generate_heat_and_jaccard.R {input.results} {output.fig} 
        """

#### Adding mod

rule generate_supple_fig_2_w_kd_mod:
    input:
        results = "decoupler_workflow/mod_results/{component}/{kd}/{cell}/decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/w_kd_plots_mod/{component}/{kd}_{cell}_supplemental_figure2.pdf"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/w_kd_plots_mod/{params.component}
        Rscript ./scripts/generate_heat_and_jaccard.R {input.results} {output.fig}
        """

rule generate_supple_fig_2_wo_kd_mod:
    input:
        results = "decoupler_workflow/mod_kd_agnostic_results/{component}/{cell}/decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/wo_kd_plots_mod/{component}/{cell}_supplemental_figure2.pdf"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/wo_kd_plots_mod/{params.component}
        Rscript ./scripts/generate_heat_and_jaccard.R {input.results} {output.fig}
        """

####


rule generate_supple_fig_3_w_kd:
    input:
        results = "decoupler_workflow/results/{component}/{kd}/{cell}/decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/w_kd_plots/{component}/{kd}_{cell}_supplemental_figure3.pdf",
        csv_out = "decoupler_workflow/w_kd_plots/{component}/{kd}_{cell}_supplemental_table.csv"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/w_kd_plots/{params.component}
        Rscript ./scripts/generate_box_and_scatter.R {input.results} {output.fig} {output.csv_out}
        """

rule generate_supple_fig_3_wo_kd:
    input:
        results = "decoupler_workflow/kd_agnostic_results/{component}/{cell}/decoupler_subset_results.rds" 
    output:
        fig = "decoupler_workflow/wo_kd_plots/{component}/{cell}_supplemental_figure3.pdf",
        csv_out = "decoupler_workflow/wo_kd_plots/{component}/{cell}_supplemental_table.csv"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/wo_kd_plots/{params.component}
        Rscript ./scripts/generate_box_and_scatter.R {input.results} {output.fig} {output.csv_out}
        """


#### Adding mod


rule generate_supple_fig_3_w_kd_mod:
    input:
        results = "decoupler_workflow/mod_results/{component}/{kd}/{cell}/decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/w_kd_plots_mod/{component}/{kd}_{cell}_supplemental_figure3.pdf",
        csv_out = "decoupler_workflow/w_kd_plots_mod/{component}/{kd}_{cell}_supplemental_table.csv"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/w_kd_plots_mod/{params.component}
        Rscript ./scripts/generate_box_and_scatter.R {input.results} {output.fig} {output.csv_out}
        """

rule generate_supple_fig_3_wo_kd_mod:
    input:
        results = "decoupler_workflow/mod_kd_agnostic_results/{component}/{cell}/decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/wo_kd_plots_mod/{component}/{cell}_supplemental_figure3.pdf",
        csv_out = "decoupler_workflow/wo_kd_plots_mod/{component}/{cell}_supplemental_table.csv"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/wo_kd_plots_mod/{params.component}
        Rscript ./scripts/generate_box_and_scatter.R {input.results} {output.fig} {output.csv_out}
        """


#### Adding signed mod

rule generate_supple_fig_2_w_kd_sign_mod:
    input:
        results = "decoupler_workflow/sign_mod_results/{component}/{kd}/{cell}/decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/w_kd_plots_sign_mod/{component}/{kd}_{cell}_supplemental_figure2.pdf"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/w_kd_plots_sign_mod/{params.component}
        Rscript ./scripts/generate_heat_and_jaccard.R {input.results} {output.fig}
        """

rule generate_supple_fig_2_wo_kd_sign_mod:
    input:
        results = "decoupler_workflow/sign_mod_kd_agnostic_results/{component}/{cell}/decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/wo_kd_plots_sign_mod/{component}/{cell}_supplemental_figure2.pdf"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/wo_kd_plots_sign_mod/{params.component}
        Rscript ./scripts/generate_heat_and_jaccard.R {input.results} {output.fig}
        """
rule generate_supple_fig_3_w_kd_sign_mod:
    input:
        results = "decoupler_workflow/sign_mod_results/{component}/{kd}/{cell}/decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/w_kd_plots_sign_mod/{component}/{kd}_{cell}_supplemental_figure3.pdf",
        csv_out = "decoupler_workflow/w_kd_plots_sign_mod/{component}/{kd}_{cell}_supplemental_table.csv"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/w_kd_plots_sign_mod/{params.component}
        Rscript ./scripts/generate_box_and_scatter.R {input.results} {output.fig} {output.csv_out}
        """

rule generate_supple_fig_3_wo_kd_sign_mod:
    input:
        results = "decoupler_workflow/sign_mod_kd_agnostic_results/{component}/{cell}/decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/wo_kd_plots_sign_mod/{component}/{cell}_supplemental_figure3.pdf",
        csv_out = "decoupler_workflow/wo_kd_plots_sign_mod/{component}/{cell}_supplemental_table.csv"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/wo_kd_plots_sign_mod/{params.component}
        Rscript ./scripts/generate_box_and_scatter.R {input.results} {output.fig} {output.csv_out}
        """

#####

rule generate_results_w_kd:
    input:
        results = "decoupler_workflow/results/{component}/{kd}/{cell}/decoupler_subset_results.rds"
    output:
        res_file = "decoupler_workflow/out_files/{component}_{kd}_{cell}_w_kd_results.tsv",
        auc_file = "decoupler_workflow/out_files/{component}_{kd}_{cell}_w_kd_summary_results.tsv"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/out_files
        Rscript ./scripts/generate_kd_results.R {input.results} {output.res_file} {output.auc_file}
        """

rule generate_results_wo_kd:
    input:
        results = "decoupler_workflow/kd_agnostic_results/{component}/{cell}/decoupler_subset_results.rds"
    output:
        res_file = "decoupler_workflow/out_files/{component}_{cell}_wo_kd_results.tsv",
        auc_file = "decoupler_workflow/out_files/{component}_{cell}_wo_kd_summary_results.tsv"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/out_files
        Rscript ./scripts/generate_wo_kd_results.R {input.results} {output.res_file} {output.auc_file}
        """


##### mod

rule generate_results_w_kd_mod:
    input:
        results = "decoupler_workflow/mod_results/{component}/{kd}/{cell}/decoupler_subset_results.rds"
    output:
        res_file = "decoupler_workflow/mod_out_files/{component}_{kd}_{cell}_w_kd_results.tsv",
        auc_file = "decoupler_workflow/mod_out_files/{component}_{kd}_{cell}_w_kd_summary_results.tsv"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/mod_out_files
        Rscript ./scripts/generate_kd_results.R {input.results} {output.res_file} {output.auc_file}
        """

rule generate_results_wo_kd_mod:
    input:
        results = "decoupler_workflow/mod_kd_agnostic_results/{component}/{cell}/decoupler_subset_results.rds"
    output:
        res_file = "decoupler_workflow/mod_out_files/{component}_{cell}_wo_kd_results.tsv",
        auc_file = "decoupler_workflow/mod_out_files/{component}_{cell}_wo_kd_summary_results.tsv"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/mod_out_files
        Rscript ./scripts/generate_wo_kd_results.R {input.results} {output.res_file} {output.auc_file}
        """


##### mod with sign

rule generate_results_w_kd_mod_sign:
    input:
        results = "decoupler_workflow/sign_mod_results/{component}/{kd}/{cell}/decoupler_subset_results.rds"
    output:
        res_file = "decoupler_workflow/sign_mod_out_files/{component}_{kd}_{cell}_w_kd_results.tsv",
        auc_file = "decoupler_workflow/sign_mod_out_files/{component}_{kd}_{cell}_w_kd_summary_results.tsv"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/sign_mod_out_files
        Rscript ./scripts/generate_kd_results.R {input.results} {output.res_file} {output.auc_file}
        """

rule generate_results_wo_kd_mod_sign:
    input:
        results = "decoupler_workflow/sign_mod_kd_agnostic_results/{component}/{cell}/decoupler_subset_results.rds"
    output:
        res_file = "decoupler_workflow/sign_mod_out_files/{component}_{cell}_wo_kd_results.tsv",
        auc_file = "decoupler_workflow/sign_mod_out_files/{component}_{cell}_wo_kd_summary_results.tsv"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://jestabrook/regulon_enrichment/decoupler_env_slot_access_v2"
    shell:
        """
        mkdir -p decoupler_workflow/sign_mod_out_files
        Rscript ./scripts/generate_wo_kd_results.R {input.results} {output.res_file} {output.auc_file}
        """
