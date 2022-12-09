
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
        "library://yasharw/priori/decoupler_env_no_delta_force_regulon_network"
    shell:
        """
        mkdir -p decoupler_workflow/results/{params.component}
        Rscript ./scripts/run_rna_benchmark.R {input.expr} {input.meta} {input.network} {output.results_out} {params.component} {output.enr_weights}
        """

rule generate_supple_fig_2:
    input:
        results = "decoupler_workflow/results/{component}/decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/baseline_RNA_experiment/{component}/supplemental_figure2.pdf"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://yasharw/priori/decoupler_env_no_delta_force_regulon_network"
    shell:
        """
        mkdir -p decoupler_workflow/baseline_RNA_experiment/
        Rscript ./scripts/generate_heat_and_jaccard.R {input.results} {output.fig}
        """

rule generate_supple_fig_3:
    input:
        results = "decoupler_workflow/results/{component}/decoupler_subset_results.rds"
    output:
        fig = "decoupler_workflow/baseline_RNA_experiment/{component}/supplemental_figure3.pdf",
        csv_out = "decoupler_workflow/baseline_RNA_experiment/{component}/supplemental_table.csv"
    params:
        component = lambda wildcards: "{}".format(wildcards.component)
    singularity:
        "library://yasharw/priori/decoupler_env_no_delta_force_regulon_network"
    shell:
        """
        mkdir -p decoupler_workflow/baseline_RNA_experiment/
        Rscript ./scripts/generate_box_and_scatter.R {input.results} {output.fig} {output.csv_out}
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
        "library://yasharw/priori/decoupler_env_no_delta_force_regulon_network"
    shell:
        """
        mkdir -p decoupler_workflow/out_files
        Rscript ./scripts/generate_results.R {input.results} {output.res_file} {output.auc_file}
        """

# Noise
rule generate_gaussian_noise:
    input:
        expr="input_data/rna_expr_shift.rds"
    params:
        gauss_sd = lambda wildcards: "{}".format(wildcards.gauss_sd),
        n_rep = lambda wildcards: "{}".format(wildcards.n_rep)
    output:
        expr_noise_out="decoupler_workflow/rna_gaussian_noise/expr_rna_noise_sd{gauss_sd}_rep{n_rep}.rds"
    singularity:
        "library://yasharw/priori/decoupler_env_no_delta_force_regulon_network"
    shell:
        """
        mkdir -p decoupler_workflow/rna_gaussian_noise
        Rscript ./scripts/generate_gaussian_noise.R {input.expr} {output.expr_noise_out} {params.gauss_sd} {params.n_rep}
        """

rule run_rna_gaussian_noise:
    input:
        expr="decoupler_workflow/rna_gaussian_noise/expr_rna_noise_sd{gauss_sd}_rep{n_rep}.rds",
        meta="input_data/rna_meta.rds",
        network="input_data/pathway_commons_decoupler.rds"
    params:
        component = lambda wildcards: "{}".format(wildcards.component),
        gauss_sd = lambda wildcards: "{}".format(wildcards.gauss_sd),
        n_rep = lambda wildcards: "{}".format(wildcards.n_rep)
    output:
        results_out="decoupler_workflow/rna_gaussian_noise/{component}/rna_gaussian_noise_sd{gauss_sd}_rep{n_rep}.rds",
        res_file="decoupler_workflow/rna_gaussian_noise/{component}/out_files/rna_gaussian_noise_sd{gauss_sd}_rep{n_rep}_results.tsv",
        auc_file="decoupler_workflow/rna_gaussian_noise/{component}/out_files/rna_gaussian_noise_sd{gauss_sd}_rep{n_rep}_summary_results.tsv"
    singularity:
        "library://yasharw/priori/decoupler_env_no_delta_force_regulon_network"
    shell:
        """
        mkdir -p decoupler_workflow/rna_gaussian_noise/{params.component}
        mkdir -p decoupler_workflow/rna_gaussian_noise/{params.component}/out_files
        Rscript ./scripts/run_rna_gaussian_noise.R {input.expr} {input.meta} {input.network} {output.results_out} {output.res_file} {output.auc_file} {params.gauss_sd} {params.n_rep} {params.component}
        """

rule randomize_edges:
    input:
        expr="input_data/rna_expr_shift.rds",
        network="input_data/pathway_commons_decoupler.rds"
    params:
        noise_perc = lambda wildcards: "{}".format(wildcards.noise_perc),
        n_rep = lambda wildcards: "{}".format(wildcards.n_rep)
    output:
        netw_noise_out="decoupler_workflow/random_edges/netw_random_edges_frac{noise_perc}_rep{n_rep}.rds",
        netw_sec_intx_noise_out="decoupler_workflow/random_edges/netw_random_edges_frac{noise_perc}_rep{n_rep}_secondary_intx.pkl"
    # singularity:
    #     "library://yasharw/priori/decoupler_env_no_delta_force_regulon_network"
    shell:
        """
        mkdir -p decoupler_workflow/random_edges
        Rscript ./scripts/randomize_edges.R {input.expr} {input.network} {output.netw_noise_out} {output.netw_sec_intx_noise_out} {params.noise_perc} {params.n_rep}
        """

rule run_random_edges:
    input:
        expr="input_data/rna_expr_shift.rds",
        meta="input_data/rna_meta.rds",
        network="decoupler_workflow/random_edges/netw_random_edges_frac{noise_perc}_rep{n_rep}.rds",
        sec_intx_network="decoupler_workflow/random_edges/netw_random_edges_frac{noise_perc}_rep{n_rep}_secondary_intx.pkl"
    params:
        component = lambda wildcards: "{}".format(wildcards.component),
        noise_perc = lambda wildcards: "{}".format(wildcards.noise_perc),
        n_rep = lambda wildcards: "{}".format(wildcards.n_rep)
    output:
        results_out="decoupler_workflow/random_edges/{component}/random_edges_frac{noise_perc}_rep{n_rep}.rds",
        res_file="decoupler_workflow/random_edges/{component}/out_files/random_edges_frac{noise_perc}_rep{n_rep}_results.tsv",
        auc_file="decoupler_workflow/random_edges/{component}/out_files/random_edges_frac{noise_perc}_rep{n_rep}_summary_results.tsv"
    singularity:
        "library://yasharw/priori/decoupler_env_no_delta_force_regulon_network"
    shell:
        """
        mkdir -p decoupler_workflow/random_edges/{params.component}
        mkdir -p decoupler_workflow/random_edges/{params.component}/out_files
        Rscript ./scripts/run_random_edges.R {input.expr} {input.meta} {input.network} {input.sec_intx_network} {output.results_out} {output.res_file} {output.auc_file} {params.noise_perc} {params.n_rep} {params.component}
        """

rule shuffle_edges:
    input:
        network="input_data/pathway_commons_decoupler.rds"
    params:
        noise_perc = lambda wildcards: "{}".format(wildcards.noise_perc),
        n_rep = lambda wildcards: "{}".format(wildcards.n_rep)
    output:
        netw_noise_out="decoupler_workflow/shuffle_edges/netw_shuffle_edges_frac{noise_perc}_rep{n_rep}.rds",
        netw_sec_intx_noise_out="decoupler_workflow/shuffle_edges/netw_shuffle_edges_frac{noise_perc}_rep{n_rep}_secondary_intx.pkl"
    # singularity:
    #     "library://yasharw/priori/decoupler_env_no_delta_force_regulon_network"
    shell:
        """
        mkdir -p decoupler_workflow/shuffle_edges
        Rscript ./scripts/shuffle_edges.R {input.network} {output.netw_noise_out} {output.netw_sec_intx_noise_out} {params.noise_perc} {params.n_rep}
        """

rule run_shuffle_edges:
    input:
        expr="input_data/rna_expr_shift.rds",
        meta="input_data/rna_meta.rds",
        network="decoupler_workflow/shuffle_edges/netw_shuffle_edges_frac{noise_perc}_rep{n_rep}.rds",
        sec_intx_network="decoupler_workflow/shuffle_edges/netw_shuffle_edges_frac{noise_perc}_rep{n_rep}_secondary_intx.pkl"
    params:
        component = lambda wildcards: "{}".format(wildcards.component),
        noise_perc = lambda wildcards: "{}".format(wildcards.noise_perc),
        n_rep = lambda wildcards: "{}".format(wildcards.n_rep)
    output:
        results_out="decoupler_workflow/shuffle_edges/{component}/shuffle_edges_frac{noise_perc}_rep{n_rep}.rds",
        res_file="decoupler_workflow/shuffle_edges/{component}/out_files/shuffle_edges_frac{noise_perc}_rep{n_rep}_results.tsv",
        auc_file="decoupler_workflow/shuffle_edges/{component}/out_files/shuffle_edges_frac{noise_perc}_rep{n_rep}_summary_results.tsv"
    singularity:
        "library://yasharw/priori/decoupler_env_no_delta_force_regulon_network"
    shell:
        """
        mkdir -p decoupler_workflow/shuffle_edges/{params.component}
        mkdir -p decoupler_workflow/shuffle_edges/{params.component}/out_files
        Rscript ./scripts/run_shuffle_edges.R {input.expr} {input.meta} {input.network} {input.sec_intx_network} {output.results_out} {output.res_file} {output.auc_file} {params.noise_perc} {params.n_rep} {params.component}
        """

rule prune_edges:
    input:
        network="input_data/pathway_commons_decoupler.rds"
    params:
        noise_perc = lambda wildcards: "{}".format(wildcards.noise_perc),
        n_rep = lambda wildcards: "{}".format(wildcards.n_rep)
    output:
        netw_noise_out="decoupler_workflow/prune_edges/netw_prune_edges_frac{noise_perc}_rep{n_rep}.rds",
        netw_sec_intx_noise_out="decoupler_workflow/prune_edges/netw_prune_edges_frac{noise_perc}_rep{n_rep}_secondary_intx.pkl"
    # singularity:
    #     "library://yasharw/priori/decoupler_env_no_delta_force_regulon_network"
    shell:
        """
        mkdir -p decoupler_workflow/prune_edges
        Rscript ./scripts/prune_edges.R {input.network} {output.netw_noise_out} {output.netw_sec_intx_noise_out} {params.noise_perc} {params.n_rep}
        """

rule run_prune_edges:
    input:
        expr="input_data/rna_expr_shift.rds",
        meta="input_data/rna_meta.rds",
        network="decoupler_workflow/prune_edges/netw_prune_edges_frac{noise_perc}_rep{n_rep}.rds",
        sec_intx_network="decoupler_workflow/prune_edges/netw_prune_edges_frac{noise_perc}_rep{n_rep}_secondary_intx.pkl"
    params:
        component = lambda wildcards: "{}".format(wildcards.component),
        noise_perc = lambda wildcards: "{}".format(wildcards.noise_perc),
        n_rep = lambda wildcards: "{}".format(wildcards.n_rep)
    output:
        results_out="decoupler_workflow/prune_edges/{component}/prune_edges_frac{noise_perc}_rep{n_rep}.rds",
        res_file="decoupler_workflow/prune_edges/{component}/out_files/prune_edges_frac{noise_perc}_rep{n_rep}_results.tsv",
        auc_file="decoupler_workflow/prune_edges/{component}/out_files/prune_edges_frac{noise_perc}_rep{n_rep}_summary_results.tsv"
    singularity:
        "library://yasharw/priori/decoupler_env_no_delta_force_regulon_network"
    shell:
        """
        mkdir -p decoupler_workflow/prune_edges/{params.component}
        mkdir -p decoupler_workflow/prune_edges/{params.component}/out_files
        Rscript ./scripts/run_prune_edges.R {input.expr} {input.meta} {input.network} {input.sec_intx_network} {output.results_out} {output.res_file} {output.auc_file} {params.noise_perc} {params.n_rep} {params.component}
        """

rule expand_edges:
    input:
        expr="input_data/rna_expr_shift.rds",
        network="input_data/pathway_commons_decoupler.rds"
    params:
        expand_perc = lambda wildcards: "{}".format(wildcards.expand_perc),
        n_rep = lambda wildcards: "{}".format(wildcards.n_rep)
    output:
        netw_noise_out="decoupler_workflow/expand_edges/netw_expand_edges_frac{expand_perc}_rep{n_rep}.rds",
        netw_sec_intx_noise_out="decoupler_workflow/expand_edges/netw_expand_edges_frac{expand_perc}_rep{n_rep}_secondary_intx.pkl"
    # singularity:
    #     "library://yasharw/priori/decoupler_env_no_delta_force_regulon_network"
    shell:
        """
        mkdir -p decoupler_workflow/expand_edges
        Rscript ./scripts/expand_edges.R {input.expr} {input.network} {output.netw_noise_out} {output.netw_sec_intx_noise_out} {params.expand_perc} {params.n_rep}
        """

rule run_expand_edges:
    input:
        expr="input_data/rna_expr_shift.rds",
        meta="input_data/rna_meta.rds",
        network="decoupler_workflow/expand_edges/netw_expand_edges_frac{expand_perc}_rep{n_rep}.rds",
        sec_intx_network="decoupler_workflow/expand_edges/netw_expand_edges_frac{expand_perc}_rep{n_rep}_secondary_intx.pkl"
    params:
        component = lambda wildcards: "{}".format(wildcards.component),
        expand_perc = lambda wildcards: "{}".format(wildcards.expand_perc),
        n_rep = lambda wildcards: "{}".format(wildcards.n_rep)
    output:
        results_out="decoupler_workflow/expand_edges/{component}/expand_edges_frac{expand_perc}_rep{n_rep}.rds",
        res_file="decoupler_workflow/expand_edges/{component}/out_files/expand_edges_frac{expand_perc}_rep{n_rep}_results.tsv",
        auc_file="decoupler_workflow/expand_edges/{component}/out_files/expand_edges_frac{expand_perc}_rep{n_rep}_summary_results.tsv"
    singularity:
        "library://yasharw/priori/decoupler_env_no_delta_force_regulon_network"
    shell:
        """
        mkdir -p decoupler_workflow/expand_edges/{params.component}
        mkdir -p decoupler_workflow/expand_edges/{params.component}/out_files
        Rscript ./scripts/run_expand_edges.R {input.expr} {input.meta} {input.network} {input.sec_intx_network} {output.results_out} {output.res_file} {output.auc_file} {params.expand_perc} {params.n_rep} {params.component}
        """