import numpy as np
import pandas as pd 
import json
import os
from itertools import compress

COMPONENTS = ['local_enrichment', 'quant_nes', 'total_enrichment']
# COMPONENTS = ['total_enrichment']

NOISE_PERC = [0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1.0] 

EXPAND_PERC = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 4.0, 9.0] 

# GAUSS_SD = [1, 2, 3, 4, 5]
GAUSS_SD = [1, 5, 10, 20, 30, 40, 50, 100]

# N_REP = [1, 2, 3, 4, 5] 
N_REP = [1] 

with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())
for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)

def message(mes):
    sys.stderr.write("|--- " + mes + "\n")


rule all:
    input:
        # Main pipeline
        # expand(["decoupler_workflow/results/{component}/decoupler_subset_results.rds", 
        #     "decoupler_workflow/out_files/{component}_results.tsv", 
        #     "decoupler_workflow/results/{component}/decoupler_priori_weights.tsv",
        #     ], component=COMPONENTS),
        # # Gaussian noise
        # expand(["decoupler_workflow/rna_gaussian_noise/expr_rna_noise_sd{gauss_sd}_rep{n_rep}.rds",
        #     "decoupler_workflow/rna_gaussian_noise/rna_gaussian_noise_sd{gauss_sd}_rep{n_rep}.rds", 
        #     "decoupler_workflow/rna_gaussian_noise/{component}/out_files/rna_gaussian_noise_sd{gauss_sd}_rep{n_rep}_results.tsv",
        #     "decoupler_workflow/rna_gaussian_noise/{component}/out_files/rna_gaussian_noise_sd{gauss_sd}_rep{n_rep}_summary_results.tsv",
        #     ], component=COMPONENTS, gauss_sd=GAUSS_SD, n_rep=N_REP),
        # Random edges
        expand(["decoupler_workflow/random_edges/netw_random_edges_frac{noise_perc}_rep{n_rep}.rds",
            "decoupler_workflow/random_edges/netw_random_edges_frac{noise_perc}_rep{n_rep}_secondary_intx.pkl",
            "decoupler_workflow/random_edges/{component}/random_edges_frac{noise_perc}_rep{n_rep}.rds", 
            "decoupler_workflow/random_edges/{component}/out_files/random_edges_frac{noise_perc}_rep{n_rep}_results.tsv",
            "decoupler_workflow/random_edges/{component}/out_files/random_edges_frac{noise_perc}_rep{n_rep}_summary_results.tsv",
            ], component=COMPONENTS, noise_perc=NOISE_PERC, n_rep=N_REP),
        # Shuffle edges
        expand(["decoupler_workflow/shuffle_edges/netw_shuffle_edges_frac{noise_perc}_rep{n_rep}.rds",
            "decoupler_workflow/shuffle_edges/netw_shuffle_edges_frac{noise_perc}_rep{n_rep}_secondary_intx.pkl",
            "decoupler_workflow/shuffle_edges/{component}/shuffle_edges_frac{noise_perc}_rep{n_rep}.rds",
            "decoupler_workflow/shuffle_edges/{component}/out_files/shuffle_edges_frac{noise_perc}_rep{n_rep}_results.tsv",
            "decoupler_workflow/shuffle_edges/{component}/out_files/shuffle_edges_frac{noise_perc}_rep{n_rep}_summary_results.tsv",
            ], component=COMPONENTS, noise_perc=NOISE_PERC, n_rep=N_REP),
        # Prune edges
        expand(["decoupler_workflow/prune_edges/netw_prune_edges_frac{noise_perc}_rep{n_rep}.rds",
            "decoupler_workflow/prune_edges/netw_prune_edges_frac{noise_perc}_rep{n_rep}_secondary_intx.pkl",
            "decoupler_workflow/prune_edges/{component}/prune_edges_frac{noise_perc}_rep{n_rep}.rds",
            "decoupler_workflow/prune_edges/{component}/out_files/prune_edges_frac{noise_perc}_rep{n_rep}_results.tsv",
            "decoupler_workflow/prune_edges/{component}/out_files/prune_edges_frac{noise_perc}_rep{n_rep}_summary_results.tsv",
            ], component=COMPONENTS, noise_perc=NOISE_PERC, n_rep=N_REP),
        # Expand edges
        expand(["decoupler_workflow/expand_edges/netw_expand_edges_frac{expand_perc}_rep{n_rep}.rds",
            "decoupler_workflow/expand_edges/netw_expand_edges_frac{expand_perc}_rep{n_rep}_secondary_intx.pkl",
            "decoupler_workflow/expand_edges/{component}/expand_edges_frac{expand_perc}_rep{n_rep}.rds", 
            "decoupler_workflow/expand_edges/{component}/out_files/expand_edges_frac{expand_perc}_rep{n_rep}_results.tsv",
            "decoupler_workflow/expand_edges/{component}/out_files/expand_edges_frac{expand_perc}_rep{n_rep}_summary_results.tsv",
            ], component=COMPONENTS, expand_perc=EXPAND_PERC, n_rep=N_REP)

include: "rules/decoupler_test_data.smk"
