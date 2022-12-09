import numpy as np
import pandas as pd 
import json
import os
from itertools import compress

COMPONENTS = ['local_enrichment', 'quant_nes', 'total_enrichment']

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
        expand(["decoupler_workflow/results/{component}/decoupler_subset_results.rds", 
            "decoupler_workflow/out_files/{component}_results.tsv", 
            "decoupler_workflow/results/{component}/decoupler_priori_weights.tsv",
            ], component=COMPONENTS)

include: "rules/decoupler_test_data.smk"
