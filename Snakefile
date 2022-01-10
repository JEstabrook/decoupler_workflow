import numpy as np
import pandas as pd 
import json
import os
from itertools import compress
from enricher.enrich import *

meta = pd.read_table('./input_data/decoupler_L1000_meta.tsv').dropna()

TARGET = meta.target.tolist()
CELL_LINE = meta.cell_id.tolist()
KD = meta.pert_id.tolist()
TIME = meta.pert_time.tolist()
DOSE = meta.pert_dose.tolist()

regulators = ["ERBB2","MAP2K1","EGFR","BCR","ERBB4","CDK1","PIK3CA","MET","CDK9","BRAF","SRC","PLK1","CDK6","AKT1","PTK2","YES1","MTOR","PRKCB","PIM1","TEK","MAPK12","PARP1","HDAC6","HDAC1","HDAC4","PTGS2","DNMT1","MAPK14","TGFB1","ALK","IKBKB","TOP2A","MAPK7","GSK3B","MAP2K5","AKT2","KIT","CSF1R","MAPK3","CDK7","ATM","SIRT1","HDAC3","EHMT2","HDAC2","ATR","PRKDC","CHEK2","RPS6KA3","RAF1","MAPK11","MAPK8","JAK2","MDM2","JAK3"]
idx = [True if x in regulators else False for x in TARGET]

TARGET = list(compress(TARGET,idx))
CELL_LINE = list(compress(CELL_LINE,idx))
KD = list(compress(KD,idx))
TIME = list(compress(TIME,idx))
DOSE = list(compress(DOSE,idx))
DOSE = np.round(DOSE,2)
print(DOSE)

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
        expand("decoupler_workflow/results/{kd}/{time}/{cell}/decoupler_subset_results.rds", zip, kd=KD,time=TIME,cell=CELL_LINE),
        expand("decoupler_workflow/w_kd_plots/{kd}_{time}_{cell}_supplemental_figure2.pdf", zip, kd=KD,time=TIME,cell=CELL_LINE),
        expand("decoupler_workflow/w_kd_plots/{kd}_{time}_{cell}_supplemental_figure3.pdf", zip, kd=KD,time=TIME,cell=CELL_LINE),
        expand("decoupler_workflow/kd_agnostic_results/{time}/{cell}/decoupler_subset_results.rds", zip, time=TIME,cell=CELL_LINE),
        expand("decoupler_workflow/wo_kd_plots/{time}_{cell}_supplemental_figure2.pdf", zip, time=TIME,cell=CELL_LINE),
        expand("decoupler_workflow/wo_kd_plots/{time}_{cell}_supplemental_figure3.pdf", zip, time=TIME,cell=CELL_LINE),
        expand("decoupler_workflow/results_dose/{kd}_{dose}/{time}/{cell}/decoupler_subset_results.rds", zip, kd=KD,dose=DOSE,time=TIME,cell=CELL_LINE),
        expand("decoupler_workflow/w_kd_dose_plots/{kd}_{dose}_{time}_{cell}_supplemental_figure2.pdf", zip, kd=KD,dose=DOSE,time=TIME,cell=CELL_LINE),
        expand("decoupler_workflow/w_kd_dose_plots/{kd}_{dose}_{time}_{cell}_supplemental_figure3.pdf", zip, kd=KD,dose=DOSE,time=TIME,cell=CELL_LINE)

include: "rules/validation.smk"
