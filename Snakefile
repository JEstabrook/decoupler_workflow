import numpy as np
import pandas as pd 
import json
import os
from itertools import compress
from enricher.enrich import *

meta = pd.read_table('./input_data/Updated_knocktf_meta_with_components.tsv',index_col=0)

TARGET = meta.target.tolist()
CELL_LINE = meta.cell.tolist()
KD = meta.knockdown_method.tolist()
COMPONENTS = meta.components.tolist()

regulators = ['STAT2','NR2C2','GATA1','U2AF2','NFATC1','CITED2','RELA','NFE2L1','SETDB1','MXI1','BCLAF1','E2F6','MAFK','MITF','SP2','STAT5A','TFDP1','MAFG','RCOR1','SNW1','NRF1','NR4A1','HSF1','USF1','CEBPZ','ERF','CTCF','MAX','NFYB','HMGN3','E2F4','HDAC8','SP1','BHLHE40','ZBTB33','BRCA1','TRIM28','ZNF143','HMGA1','USF2','GTF2F1','JUND','SMAD5','TAL1','KAT2B','TBL1XR1','STAT6','ATF3','MAZ','SRF','LMNA','PCBP1','NFE2L2','RFX5','SIX5','NR2F2','GATA2','RAD21','FOXM1','CTBP1','STAT1','POLR2G','CHD2','SMARCA4']
idx = [True if x in regulators else False for x in TARGET]


TARGET = list(compress(TARGET,idx))
CELL_LINE = list(compress(CELL_LINE,idx))
KD = list(compress(KD,idx))
COMPONENTS = list(compress(COMPONENTS,idx))

print(KD)
print(CELL_LINE)
print(COMPONENTS)

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
        expand("decoupler_workflow/results/{component}/{kd}/{cell}/decoupler_subset_results.rds", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        expand("decoupler_workflow/w_kd_plots/{component}/{kd}_{cell}_supplemental_figure2.pdf", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        expand("decoupler_workflow/w_kd_plots/{component}/{kd}_{cell}_supplemental_figure3.pdf", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        expand("decoupler_workflow/kd_agnostic_results/{component}/{cell}/decoupler_subset_results.rds", cell=set(CELL_LINE), component=COMPONENTS),
        expand("decoupler_workflow/wo_kd_plots/{component}/{cell}_supplemental_figure2.pdf",cell=set(CELL_LINE), component=COMPONENTS),
        expand("decoupler_workflow/wo_kd_plots/{component}/{cell}_supplemental_figure3.pdf",cell=set(CELL_LINE), component=COMPONENTS),
        expand("decoupler_workflow/out_files/{component}_{cell}_wo_kd_results.tsv",cell=set(CELL_LINE), component=COMPONENTS),
        expand("decoupler_workflow/out_files/{component}_{kd}_{cell}_w_kd_results.tsv", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        expand("decoupler_workflow/sign_mod_out_files/{component}_{cell}_wo_kd_results.tsv",cell=set(CELL_LINE), component=COMPONENTS),
        expand("decoupler_workflow/sign_mod_out_files/{component}_{kd}_{cell}_w_kd_results.tsv", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        expand("decoupler_workflow/w_kd_plots_sign_mod/{component}/{kd}_{cell}_supplemental_figure2.pdf", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        expand("decoupler_workflow/w_kd_plots_sign_mod/{component}/{kd}_{cell}_supplemental_figure3.pdf", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        expand("decoupler_workflow/wo_kd_plots_sign_mod/{component}/{cell}_supplemental_figure2.pdf",cell=set(CELL_LINE), component=COMPONENTS),
        expand("decoupler_workflow/wo_kd_plots_sign_mod/{component}/{cell}_supplemental_figure3.pdf",cell=set(CELL_LINE), component=COMPONENTS),

include: "rules/validation.smk"
