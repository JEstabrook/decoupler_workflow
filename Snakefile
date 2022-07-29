import numpy as np
import pandas as pd 
import json
import os
from itertools import compress
from enricher.enrich import *

meta = pd.read_table('./input_data/KnockTF_GEO_w_ko_expanded_meta.tsv',index_col=0)
meta['cell'] = meta['cell'].str.replace('(','')
meta['cell'] = meta['cell'].str.replace(')','')
meta['cell'] = meta['cell'].str.replace('/','_')
#meta = meta[meta.cell.isin(meta['cell'].value_counts()[meta['cell'].value_counts() >24].index.tolist())]
meta = meta[meta['pct_ko'] > .25]
subset_kd = ['siRNA','shRNA','Drug']
meta = meta[meta.knockdown_method.isin(subset_kd)]


TARGET = meta.target.tolist()
CELL_LINE = meta.cell.tolist()
KD = meta.knockdown_method.tolist()
COMPONENTS = meta.components.tolist()

regulators = ['LEF1', 'PITX2', 'AR', 'ESRRA', 'ESR2', 'JUN', 'AHR', 'ZEB1', 'NR3C1', 'MYC', 'MYB', 'NFE2L2', 'ESR1', 'TFAP4', 'MSX1', 'GATA6', 'IRF7', 'STAT3', 'STAT5A', 'HIF1A', 'ELK1', 'ETS2', 'FOXA1', 'ZIC2', 'TP53', 'RELA', 'YY1', 'GTF3A', 'ATF2', 'MTF1', 'TCF4', 'EP300', 'CREB1', 'ATF3', 'SMAD3', 'SMAD4', 'FOXO1', 'RUNX1', 'ETS1', 'FOXO3', 'FOXM1', 'HSF1', 'PAX2', 'SOX9', 'TFAP2C', 'STAT6', 'HMGA1', 'RUNX2', 'PAX3', 'STAT1', 'MYCN', 'NR2F2', 'NFKB1', 'GATA3', 'HSF2', 'SF1', 'BRCA1', 'EGR3', 'EPAS1', 'CDX2', 'ETV4', 'BACH1', 'MAF', 'SMAD2', 'TCF7L1', 'IRF2', 'TP63', 'ARNT', 'ATG5', 'PGR', 'REST', 'XBP1', 'PTEN', 'SP3', 'KLF4', 'GATA2', 'HOXA5', 'ATM', 'ATR', 'NFATC3', 'HES1', 'POU5F1', 'NANOG', 'SOX2', 'TYK2', 'HNF1B', 'BRD4', 'NFYC', 'NFYA', 'SNAI1', 'GLI2', 'CHEK1', 'RELB', 'MITF', 'POSTN', 'PROX1', 'EREG', 'THRA', 'THRB', 'TWIST1', 'FOXP1', 'MYBL2', 'NFKB2', 'ACTL6A']
idx = [True if x in regulators else False for x in TARGET]


TARGET = list(compress(TARGET,idx))
CELL_LINE = list(compress(CELL_LINE,idx))
KD = list(compress(KD,idx))
COMPONENTS = list(compress(COMPONENTS,idx))


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
        expand("decoupler_workflow/mod_results/{component}/{kd}/{cell}/decoupler_subset_results.rds", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        #expand("decoupler_workflow/w_kd_plots/{component}/{kd}_{cell}_supplemental_figure2.pdf", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        #expand("decoupler_workflow/w_kd_plots_mod/{component}/{kd}_{cell}_supplemental_figure2.pdf", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        #expand("decoupler_workflow/w_kd_plots/{component}/{kd}_{cell}_supplemental_figure3.pdf", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        #expand("decoupler_workflow/w_kd_plots_mod/{component}/{kd}_{cell}_supplemental_figure3.pdf", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        expand("decoupler_workflow/kd_agnostic_results/{component}/{cell}/decoupler_subset_results.rds", cell=set(CELL_LINE), component=COMPONENTS),
        expand("decoupler_workflow/mod_kd_agnostic_results/{component}/{cell}/decoupler_subset_results.rds", cell=set(CELL_LINE), component=COMPONENTS),
        #expand("decoupler_workflow/wo_kd_plots/{component}/{cell}_supplemental_figure2.pdf",cell=set(CELL_LINE), component=COMPONENTS),
        #expand("decoupler_workflow/wo_kd_plots_mod/{component}/{cell}_supplemental_figure2.pdf",cell=set(CELL_LINE), component=COMPONENTS),
        #expand("decoupler_workflow/wo_kd_plots/{component}/{cell}_supplemental_figure3.pdf",cell=set(CELL_LINE), component=COMPONENTS),
        #expand("decoupler_workflow/wo_kd_plots_mod/{component}/{cell}_supplemental_figure3.pdf",cell=set(CELL_LINE), component=COMPONENTS),
        #expand("decoupler_workflow/w_kd_plots_sign_mod/{component}/{kd}_{cell}_supplemental_figure2.pdf", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        #expand("decoupler_workflow/w_kd_plots_sign_mod/{component}/{kd}_{cell}_supplemental_figure3.pdf", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        #expand("decoupler_workflow/wo_kd_plots_sign_mod/{component}/{cell}_supplemental_figure2.pdf",cell=set(CELL_LINE), component=COMPONENTS),
        #expand("decoupler_workflow/wo_kd_plots_sign_mod/{component}/{cell}_supplemental_figure3.pdf",cell=set(CELL_LINE), component=COMPONENTS),
        expand("decoupler_workflow/out_files/{component}_{cell}_wo_kd_results.tsv",cell=set(CELL_LINE), component=COMPONENTS),
        expand("decoupler_workflow/out_files/{component}_{kd}_{cell}_w_kd_results.tsv", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        expand("decoupler_workflow/mod_out_files/{component}_{cell}_wo_kd_results.tsv",cell=set(CELL_LINE), component=COMPONENTS),
        expand("decoupler_workflow/mod_out_files/{component}_{kd}_{cell}_w_kd_results.tsv", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        expand("decoupler_workflow/sign_mod_out_files/{component}_{cell}_wo_kd_results.tsv",cell=set(CELL_LINE), component=COMPONENTS),
        expand("decoupler_workflow/sign_mod_out_files/{component}_{kd}_{cell}_w_kd_results.tsv", zip, kd=KD,cell=CELL_LINE, component=COMPONENTS),
        expand("decoupler_workflow/sign_mod_kd_agnostic_results/{component}/{cell}/decoupler_priori_weights.tsv", zip, cell=CELL_LINE, component=COMPONENTS)

include: "rules/validation.smk"
