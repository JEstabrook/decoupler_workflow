import numpy as np
import pandas as pd 
import json
import os
from itertools import compress,cycle
from enricher.enrich import *

meta = pd.read_table('./input_data/KnockTF_joined_meta.tsv',index_col=0)

TARGET = meta.target.tolist()
CELL_LINE = meta.cell.tolist()
KD = meta.knockdown_method.tolist()

regulators = ['STAT2','NR2C2','GATA1','U2AF2','NFATC1','CITED2','RELA','NFE2L1','SETDB1','MXI1','BCLAF1','E2F6','MAFK','MITF','SP2','STAT5A','TFDP1','MAFG','RCOR1','SNW1','NRF1','NR4A1','HSF1','USF1','CEBPZ','ERF','CTCF','MAX','NFYB','HMGN3','E2F4','HDAC8','SP1','BHLHE40','ZBTB33','BRCA1','TRIM28','ZNF143','HMGA1','USF2','GTF2F1','JUND','SMAD5','TAL1','KAT2B','TBL1XR1','STAT6','ATF3','MAZ','SRF','LMNA','PCBP1','NFE2L2','RFX5','SIX5','NR2F2','GATA2','RAD21','FOXM1','CTBP1','STAT1','POLR2G','CHD2','SMARCA4']
idx = [True if x in regulators else False for x in TARGET]


TARGET = list(compress(TARGET,idx))
CELL_LINE = list(compress(CELL_LINE,idx))
KD = list(compress(KD,idx))

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

CONTROL_COND = ['kdallcontrol', 'kdonlycontrol']
no_kd_control =  ['nokdallcontrol', 'nokdonlycontrol']
DOWNSAMPLE = ['TRUE', 'FALSE']
zip_list = zip(CELL_LINE, KD, cycle(DOWNSAMPLE),cycle(CONTROL_COND))
cell,kd,down,control = list(zip(*zip_list))

zip_list = zip(CELL_LINE, KD, cycle(DOWNSAMPLE),cycle(no_kd_control))
_,_,_,nokdcontrol = list(zip(*zip_list))


rule all:
    input:
        expand("decoupler_workflow/results/{kd}/{cell}/decoupler_subset_results.rds", zip, kd=KD,cell=CELL_LINE),
        expand("decoupler_workflow/w_kd_plots/{kd}_{cell}_supplemental_figure2.pdf", zip, kd=KD,cell=CELL_LINE),
        expand("decoupler_workflow/w_kd_plots/{kd}_{cell}_supplemental_figure3.pdf", zip, kd=KD,cell=CELL_LINE),
        expand("decoupler_workflow/kd_agnostic_results/{cell}/decoupler_subset_results.rds", cell=set(CELL_LINE)),
        expand("decoupler_workflow/wo_kd_plots/{cell}_supplemental_figure2.pdf",cell=set(CELL_LINE)),
        expand("decoupler_workflow/wo_kd_plots/{cell}_supplemental_figure3.pdf",cell=set(CELL_LINE)),
        expand("decoupler_workflow/downsample/kdallcontrol_{downsample}_{kd}_{cell}_decoupler_subset_results.rds", zip, kd=KD, cell=CELL_LINE, downsample=down),
        expand("decoupler_workflow/downsample/kdonlycontrol_{downsample}_{kd}_{cell}_decoupler_subset_results.rds", zip, kd=KD, cell=CELL_LINE, downsample=down),
        expand("decoupler_workflow/downsample/nokdallcontrol_{downsample}_{cell}_decoupler_subset_results.rds", zip, cell=CELL_LINE, downsample=down),
        expand("decoupler_workflow/downsample/nokdonlycontrol_{downsample}_{cell}_decoupler_subset_results.rds", zip, cell=CELL_LINE, downsample=down),
        expand("decoupler_workflow/downsample/w_kd_plots/kdallcontrol_{downsample}_{kd}_{cell}_supplemental_figure2.pdf",zip, kd=KD,cell=CELL_LINE, downsample=down),
        expand("decoupler_workflow/downsample/w_kd_plots/kdonlycontrol_{downsample}_{kd}_{cell}_supplemental_figure2.pdf",zip, kd=KD,cell=CELL_LINE, downsample=down),
        expand("decoupler_workflow/downsample/wo_kd_plots/nokdallcontrol_{downsample}_{cell}_supplemental_figure2.pdf",zip, cell=CELL_LINE, downsample=down),
        expand("decoupler_workflow/downsample/wo_kd_plots/nokdonlycontrol_{downsample}_{cell}_supplemental_figure2.pdf",zip, cell=CELL_LINE, downsample=down),
        expand("decoupler_workflow/downsample/w_kd_plots/kdallcontrol_TRUE_{kd}_{cell}_supplemental_figure3.pdf",zip, kd=KD,cell=CELL_LINE),
        expand("decoupler_workflow/downsample/w_kd_plots/kdonlycontrol_TRUE_{kd}_{cell}_supplemental_figure3.pdf",zip, kd=KD,cell=CELL_LINE),
        expand("decoupler_workflow/downsample/wo_kd_plots/nokdallcontrol_TRUE_{cell}_supplemental_figure3.pdf",zip, cell=CELL_LINE),
        expand("decoupler_workflow/downsample/wo_kd_plots/nokdonlycontrol_TRUE_{cell}_supplemental_figure3.pdf",zip, cell=CELL_LINE),
        expand("decoupler_workflow/wo_kd_plots/{cell}_base_plot_PR_HEAT.pdf",zip, cell=CELL_LINE),
        expand("decoupler_workflow/w_kd_plots/{kd}_{cell}_base_plot_PR_HEAT.pdf",zip, kd=KD, cell=CELL_LINE),
        expand("decoupler_workflow/nokddownsample/baseplots/{nokdcontrol}__{downsample}_{cell}_base_plot_PR_HEAT.pdf",zip, nokdcontrol=nokdcontrol, downsample=down, cell=cell),
        expand("decoupler_workflow/wkddownsample/baseplots/{control}_{kd}_{downsample}_{cell}_base_plot_PR_HEAT.pdf",zip, control=control, kd=KD, downsample=down, cell=cell)


include: "rules/validation.smk"
