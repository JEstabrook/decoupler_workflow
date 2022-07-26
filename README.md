# decouplerRBench snakemake workflow
`decoupleRBench snakemake workflow` evaluates the performance of biological activity inference methods using perturbation experiments. This workflow incorporates the regulon-enrichment package enricher. It builds on `decoupleR`, and more specifically the `decouple` wrapper function. As such, it requires the decoupleR package to be installed and it is recommended that the user is familiar with the [basics of decoupleR](https://saezlab.github.io/decoupleR/articles/decoupleR.html#basics-1). The benchmark pipeline requires an input tibble with user-specified settings, benchmark data in the form of a count table and a corresponding metadata table.

## Set-up
To set-up necessary folders and environments to run the workflow, please run:

### Clone decoupler workflow page
Clone the contents of this Github onto your system:
```
git init
git clone https://github.com/JEstabrook/decoupler_workflow.git
cd decoupler_workflow/
```


### Create logs folder
Create a folder to store the logs files from the analysis:
```
mkdir -p logs/
```
     
### Prepare regulon-enrichment environment
```
conda create -n enrich_env
conda activate enrich_env 
conda install -c estabroj89 regulon-enrichment
```

### Snakemake
Ensure that snakemake is installed on your system and that it is at least version 5.0.0. To install snakemake in your environment, run:
```
conda install -c conda-forge -c bioconda snakemake
```

### Initialize Singularity
The environments for the workflow are containerized in Singularity. If Singularity is not installed on your system, please follow the steps in the [Singularity documentation](https://docs.sylabs.io/guides/3.0/user-guide/installation.html#). Note that Singularity may be installed on the compute nodes in your system and, therefore, can only be accessed in an interactive session. To initialize Singularity, please run: 
```
module load /path/to/singularity
```

## Run workflow
Use the following command to run the workflow:
```
sbatch submit_snakemake.sh
```

## Evaluation
For a given `decoupleR` method, activities are inferred for each regulator and 
experiment. To evaluate their performance, all experiments are concatenated 
together to generate a response vector (whether a regulator is perturbed or not)
and a predictor vector (the regulator activities). Then, using different 
thresholds we can calculate AUROC and AUPRC for each method. Given that the true 
positive classes are limited by the regulators covered in the perturbation 
experiments, we use a downsampling strategy, where for each permutation an 
equal number of negative classes are randomly sampled.
