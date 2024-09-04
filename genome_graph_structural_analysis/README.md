<h1>Structural analysis of Genome Graph</h1>

This repo contains codes to construct genome graphs and analyse its structure in-depth. <br>
Kindly execute the steps linearly and conserve the directory format. <br>
*Tested on Linux/Unix machines*

<br>

**Creating and activating the conda environnment**

The tools required for performing genome graph related analyses can be installed using conda [Install Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html). <br>
To install all the tools in one go, create a conda environment using the yml file provied with this repo at ```envs/structure_gg.yaml```.

1) Create the conda virtual environment: 
```
conda env create -f structure_gg.yaml
```

To verify if above step ran smoothly, run ```conda info --envs``` and check if ```structure_gg``` gets listed. 

2) Activate the environment created: 
```
conda activate structure_gg
```

<br>
<br> 

**Analysis of WGS w.r.t. genome graphs**

The directory contains a portable Snakemake based workflow that can be used for the structural analysis of genome graphs.

1) Move the processed bgzipped and indexed VCF files to ```datasets/input_files/processed_vcf_files/```. The files should have the extension .vcf.gz and the correspnding indices should have the extension .vcf.gz.tbi

2) Edit the filepaths of the linear reference genome and the karyotype file in ```config/config.yaml```.

3) Move the karyotype file describing the chromosomes into ```datasets/input_files/```. The current directory contains the karyotype file for human reference genome GRCH38.

4) Activate the ```structure_gg``` conda environment.
```
conda activate structure_gg
```

5) Perform a Snakemake dry run to ensure that the workflow works.
```
snakemake -npr
``` 


6) Once the dry run has run successfully, the structural analysis can be started with the command below:
```
snakemake --cores [integer]--jobs [integer] -p -r 
```

> Note: Edit the parameters accordingly <br>
> --cores : Use at most N CPU cores/jobs in parallel
> --jobs : max number of jobs to be performed parallely. <br>
> --use-conda: to enable conda through job scheduler <br>
>
> *Optional parameters:* <br>
> -p & -r: print the command and the reason to run the command
<br> 
> --cluster: the qsub command that snakemake uses to submit individual jobs. <br>
> --forceall: Mandatorily run all the rules in the Snakefile <br>
> --rerun-incomplete: If the analysis stops before completion, this parameter can be added to recreate any incomplete files. <br>


Kindly look into [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/) for more details.