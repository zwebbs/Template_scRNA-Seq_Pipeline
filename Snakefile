# File Name: Zach Weber
# Created By: ZW
# Created On: 2023-06-29
# Purpose: runs the cellranger scRNAseq pipeline
#  and downstream analysis in Seurat for the 
#  EHT single cell datasets.

# Module Imports
# -----------------------------------------------------------------------------
from gardnersnake import Configuration
from pathlib import Path


# Global Configuration
# -----------------------------------------------------------------------------
yaml_config_filepath = Path(config["yaml_config"])
cfg = Configuration(filepath=yaml_config_filepath)
cfg.load()

# define global configurations
GLOBALS = cfg.global_params
WORKDIR = GLOBALS.working_directory
workdir: WORKDIR

# define directories and files objects
DIRECTORIES = GLOBALS.misc.directories
FILES = GLOBALS.files

# define specific directories that are used often
SCRATCH_DIR = Path(DIRECTORIES.scratch_directory)
FASTQ_DIR = Path(DIRECTORIES.fastq_directory)
ALIGN_DIR = Path(DIRECTORIES.align_directory)
DATA_DIR = Path(DIRECTORIES.data_directory)
SCRIPTS_DIR = Path(DIRECTORIES.scripts_directory)
PLOT_DIR = Path(DIRECTORIES.plot_directory)
REF_DIR = Path(DIRECTORIES.reference_directory)
LOGS_DIR = Path(DIRECTORIES.logs_directory)

# get run_ids from sequencing data
SEQ = FILES.sequencing
RUN_IDS = [s["run_id"] for s in SEQ]
print(RUN_IDS)

# RULE 0. Global Outputs
# -----------------------------------------------------------------------------
rule All:
    input:
        expand(ALIGN_DIR / "{run_id}.cellranger.rc.out", run_id=RUN_IDS),
        f"{GLOBALS.analysis_name}.integrated.RData"

# RULE 1. Align Data using 10X CellRanger Count
# ----------------------------------------------------------------------------- 
cellranger_count_rp = cfg.get_rule_params("CellRanger_Count")
rule CellRanger_Count:
    input:
        manifest = (FASTQ_DIR / "{run_id}/manifest.txt")
    params: **(cellranger_count_rp.parameters),
        fastq_dir = lambda wildcards: (FASTQ_DIR / f"{wildcards.run_id}/"),
        align_dir = ALIGN_DIR,
        reference_dir = REF_DIR 
    resources: **(cellranger_count_rp.resources),
        job_id = lambda wildcards: f"{wildcards.run_id}",
        logs = str(LOGS_DIR)
    output:
        cellranger_rc = (ALIGN_DIR / "{run_id}.cellranger.rc.out")
    shell:
        "cd {params.align_dir} && "
        "cellranger count"
        " --id {resources.job_id}"
        " --transcriptome {params.reference_dir}"
        " --fastqs {params.fastq_dir}"
        " --sample {resources.job_id}"
        " --localcores {params.localcores}"
        " --localmem {params.localmemgb} && "
        "check_directory -o {output.cellranger_rc}"
        " {params.checkfiles} {resources.job_id}/outs/"


# RULE 2. Detect Doublets using scrublet
# ------------------------------------------------------------------------------
scrublet_detection_rp = cfg.get_rule_params("Scrublet_Doublet_Detect")
rule Scrublet_Doublet_Detect:
    input: 
        cellranger_rc = rules.CellRanger_Count.output.cellranger_rc
    params: **(scrublet_detection_rp.parameters),
        data_dir = lambda wildcards: (ALIGN_DIR / f"{wildcards.run_id}/outs/filtered_feature_bc_matrix/"),
        plot_dir = PLOT_DIR,
        out_dir = (DATA_DIR / "scrublet/"),
        sample_name = lambda wildcards: f"{wildcards.run_id}"
    resources: **(scrublet_detection_rp.resources),
        job_id = lambda wildcards: f"{wildcards.run_id}",
        logs = str(LOGS_DIR)
    output:
        (PLOT_DIR / "{run_id}_scrublet_hist.png"),
        (PLOT_DIR / "{run_id}_scrublet_umap.png"),
        (DATA_DIR / "scrublet/{run_id}_scrublet_doublets.csv")
    shell:
        "mkdir -p {params.plot_dir} {params.out_dir} && "
        "python scripts/doublet_detection.py -r {params.doublet_prior_rate}"
        " --datadir {params.data_dir} --sample-name {params.sample_name}"
        " --outdir {params.out_dir} --plotdir {params.plot_dir}"


# RULE 2. Integrate Cellranger Outputs using Seurat
# ------------------------------------------------------------------------------
integrate_datasets_rp = cfg.get_rule_params("Integrate_Datasets")
rule Integrate_Datasets:
    input:
        scr_files = expand(DATA_DIR / "scrublet/{run_id}_scrublet_doublets.csv", run_id=RUN_IDS),
        rcs = expand(ALIGN_DIR / "{run_id}.cellranger.rc.out", run_id=RUN_IDS),
        manifest = FILES.integration_manifest
    params: **(integrate_datasets_rp.parameters),
        scripts_dir = SCRIPTS_DIR
    resources: **(integrate_datasets_rp.resources),
        job_id = GLOBALS.analysis_name,
        logs = str(LOGS_DIR)
    output:
        rdata_file = f"{GLOBALS.analysis_name}.integrated.RData"
    shell:
        "Rscript {params.scripts_dir}/Integrate_Data.R"
        " --manifest {input.manifest}"
        " --outfile {output.rdata_file}"

