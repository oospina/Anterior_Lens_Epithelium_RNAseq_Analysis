#!/bin/bash
#SBATCH --job-name=nfcore_rnaseq
#SBATCH --output=nfcore_rnaseq_run/slurm_output.out
#SBATCH --error=nfcore_rnaseq_run/slurm_error.err
#SBATCH --partition=cpu
#SBATCH --mem=256GB
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00

# Load necessary modules
module load Nextflow
module load sbgrid/x86_64_Linux
module load sbgrid/java/jdk1.8.0_144
module load jdk/21.0.6

echo "Slurm - Start analysis: "
echo $(date)

# Run pipeline
nextflow run nf-core/rnaseq \
    --input data/eye_rnaseq_nfcore_samplesheet.csv \
    --outdir nfcore_rnaseq_run/nfcore_rnaseq_output \
    --fasta data/GRCh38.p14_from_NCBI/GRCh38.p14_ncbi_dataset/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna \
    --gtf data/GRCh38.p14_from_NCBI/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz \
    --star_index data/genome_idx_nfcore_rnaseq/index/star \
    --featurecounts_group_type gene_biotype \
    --skip_pseudo_alignment \
    --remove_ribo_rna \
    --skip_bigwig \
    -profile singularity

echo "Slurm - End analysis: "
echo $(date)


