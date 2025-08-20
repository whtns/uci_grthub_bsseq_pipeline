# RNA-seq Analysis Snakemake Workflow

This Snakemake workflow has been converted from the original SLURM batch script `test.sub` and performs RNA-seq analysis including:

1. **Trimming** with Trimmomatic
2. **Alignment** with HISAT2
3. **Sorting and indexing** BAM files with samtools
4. **Feature counting** with featureCounts
5. **Quantification** with Salmon

## Files Created

- `Snakefile` - Main workflow definition
- `config.yaml` - Configuration file with paths and parameters
- `cluster.yaml` - SLURM cluster configuration
- `submit_snakemake.sh` - Script to submit the workflow to SLURM

## Usage

### 1. Prerequisites

Make sure you have Snakemake installed. You can install it using conda:

```bash
conda install -c conda-forge -c bioconda snakemake
```

### 2. Configuration

Edit `config.yaml` to modify:
- Sample names
- Input/output paths
- Reference file locations
- Tool parameters

### 3. Running the workflow

#### Option A: Submit to SLURM cluster
```bash
sbatch submit_snakemake.sh
```

#### Option B: Run locally (for testing)
```bash
snakemake --cores 8 --use-conda
```

#### Option C: Dry run (to check workflow)
```bash
snakemake --dry-run
```

### 4. Workflow visualization

Generate a workflow diagram:
```bash
snakemake --dag | dot -Tpng > workflow.png
```

## Key Differences from Original Script

1. **Modular design**: Each step is a separate rule
2. **Dependency management**: Snakemake automatically handles job dependencies
3. **Parallel execution**: Multiple samples can be processed simultaneously
4. **Configuration-driven**: Easy to modify parameters without editing the main workflow
5. **Resource management**: Better integration with SLURM scheduler
6. **Reproducibility**: Workflow tracks input/output dependencies

## Customization

### For different library types:
- **Non-stranded libraries**: Change `rna_strandness` to "unstranded" and `library_type` to "IU" in `config.yaml`
- **Different strand orientation**: Modify the strandness parameters accordingly

### For multiple samples:
Modify the workflow to accept a samples list and use wildcards for processing multiple samples in parallel.

## Output Files

The workflow generates:
- Trimmed FASTQ files (`*_trimmed_1P.fq.gz`, `*_trimmed_2P.fq.gz`)
- Aligned BAM files (`*_align_sorted.bam`)
- Feature count tables (`*_counts.txt`)
- Salmon quantification files (`*_quant.sf`)

## Troubleshooting

1. Check SLURM job status: `squeue -u $USER`
2. View workflow status: `snakemake --summary`
3. Check individual rule logs in the SLURM output files


# TODO
1. retrieve counts from featureCounts and Salmon quantification files, and summarize them in a final report
2. Run DESeq2 or edgeR for differential expression analysis
3. 