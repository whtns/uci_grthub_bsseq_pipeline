# WGBS Snakemake Workflow

This workflow is designed to run the basic steps for a whole-genome bisulfite 
sequencing experiment. It's intended to automate the workflow for future-use and
reproducibility. Its design is explicitly simple to make it easy for users to not
only understand the order and purpose of each step, but to be able to look at the
code, figure out how it works and get it running extremely easily.

One advantage of running this through Snakemake is that it intelligently handles
threading and replaces completed processes up to the number of cores specified
at run-time. However, options for the thread count for each step are configurable
in the .yaml file.

## Dependencies
Most recent tested versions indicated. Though, more recent versions and slightly
older ones are likely a-okay!

1. [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) v0.6.4
2. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v0.11.8
3. [bwa-meth](https://github.com/brentp/bwa-meth) v0.2.2
4. [samtools](https://www.htslib.org/) v1.9
5. [Picard Tools](https://broadinstitute.github.io/picard/) v2.22.3
6. [MethylDackel](https://github.com/dpryan79/MethylDackel) v0.4
7. [Mosdepth](https://github.com/brentp/mosdepth) v0.2.9
8. [Snakemake](https://snakemake.readthedocs.io) v5.14.0
9. [Python3](https://www.python.org/) v3.5

## Quick-Start Guide
Firstly, download all the dependencies and make sure they're in your $PATH (that
you can run them from a BASH prompt). Then clone the github repo:

```shell
git clone https://github.com/groverj3/wgbs_snakemake.git
```

Edit the config.yaml file to include your sample IDs (fastq filenames,
excluding extensions, pair numbers, lane info, etc.) and a reference genome
(which may be pre-indexed). You'll definitely want to make sure that the adapter
sequence in there matches what's in your samples.

Currently, the workflow expects an R1 and R2 file for each sample. Place the
individual .fastq.gz files for R1 and R2 into the input_data directory. Once
you have all the required dependencies installed run the workflow with:

```shell
snakemake --cores {cores_here}
```

## Configuration: output directory and FASTQ locations

This workflow now supports a top-level output directory that centralizes all
generated files. Set it in `config.yaml` with the key `output_dir`. If not set,
the workflow defaults to `output` and will create directories such as
`output/trimmed`, `output/1_fastqc_raw`, `output/3_aligned_sorted_markdupes`,
etc.

Example `config.yaml` fragment:

```yaml
output_dir: "output"
paths:
    fastqs: "data/FASTQ"  # optional; where to look for FASTQ files
# other keys...
```

Notes:
- The Snakefile will look for FASTQs in `paths.fastqs`, then `paths.data`, and
    finally a default `data/FASTQ` directory. You can change that by editing
    `config.yaml`.
- To inspect what Snakemake will do without running commands, use a dry-run:

```shell
snakemake -n -p --cores {cores_here}
```


## Workflow
1. Index the reference genome with bwameth and samtools faidx
2. Quality checking, and output of sample information with FastQC
3. Adapter and quality trimming with Trim Galore!
4. Alignment to a reference genome with bwa-meth
5. Marking PCR duplicates with Picard Tools MarkDuplicates
6. Detecting methylation bias per read position with MethylDackel
7. Extracting methylation calls per position into bedGraph and methylKit formats with MethylDackel
8. Calculating depth and coverage with Mosdepth

The output from the workflow is suitable for DMR-calling or aggregation of calls
to determine % methylation per feature.

![DAG](dag.png)

## All About the Workflow

Whole-genome bisulfite sequencing is a modification of whole-genome shotgun
sequencing designed to convert unmethylated cytosines into uracil. These uracils
are then sequenced as thymine. By tallying up the number of cytosines and
thymines for each cytosine in the reference genome you can then calculate a
percentage methylation for each annotated cytosine in your species' reference
assembly.

While it is possible to determine this by calling C -> T SNPs against a reference
there is purpose-built software for this task. The most common aligner for WGBS
is currently [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/),
and in our testing it performed well. However, we decided to use a slightly
different pipeline for the purposes of our work in the
[Mosher Lab](https://cals.arizona.edu/research/mosherlab/Mosher_Lab/Home.html).
The pipeline we settled on is a combination of open source tools built around
bwameth for alignment and MethylDackel for methylation calling. In our experience
this pipeline was many times faster and resulted in a marginally higher mapping
rate. Additionally, it uses Picard Tools to mark potential PCR duplicates, and
its method for doing so is not as conservative as Bismark's internal version of
the same process. Bismark's speed has improved in more recent versions, and is
under more active development but bwameth still produces comparable results in
less time.

Use of MethylDackel allows us to determine per-position biases in terms of
methylation calls on the reads, and different ones based on read orientation.
Using some shell script hacking we can extract its recommendations for inclusion
bounds for methylation calling based on these biases and use in the methylation
calling step. This should reduce false positive or negatives based on effects of
cytosines being too close to adapters, interference from end-repair, or simply
incomplete trimming.

At the conclusion of the pipeline overall fold-coverage is calculated using the
very fast Mosdepth tool from Brent Pedersen and some python scripts.

## Citing the Workflow
Please do cite us! The included Zenodo DOI is the easiest way. Additionally, you
should consider citing the paper in which we first used this workflow:

Grover JW *et al*. Abundant expression of maternal siRNAs is a conserved feature of seed development. 2020.
PNAS. https://doi.org/10.1073/pnas.2001332117
