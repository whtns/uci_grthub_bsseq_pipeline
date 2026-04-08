# Author: Jeffrey Grover
# Purpose: Run the whole-genome bisulfite sequencing workflow
# Created: 2019-05-22


# Get overall workflow parameters from config.yaml
configfile: 'config.yaml'

# SAMPLES = config['samples']
REFERENCE_GENOME = config['reference_genome']
# Base output directory (can be set in config.yaml as `output_dir: "output"`)
OUTPUT_DIR = config.get('output_dir', 'output')
import glob, os, re

FASTQ_DIR = 'data/FASTQ'

# Auto-detect samples from FASTQ directory
def get_samples_from_fastq_dir():
    """
    Detect sample names by inspecting the FASTQ directory structure.

    This project uses FASTQ filenames like
    "mR480-L1-P01-TAAGGCGA-CTCTCTAT-READ1-Sequences.fastq.gz".
    We therefore look for files containing "READ1" (or "READ2") and
    extract the sample prefix to the left of the "-READ" marker.

    Supported layouts:
    - per-sample subfolders: <FASTQ_DIR>/<sample>/*READ1*.fastq.gz
    - flat files in FASTQ_DIR: <FASTQ_DIR>/*READ1*.fastq.gz

    Returns a sorted list of unique sample names, excluding 'Undetermined'.
    """
    # Prefer configured keys; fall back to default
    paths_cfg = config.get("paths", {})
    fastq_path = paths_cfg.get("fastqs") or paths_cfg.get("data") or FASTQ_DIR
    if not os.path.exists(fastq_path):
        return []

    samples = set()

    # Find all R1 files recursively (covers both flat and subdir layouts)
    r1_files = glob.glob(os.path.join(fastq_path, "**", "*READ1*.fastq.gz"), recursive=True)

    for path in r1_files:
        basename = os.path.basename(path)
        if "Undetermined" in basename:
            continue

        # Prefer the pattern <sample>-...-READ1-...; capture what's before '-READ1'
        if "-READ1" in basename:
            name = basename.split("-READ1", 1)[0]
        elif "_R1_" in basename:
            # fallback to Illumina-style names
            name = basename.split("_R1_", 1)[0]
        else:
            # last resort: strip lane/read suffixes
            name = re.sub(r"(.*)(_R?1).*", r"\1", basename)

        # Clean up any trailing separators
        name = name.rstrip("-_.")
        if name and name != "Undetermined":
            samples.add(name)

    return sorted(samples)


# Helper used by the fastqc_raw rule to validate sample detection and locate FASTQ
def fastqc_raw_input(wildcards):
    detected = get_samples_from_fastq_dir()
    if detected and wildcards.sample not in detected:
        raise ValueError(
            f"Sample '{wildcards.sample}' not found in FASTQ directory. Detected: {detected}"
        )
    # locate the FASTQ file for the sample/mate using the detected FASTQ layout
    return find_fastq_for_sample(wildcards.sample, wildcards.mate)


# Locate a FASTQ file for a sample/mate by searching the FASTQ directory.
def find_fastq_for_sample(sample, mate):
    paths_cfg = config.get("paths", {})
    fastq_path = paths_cfg.get("fastqs") or paths_cfg.get("data") or FASTQ_DIR
    # Try a few variants of the sample name: raw, and cleaned (strip common suffixes
    # that may have been appended by downstream rules, e.g. '.sorted' or '.sorted.markdupes')
    candidates = [sample]
    # remove trailing known suffixes (like .sorted, .sorted.markdupes, .sorted.markdupes.bam etc)
    cleaned = re.sub(r"(\.sorted(?:\..*)?$)|(\.sorted\..*)$", "", sample)
    if cleaned != sample:
        candidates.append(cleaned)

    # Also try stripping any file-extension-like suffix (after first dot)
    if '.' in sample:
        base = sample.split('.', 1)[0]
        if base not in candidates:
            candidates.append(base)

    for s in candidates:
        # Use the same matching strategy as sample detection: look for READ1/READ2 markers
        pattern = os.path.join(fastq_path, "**", f"{s}*READ{mate}*fastq.gz")
        matches = glob.glob(pattern, recursive=True)
        if not matches:
            # fall back to Illumina-style pattern
            pattern2 = os.path.join(fastq_path, "**", f"{s}*_R{mate}_*.fastq.gz")
            matches = glob.glob(pattern2, recursive=True)
        if matches:
            return sorted(matches)[0]

    # If we reach here, nothing was found; provide a helpful error listing attempted candidates
    raise ValueError(
        f"No FASTQ found for sample '{sample}' (tried: {candidates}) mate {mate} in {fastq_path}"
    )

# Extract sample list and configuration
# Use auto-detected samples if available, otherwise fall back to config
auto_samples = get_samples_from_fastq_dir()
SAMPLES = auto_samples if auto_samples else config.get("samples", [])


rule all:
    input:
        # expand(f'{OUTPUT_DIR}/mosdepth/{{sample}}.sorted.coverage.txt', sample=SAMPLES),
        expand(f'{OUTPUT_DIR}/bismark/{{sample}}.sorted.markdupes.bai', sample=SAMPLES),
        expand(f'{OUTPUT_DIR}/trimmed/{{sample}}_R{{mate}}_val_{{mate}}_fastqc.{{ext}}', sample=SAMPLES, mate=[1, 2], ext=['html', 'zip']),
        expand(f'{OUTPUT_DIR}/methyldackel/{{sample}}.sorted_{{context}}.{{ext}}',
        #  sample=SAMPLES, context=['CpG'],
         sample=SAMPLES, context=['CpG', 'CHG', 'CHH'],
         ext=['bedGraph', 'methylKit'])


# Run fastqc on the raw .fastq.gz files
rule fastqc_raw:
    input:
        fastqc_raw_input
    output:
        f'{OUTPUT_DIR}/fastqc_raw/{{sample}}-READ{{mate}}-Sequences_fastqc.html',
        f'{OUTPUT_DIR}/fastqc_raw/{{sample}}-READ{{mate}}-Sequences_fastqc.zip'
    conda: "bsseq"
    resources:
        mem_mb = 4000,
        cpus = config["params"]["cpus"],
        partition = "standard",
        account = "sbsandme_lab"
    params:
        fastqc_path = config['tools'].get('fastqc_path', 'fastqc'),
        out_dir = f'{OUTPUT_DIR}/fastqc_raw/'
    shell:
        '''
        module load fastqc/0.11.9
        fastqc -o {params.out_dir} {input}
        module unload fastqc/0.11.9
        '''


# Trim the read pairs using Trimmomatic (replaces Trim Galore)
rule trimmomatic:
    input:
        f'{OUTPUT_DIR}/fastqc_raw/{{sample}}-READ1-Sequences_fastqc.html',
        f'{OUTPUT_DIR}/fastqc_raw/{{sample}}-READ1-Sequences_fastqc.zip',
        f'{OUTPUT_DIR}/fastqc_raw/{{sample}}-READ2-Sequences_fastqc.html',
        f'{OUTPUT_DIR}/fastqc_raw/{{sample}}-READ2-Sequences_fastqc.zip',
        R1 = lambda wildcards: find_fastq_for_sample(wildcards.sample, 1),
        R2 = lambda wildcards: find_fastq_for_sample(wildcards.sample, 2)
    output:
        pair1 = f'{OUTPUT_DIR}/trimmed/{{sample}}_R1_val_1.fq.gz',
        report1 = f'{OUTPUT_DIR}/trimmed/{{sample}}_R1.fastq.gz_trimming_report.txt',
        pair2 = f'{OUTPUT_DIR}/trimmed/{{sample}}_R2_val_2.fq.gz',
        report2 = f'{OUTPUT_DIR}/trimmed/{{sample}}_R2.fastq.gz_trimming_report.txt'
    conda: "bsseq"
    threads: 8
    resources:
        mem_mb = 4000,
        cpus = config["params"]["cpus"],
        partition = "standard",
        account = "sbsandme_lab"
    params:
        trimmomatic_jar = config['tools'].get('trimmomatic'),
        adapters = config.get('references', {}).get('adapters'),
        illumina_clip = config.get('params', {}).get('trimmomatic', {}).get('illuminaclip'),
        sliding_window = config.get('params', {}).get('trimmomatic', {}).get('sliding_window'),
        min_length = config.get('params', {}).get('trimmomatic', {}).get('min_length')
    shell:
        '''
        # Run Trimmomatic PE and capture its stderr (which includes the summary) into report1
        java -jar {params.trimmomatic_jar} PE -threads {threads} \
            {input.R1} {input.R2} \
            {output.pair1}.tmp {output.pair1}.unpaired.tmp {output.pair2}.tmp {output.pair2}.unpaired.tmp \
            ILLUMINACLIP:{params.adapters}:{params.illumina_clip} \
            SLIDINGWINDOW:{params.sliding_window} MINLEN:{params.min_length} \
            2> {output.report1}

        # Compress paired outputs to match expected .fq.gz filenames
        gzip -c {output.pair1}.tmp > {output.pair1}
        gzip -c {output.pair2}.tmp > {output.pair2}

        # Mirror the report for R2 so downstream rules that expect two report files still work
        cp {output.report1} {output.report2}

        # Clean up temporary files
        rm -f {output.pair1}.tmp {output.pair2}.tmp {output.pair1}.unpaired.tmp {output.pair2}.unpaired.tmp
        '''


# Run fastqc on the trimmmed reads
rule fastqc_trimmmed:
    input:
        f'{OUTPUT_DIR}/trimmed/{{sample}}_R{{mate}}.fastq.gz_trimming_report.txt',
        fq_gz = f'{OUTPUT_DIR}/trimmed/{{sample}}_R{{mate}}_val_{{mate}}.fq.gz'
    output:
        f'{OUTPUT_DIR}/trimmed/{{sample}}_R{{mate}}_val_{{mate}}_fastqc.html',
        f'{OUTPUT_DIR}/trimmed/{{sample}}_R{{mate}}_val_{{mate}}_fastqc.zip'
    conda: "bsseq"
    resources:
        mem_mb = 4000,
        cpus = config["params"]["cpus"],
        partition = "standard",
        account = "sbsandme_lab"
    params:
        fastqc_path = config['tools'].get('fastqc_path', 'fastqc'),
        out_dir = f'{OUTPUT_DIR}/trimmed/'
    shell:
        '''
        module load fastqc/0.11.9
        fastqc -o {params.out_dir} {input.fq_gz}
        module unload fastqc/0.11.9
        '''


# Align to the reference
rule bismark_align:
    input:
        R1 = f'{OUTPUT_DIR}/trimmed/{{sample}}_R1_val_1.fq.gz',
        R2 = f'{OUTPUT_DIR}/trimmed/{{sample}}_R2_val_2.fq.gz'
    output:
        f'{OUTPUT_DIR}/bismark/{{sample}}_R1_val_1_bismark_bt2_pe.bam'
    threads: 8
    conda: "bsseq"
    resources:
        mem_mb = 32000,
        cpus = config["params"]["cpus"],
        partition = "standard",
        account = "sbsandme_lab"
    params:
        bismark_path = config['tools'].get('bismark_path', 'bismark'),
        genome_dir = config["bismark_genome"],
        output_dir = f'{OUTPUT_DIR}/bismark/'
    shell:
        '''
        module load samtools/1.15.1
        module load bismark/0.23.1
        # Run Bismark paired-end alignment using Bowtie2 and write output to a temporary folder
        bismark -p {threads} {params.genome_dir} \
            -1 {input.R1} -2 {input.R2} --output_dir {params.output_dir}
        module unload samtools/1.15.1
        module unload bismark/0.23.1
        '''


# Sort the output files
rule samtools_sort:
    input:
        f'{OUTPUT_DIR}/bismark/{{sample}}_R1_val_1_bismark_bt2_pe.bam'
    output:
        f'{OUTPUT_DIR}/bismark/{{sample}}.sorted.bam'
    threads:
        config['samtools_sort']['threads']
    conda: "bsseq"
    resources:
        mem_mb = 4000,
        cpus = config["params"]["cpus"],
        partition = "standard",
        account = "sbsandme_lab"
    params:
        mem = config['samtools_sort']['mem']
    shell:
        '''
        module load samtools/1.15.1
        samtools sort \
        -@ {threads} \
        -m {params.mem} \
        -O BAM \
        -T {input}.samtools_sort.tmp \
        -o {output} \
        {input}
        module unload samtools/1.15.1
        '''


# Add or replace read groups using Picard so downstream tools have RG tags
rule add_read_groups:
    input:
        f'{OUTPUT_DIR}/bismark/{{sample}}.sorted.bam'
    output:
        f'{OUTPUT_DIR}/bismark/{{sample}}.sorted.rg.bam'
    conda: "bsseq"
    resources:
        mem_mb = 4000,
        cpus = config["params"]["cpus"],
        partition = "standard",
        account = "sbsandme_lab"
    params:
        PICARD_BINARY = config['tools'].get('picard_path')
    shell:
        '''
        java -jar {params.PICARD_BINARY} AddOrReplaceReadGroups \
        -I {input} \
        -O {output} \
        -RGID 1 \
        -RGLB lib1 \
        -RGPL ILLUMINA \
        -RGPU unit1 \
        -RGSM {wildcards.sample}
        '''


# Mark potential PCR duplicates with Picard Tools
rule mark_dupes:
    input:
        f'{OUTPUT_DIR}/bismark/{{sample}}.sorted.rg.bam'
    output:
        f'{OUTPUT_DIR}/bismark/{{sample}}.sorted.markdupes.bam'
    conda: "bsseq"
    threads: 8
    resources:
        mem_mb = 48000,
        cpus = config["params"]["cpus"],
        partition = "standard",
        account = "sbsandme_lab"
    log:
        f'{OUTPUT_DIR}/bismark/{{sample}}.markdupes.log'
    params:
        PICARD_BINARY = config['tools']['picard_path']
    shell:
        '''
        module load samtools/1.15.1
        java -jar {params.PICARD_BINARY} MarkDuplicates \
        -I {input} \
        -O {output} \
        -M {log}
        module unload samtools/1.15.1
        '''


rule samtools_index:
    input:
        f'{OUTPUT_DIR}/bismark/{{sample}}.sorted.markdupes.bam'
    output:
        f'{OUTPUT_DIR}/bismark/{{sample}}.sorted.markdupes.bai'
    threads:
        config['samtools_index']['threads']
    conda: "bsseq"
    resources:
        mem_mb = 4000,
        cpus = config["params"]["cpus"],
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        '''
        module load samtools/1.15.1
        samtools index \
        -@ {threads} \
        -b \
        {input} \
        {output}
        module unload samtools/1.15.1
        '''


# Run MethylDackel to get the inclusion bounds for methylation calling
rule methyldackel_mbias:
    input:
        f'{OUTPUT_DIR}/bismark/{{sample}}.sorted.markdupes.bai',
        bam = f'{OUTPUT_DIR}/bismark/{{sample}}.sorted.markdupes.bam',
    output:
        f'{OUTPUT_DIR}/methyldackel/{{sample}}.sorted_OB.svg',
        # f'{OUTPUT_DIR}/methyldackel/{{sample}}.sorted_OT.svg',
        mbias = f'{OUTPUT_DIR}/methyldackel/{{sample}}.sorted.mbias'
    threads:
        config['methyldackel']['threads']
    conda: "bsseq"
    resources:
        mem_mb = 32000,
        cpus = config["params"]["cpus"],
        partition = "standard",
        account = "sbsandme_lab"
    params:
        out_prefix = f'{OUTPUT_DIR}/methyldackel/{{sample}}.sorted',
        genome = REFERENCE_GENOME
    shell:
        '''
        MethylDackel mbias \
        --CHG \
        --CHH \
        -@ {threads} \
        {params.genome} \
        ./{input.bam} \
        ./{params.out_prefix} \
        2> ./{output.mbias}
        '''


# Run MethylDackel to extract cytosine stats
rule methyldackel_extract:
    input:
        f'{OUTPUT_DIR}/bismark/{{sample}}.sorted.markdupes.bai',
        bam = f'{OUTPUT_DIR}/bismark/{{sample}}.sorted.markdupes.bam',
        mbias = f'{OUTPUT_DIR}/methyldackel/{{sample}}.sorted.mbias'
    output:
        f'{OUTPUT_DIR}/methyldackel/{{sample}}.sorted_CpG.bedGraph',
        f'{OUTPUT_DIR}/methyldackel/{{sample}}.sorted_CHG.bedGraph',
        f'{OUTPUT_DIR}/methyldackel/{{sample}}.sorted_CHH.bedGraph',
        f'{OUTPUT_DIR}/methyldackel/{{sample}}.sorted_CpG.methylKit',
        f'{OUTPUT_DIR}/methyldackel/{{sample}}.sorted_CHG.methylKit',
        f'{OUTPUT_DIR}/methyldackel/{{sample}}.sorted_CHH.methylKit'
    threads:
        config['methyldackel']['threads']
    conda: "bsseq"
    resources:
        mem_mb = 4000,
        cpus = config["params"]["cpus"],
        partition = "standard",
        account = "sbsandme_lab"
    params:
        out_prefix = f'{OUTPUT_DIR}/methyldackel/{{sample}}.sorted',
        genome = REFERENCE_GENOME
    shell:
        '''
        # Get bounds for inclusion from the .mbias files

        OB=$(cut -d ' ' -f 5 {input.mbias})
        # OT=$(cut -d ' ' -f 7 {input.mbias})

        # Get a MethylKit compatible file

        MethylDackel extract \
        --CHG \
        --CHH \
        --OB $OB \
        --methylKit \
        -@ {threads} \
        -o {params.out_prefix} \
        {params.genome} \
        {input.bam}

                MethylDackel extract \
        --CHG \
        --CHH \
        --OB $OB \
        -@ {threads} \
        -o {params.out_prefix} \
        {params.genome} \
        {input.bam}
        '''


# Get the depth for each sample
rule mosdepth:
    input:
        f'{OUTPUT_DIR}/bismark/{{sample}}.sorted.bai',
        bam = f'{OUTPUT_DIR}/bismark/{{sample}}.sorted.bam'
    output:
        f'{OUTPUT_DIR}/mosdepth/{{sample}}.sorted.mosdepth.global.dist.txt',
        f'{OUTPUT_DIR}/mosdepth/{{sample}}.sorted.mosdepth.summary.txt',
        f'{OUTPUT_DIR}/mosdepth/{{sample}}.sorted.per-base.bed.gz',
        f'{OUTPUT_DIR}/mosdepth/{{sample}}.sorted.per-base.bed.gz.csi'
    threads:
        config['mosdepth']['threads']
    conda: "bsseq"
    resources:
        mem_mb = 4000,
        cpus = config["params"]["cpus"],
        partition = "standard",
        account = "sbsandme_lab"
    params:
        mapping_quality = config['mosdepth']['mapping_quality'],
        mosdepth_path = config['tools']['mosdepth_path'],
        out_prefix = f'{OUTPUT_DIR}/mosdepth/{{sample}}.sorted'
    shell:
        '''
        {params.mosdepth_path} \
        -x \
        -t {threads} \
        -Q {params.mapping_quality} \
        {params.out_prefix} \
        {input.bam}
        '''


# Calculate the coverage from the mosdepth output
rule calc_coverage:
    input:
        bed = f'{OUTPUT_DIR}/mosdepth/{{sample}}.sorted.per-base.bed.gz'
    output:
        f'{OUTPUT_DIR}/mosdepth/{{sample}}.sorted.coverage.txt'
    conda: "bsseq"
    resources:
        mem_mb = 4000,
        cpus = config["params"]["cpus"],
        partition = "standard",
        account = "sbsandme_lab"
    params:
        genome = REFERENCE_GENOME
    shell:
        '''
        src/mosdepth_to_x_coverage.py \
        -f {params.genome} \
        -m {input.bed} \
        > {output}
        '''
