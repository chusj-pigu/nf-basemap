# nf-basemap

[Nextflow] workflow for basecalling and mapping of whole genome or whole transcriptome Nanopore fastq reads using dorado, minimap2, samtools and multiqc.

## Dependencies

Requires either [Docker] or [Apptainer] installed.

Usage:

```sh
nextflow run chusj-pigu/nf-basemap -r main --sample_sheet [SAMPLESHEET.csv] [--ref REF_GENOME] [OPTIONS]
```

To list available options and parameters for this workflow, run :

``` sh
nextflow run chusj-pigu/nf-basemap -r main --help
```

## Overview

This workflow can be run locally or on Compute Canada. To use on Compute Canada, use the `-profile drac` option. For testing purposes, basecalling can be run on cpu only with `-profile test`.

## Inputs

- `--sample_sheet`: Path to the samplesheet with first column containing sample_id and second column containing the absolute path to fastq or pod5 files (MANDATORY).
- `--ref`: Path to the reference fasta file for alignment.
- `--out_dir`: Output directory to place mapped files and reports in [default: output].
When running on Digital Research Alliance of Canada Narval:
- `--model_path`: Path for the basecalling model, required when running with drac profile [default: path to sup@v5.0.0].
- `--m_bases_path`: Path for the modified basecalling model, required when running with drac profile [default: path to sup@v5.0.0_5mC_5hmC].
- `--ubam`: Path to the samplesheet with first column containing sample_id and second column containing the absolute path to ubam file to append to [default: none].
- `--demux_sheet`: If demux option is selected, pass on a samplesheet specifying sample_id for each barcode. First column must be barcode, second column must be sample_id [default: none].

## Outputs

This workflow will output:

| File Path             | Description | Condition        |
| --------------------- | ----------- | ---------------- |
| {sample_id}_passed.fq.gz | Merged raw reads that have passed filter of average QS >= `--minqs` | If --skip_basecall is not used |
| {sample_id}_failed.fq.gz | Merged raw reads that have failed filter of average QS >= `--minqs` | If --skip_basecall is not used |
| {sample_id}.{ref}.bam<br>{sample_id}.{ref}.bam.bai | Aligned and sorted bam file mapped to reference `--ref` along with it's index | If --skip_mapping is not used |
| multiqc_report.html | [multiqc] report containing [Nanoplot] and [mosdepth] outputs | Always |

## Parameters

- `--model`: Basecalling model to use [default: sup@v5.0.0].
- `--m_bases`: Modified bases to be called, separated by commas if more than one is desired. Requires path to model if run with drac profile [default: 5mC_5hmC].
- `--device`: Parameter to choose device when basecalling. Specify CPU or GPU device: 'auto', 'cpu', 'cuda:all' or 'cuda:<device_id>[,<device_id>...]'. Specifying 'auto' will choose either 'cpu', 'metal' or 'cuda:all' depending on the presence of a GPU device. [nargs=0..1] [default: "auto"]
- `--skip_basecall`: Basecalling step will be skipped; input must be in fastq [default: false].
- `--skip_mapping`: Mapping will be skipped [default: false].
- `--demux`: Option to demultiplex fastq reads, kit used has to be provided (ex. SQK-PCB114-24) [default: null].
- `--duplex`: Dorado will basecall in duplex mode instead of simplex [default: false].
- `--resume`: Use when basecalling has stopped without finishing, and you want to resume dorado from where it stopped [default: false]
- `--no_mod`: Basecalling without base modification [default: false].
- `--minqs`: Minimum average read QS score to keep [default: 10]
- `-profile`: Use standard for running locally, test when running a small local test, drac when running on Digital Research Alliance of Canada Narval and test_drac when running tests on Narval [default: standard].
- `--batch`: Batch size for basecalling; if 0, optimal batch size will be automatically selected [default: 0].

## Steps

### 1. Basecalling

Basecalling is done with [Dorado] in simplex mode with modified base calling by default, using the Super accurate algorithm. The pod5 files are first split by channel to facilitate duplex calling. To use without modified basecalling, use the `--no_mod` option. This can be done in duplex mode using the `--duplex` option, or entirely skipped using the `--skip_basecall` option.
If `--demux KIT` is used, `dorado demux` will be carried out just after basecalling.

### 2. QC stats and filtering

QC plots are generated by [NanoPlot]. Reads that have a QS score <10 will be classified as failed, and subsequent steps are only carried out on reads with QS score >= 10.

### 3. Mapping

If using `--skip_basecall`, the workflow will start here using `--fastq` as input. Mapping of the reads is done with [minimap2] using `-ax -map-ont` options for long reads generated by Nanopore sequencing. Mapping can be skipped with `--skip_mapping`.

### 4. Sorting and indexing

Using [samtools], the .sam file is first converted to .bam using `samtools view`, sorted with `samtools sort` and indexed using `samtools index`. Mapping statistics are computed with [mosdepth].

[Docker]: https://www.docker.com
[Apptainer]: https://apptainer.org
[Nextflow]: https://www.nextflow.io/docs/latest/index.html
[Dorado]: https://github.com/nanoporetech/dorado
[minimap2]: https://lh3.github.io/minimap2/minimap2.html
[samtools]: http://www.htslib.org
[multiqc]: https://multiqc.info
[mosdepth]: https://github.com/brentp/mosdepth
[NanoPlot]: https://github.com/wdecoster/NanoPlot
