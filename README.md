# qcflow

A [nextflow](www.nextflow.io) pipeline for running multiple illumina
short-read filtering, trimming, and QC steps. QC is performed at each stage.

This is intended for my own use, but feel free to do what ever you like with it.


## Usage

```bash
nextflow run -resume darcyabjones/qcflow \
--fastq "fastq/*_R{1,2}.fastq.gz" \
--adapter1 "data/truseq_fwd.fasta" \
--adapter2 "data/truseq_rev.fasta" \
--scontaminants "data/synth_cont.fasta" \
--references "my_genome*.fasta" \
--map
```

Note that the quotes around globbing patterns are important for some
shells that automatically expand globs (e.g. zsh).

## Mandatory Arguments

```
param                   | description
---------------------------------------------------------------------------
`--fastq <fastq pairs>` | The fastq pairs to process. This must be
                        | provided as a glob pattern to capture the pairs.
                        | E.G. `sample{1,2}.fastq` will capture
                        | `sample1.fastq sample2.fastq`.
```

## Options

```
param                      | default          | description
---------------------------------------------------------------------------
`--outdir <path>`          | ./               | The base directory to
                           |                  | publish the results in.

`--references <*fasta>`    | none             | A glob pattern of reference
                           |                  | genomes to use for mapping
                           |                  | and for masking contaminant
                           |                  | filtering.

`--adapter1 <*fasta>`      | data/            | Adapter sequences to trim
                           | truseq_fwd.fasta | from the 5' end.

`--adapter2 <*fasta>`      | data/            | Adapter sequences to trim
                           | truseq_rev.fasta | from the 3' end.

`--scontaminants <*fasta>` | data/            | Synthetic contaminants to
                           | synth_cont.fasta | filter from the reads.
                           |                  | This is for things like
                           |                  | PHiX or primer dimer, or
                           |                  | common lab vectors that are
                           |                  | likely contaminants of
                           |                  | sequencing rather than of
                           |                  | the samples themselves.

`--contaminants <*fasta>   | none             | Filter reads that match
                           |                  | these sequences out.
                           |                  | This is for filtering out
                           |                  | sample contaminants, e.g.
                           |                  | bacteria or endophytes.
                           |                  | If a reference genome is
                           |                  | provided, regions in this
                           |                  | database matching the
                           |                  | reference will be masked.
                           |                  | Generally I would run kraken
                           |                  | to check for contamination
                           |                  | before using this.

`--filter_phred <int>`     | 5                | Filter out reads that have
                           |                  | lower average phred scores
                           |                  | than this. Keep this low.
                           |                  | bbduk seems to be a bit
                           |                  | "filter happy" for this.

`--trim_phred <int>`       | 2                | Trim bases with phred
                           |                  | qualities lower than this
                           |                  | off the end of reads.
                           |                  | Generally quality trimming
                           |                  | is not useful unless you
                           |                  | have very poor data.

`--use_bbduk_trim`         | false            | Use bbduk for trimming
                           |                  | instead of cutadapt.

`--min_read_length <int>`  | 50               | Filter out reads with
                           |                  | lengths less than this
                           |                  | after trimming.

`--map`                    | false            | Align the raw reads to the
                           |                  | reference genomes and get
                           |                  | qc stats. Useful for
                           |                  | estimating insert/fragment
                           |                  | size, or error rates.

`--merge`                  | false            | Merge paired end reads using
                           |                  | bbmerge. This is intended
                           |                  | for fragment size estimation
                           |                  | when a reference is
                           |                  | unavailable. The parameters
                           |                  | are not appropriate for
                           |                  | merging prior to assembly.

`--krakendb <dir>`         | none             | Search for matches to
                           |                  | potential contaminants in
                           |                  | this kraken database.
                           |                  | This should be seen as a
                           |                  | first pass, check. If you
                           |                  | see something more
                           |                  | substantial you might want
                           |                  | to do a more accurate
                           |                  | alignment.

`--kraken_low_mem`         | false            | Prevents kraken from
                           |                  | memory mapping the file.
                           |                  | This is used in two contexts.
                           |                  | 1) preserve ram, run slow.
                           |                  | 2) See "running kraken".
```

## Running kraken

Kraken is useful as a first pass to detect potential contaminants in your
sequencing. Which you might then choose to filter out using something like
bbduk or bbsplit. It should be noted that the taxonomic assignments given
by kraken are highly dependent of the database you provide, so you should
usually try aligning with a more sensitive aligner (e.g. bma-mem, bbmap),
or database search tool (e.g. blast+) before pushing the contamination red
button.

There is a slurm script to download and prepare a kraken database in the
`batch_scripts` directory in the repo.

Kraken is by far the slowest step in the pipeline if it isn't set up well.
The issue is that for every fastq pair, it tries to load the database into
ram again, which is slow. The best thing to do is to copy the database to
a memory mounted filesystem before running the pipelineand use the
`--kraken_low_mem` option. This way the file stays in memory for the full
time and kraken won't try to load it into memory again.
Many super computers already have memory mounted drives on nodes.
E.G. pawsey mounts `/tmp` in ram, so you can just copy it there as part
of the job submission script.

## Outputs

There's a lot. I'll describe it one day.

## Requirements

* `BBMap` <https://sourceforge.net/projects/bbmap/>.
  Developed with v38.39.
* `fastqc` <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>.
  Developed with v0.11.8.
* `cutadapt` <https://cutadapt.readthedocs.io/en/stable/guide.html>.
  Developed with v1.18
* `multiqc` <https://multiqc.info/>
  Developed with v1.7
* `kraken2` <https://ccb.jhu.edu/software/kraken2/>.
  Developed with v2.0.7-beta.
  Optional, only required for Kraken steps.
