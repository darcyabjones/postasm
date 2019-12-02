#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

def helpMessage() {
    log.info"""
    # qcflow

    """.stripIndent()
}

params.reference = false
params.genes = false
params.assemblies = false
params.reads = false
params.bams = false
params.map = false
params.nobusco = false
params.busco_lineage = false
params.augustus_species = false
params.augustus_params = false
params.nosourmash = false
params.sourmashdb = false


// Find common prefix between 2 strings.
// Used here to emulate basename given from the Channel.fromPairs glob
static String lcp(String r1, String r2){
    def l = [r1.toList(), r2.toList()]
        .transpose()
        .takeWhile {it[0] == it[1]}
        .collect {it[0]}
        .join("")
    return l.replaceAll(/[-\._]$/, "")
}


def checkNotNull = { field, parameter, row  ->
    if (row.name == null || row.name == '') {
        log.info "Encountered empty value in ${field} column of ${parameter} table."
        log.info "Please make sure no required values are empty."
        exit 1
    }
    row
}


if ( params.busco_lineage ) {
    buscoLineage = Channel.fromPath(
        params.busco_lineage,
        checkIfExists: true,
        type: "file"
    )
}

if ( params.sourmashdb ) {
    sourmashDB = Channel
        .fromPath( params.sourmashdb, checkIfExists: true, type: "file" )
        .first()
} else if ( !params.nosourmash ) {
    process getSourmashDB {
        label "download"
        storeDir "${params.outdir}/databases"

        output:
        file "genbank-k31.lca.json.gz" into sourmashDB

        """
        wget -O genbank-k31.lca.json.gz https://osf.io/4f8n3/download
        """
    }
}

if ( params.assemblies ) {
    Channel.fromPath(params.assemblies, checkIfExists: true)
        .splitCsv(by: 1, sep: '\t', header: true)
        .into {
            assemblies4Scaffolds;
            assemblies4Contigs;
            assemblies4ScaffoldsAndContigs;
            assemblies4Mitochondria;
        }

    scaffolds = assemblies4Scaffolds
        .filter { it.scaffolds != null }
        .map { checkNotNull("name", "assemblies", it) }
        .map { [it.name, file(it.scaffolds, checkIfExists: true)] }
        .unique()

    contigs = assemblies4Contigs
        .filter { it.contigs != null }
        .map { checkNotNull("name", "assemblies", it) }
        .map { [it.name, file(it.contigs, checkIfExists: true)] }
        .unique()

    scaffoldsAndContigs = assemblies4ScaffoldsAndContigs
        .filter { it.scaffolds != null && it.contigs != null }
        .map { checkNotNull("name", "assemblies", it) }
        .map {[
            it.name,
            file(it.scaffolds, checkIfExists: true),
            file(it.contigs, checkIfExists: true)
        ]}
        .unique()

    mitochondria = assemblies4Mitochondria
        .filter { it.mitochondria != null }
        .map { checkNotNull("name", "assemblies", it) }
        .map { [it.name, file(it.mitochondria, checkIfExists: true)] }
        .unique()
} else {
    log.info "Hey I need some assemblies to assess please."
    exit 1
}


if ( params.reads ) {
    reads = Channel.fromPath(params.reads, checkIfExists: true)
        .splitCsv(by: 1, sep: '\t', header: true)
        .filter { it.read1_file != null && it.read2_file != null }
        .map { checkNotNull("name", "reads", it) }
        .map {[
            it.name,
            file(it.read1_file, checkIfExists: true),
            file(it.read2_file, checkIfExists: true)
        ]}
        .unique()
} else {
    reads = Channel.empty()
}


if ( params.bams ) {
    bams = Channel.fromPath(params.bams, checkIfExists: true)
        .splitCsv(by: 1, sep: '\t', header: true)
        .filter { it.bam != null }
        .map { checkNotNull("name", "bams", it) }
        .map { [it.name, file(it.bam, checkIfExists: true)] }
        .unique()
} else {
    bams = Channel.empty()
}

// End of input validation


scaffolds.into {
    scaffolds4JoinWithContigs;
    scaffolds4FilterLength;
}

// This is different to scaffoldsAndContigs because it is long format
// rather than wide format.
scaffolds4JoinWithContigs
    .map {n, f -> [n, "scaffolds", f]}
    .concat( contigs.map {n, f -> [n, "contigs", f]} )
    .into {
        scafsCatContigs4AssemblyStats;
        scafsCatContigs4RunBusco;
        scafsCatContigs4AlignReads;
        scafsCatContigs4FindMitochondrial;
    }


mitochondria.into {
    mitochondria4AssemblyStats;
    mitochondria4FindMitochondrial;
}

reads.into {
    reads4AssemblySpectra;
    reads4AlignReads;
}


process assemblyStats {

    label "bbmap"
    tag "${name} - ${type}"
    publishDir "${params.outdir}/asm_stats"

    input:
    set val(name), val(type), file(fasta) from scafsCatContigs4AssemblyStats
        .concat (
            mitochondria4AssemblyStats.map {n, f -> [n, "mitochondria", f]}
        )

    output:
    set val(name), val(type),
        file("${fasta.simpleName}_${type}_gc.txt"),
        file("${fasta.simpleName}_${type}_gchist.txt"),
        file("${fasta.simpleName}_${type}_shist.txt"),
        file("${fasta.simpleName}_${type}_stats.txt") into assemblyStatsResults

    // TODO add extra columns for sample and concatenation step.
    """
    stats.sh \
      in=${fasta} \
      gc=${fasta.simpleName}_${type}_gc.txt \
      gchist=${fasta.simpleName}_${type}_gchist.txt \
      shist=${fasta.simpleName}_${type}_shist.txt \
      extended=t \
      score=t \
      format=3 \
      gcformat=1 \
    | awk -v fname="${name}" \
      'BEGIN {OFS="\t"} NR == 1 {print "name", $0} NR == 2 {print fname, $0}'
    > ${fasta.simpleName}_${type}_stats.txt
    """
}

process combineAssemblyStats {
    label "posix"
    label { type }
    publishDir "${params.outdir}/asm_stats"

    input:
    set val(type), file("*_stats.txt") from assemblyStatsResults
        .map { n, t, gc, gch, sh, st -> [t, st] }
        .groupTuple(by: 0)

    output:
    set val(type), file("${type}_stats.tsv") into combinedAssemblyStatsResults

    """
    array=( *_stats.txt )
    { cat \${array[@]:0:1}; tail -n+2 --quiet \${array[@]:1}; } > ${type}_stats.tsv
    """
}

process alignContigs {
    label "mummer"
    tag { name }
    publishDir "${params.outdir}/asm_align_contigs_to_scaffolds"

    input:
    set val(name), file(scaffolds), file(contigs) from scaffoldsAndContigs

    output:
    set val(name), file("${name}.delta"), file("${name}.coords"), file("${name}.sam") into alignedContigs

    """
    nucmer \
      --maxmatch \
      --threads ${task.cpus} \
      --delta "${name}.delta" \
      --sam-long "${name}.sam" \
      ${scaffolds} \
      ${contigs}

    show-coords -T -c -l -H ${name}.delta > ${name}.coords
    """
}

if ( params.reads ) {

    // We run this several times to see if filtering by contig
    // length results in loss of information.
    // Generally we need to filter out contigs < 200 bp for submission to
    // genbank, beyond that is judgement.
    thresholds = [0, 200, 300, 400, 600, 1000, 50000]

    process filterShortContigs {
        label "seqkit"
        tag { "${name} - ${threshold}" }
        publishDir "${params.outdir}/asm_spectra/filtered_fastas"

        input:
        set val(name), file(asm) from scaffolds4FilterLength
        each threshold from thresholds

        output:
        set val(name), val(threshold),
            file("${asm.simpleName}_filtered_gt${threshold}.fasta") into filtered4KatSpectra

        """
        seqkit seq \
          --min-len ${threshold} \
          "${asm}" \
          > "${asm.simpleName}_filtered_gt${threshold}.fasta"
        """
    }

    joined4KatSpectra = filtered4KatSpectra
        .combine(reads4AssemblySpectra, by: 0)

    // Analyse kmer content in different spectra
    //
    // We could also look within spectra.
    // https://kat.readthedocs.io/en/latest/walkthrough.html#in-assemblies
    process assemblySpectra {
        label "kat"
        tag "${name} - ${thres}"
        publishDir "${params.outdir}/asm_spectra"

        input:
        set val(name), val(thres), file(asm),
            file("*R1.fastq"), file("*R2.fastq") from joined4KatSpectra
                .groupTuple(by: [0, 1, 2])

        output:
        set val(name),
            file("${asm.simpleName}-main.mx"),
            file("${asm.simpleName}-main.mx.spectra-cn.png"),
            file("${asm.simpleName}.stats"),
            file("${asm.simpleName}.dist_analysis.json"),
            file("${asm.simpleName}.1.hist"),
            file("${asm.simpleName}.1.hist.png"),
            file("${asm.simpleName}.2.hist"),
            file("${asm.simpleName}.2.hist.png") into assemblySpectraResults

        """
        # KAT figures out if gzipped or not from the binary, not extension!
        kat comp \
          --threads ${task.cpus} \
          --output_prefix ${asm.simpleName} \
          --output_hists \
          '*R?.fastq' \
          "${asm}"
        """
  }
}

if ( !params.nobusco && params.busco_lineage ) {

    // Be careful using busco to compare assemblies.
    // The results are heavily influenced by the blast database evalues,
    // which in turn depend on the input (i.e. assemblies).
    process runBusco {
        label "busco"
        tag "${name} - ${type}"
        publishDir "${params.outdir}/asm_busco/${type}"

        input:
        set val(name), val(type), file(fasta) from scafsCatContigs4RunBusco
        file "lineage" from buscoLineage

        output:
        file "${fasta.baseName}" into buscoResults

        script:
        if ( params.augustus_species ) {
            species="--species ${params.augustus_species}"
        } else {
            species=""
        }

        if ( params.augustus_params ) {
            aug_params="--augustus_options='${params.augustus_params}'"
        } else {
            aug_params=""
        }

        """
        run_BUSCO.py \
          --in "${fasta}" \
          --out "${fasta.baseName}" \
          --cpu ${task.cpus} \
          --mode "genome" \
          --lineage_path "lineage" \
          ${species} \
          ${aug_params}

        mv "run_${fasta.baseName}" "${fasta.baseName}"
        """
    }
}


if ( params.reads && params.map) {

    /*
     * Align reads to assemblies and collect stats.
     */
    process genomeAlign {
        label "bbmap"
        label "biggish_task"

        publishDir "${params.outdir}/asm_read_aligned/${type}"

        tag "${name} - ${type}"

        input:
        set val(name), val(type), file(asm),
            file("*R1.fastq"), file("*R2.fastq") from scafsCatContigs4AlignReads
            .combine(reads4AlignReads, by: 0)
            .groupTuple(by: [0, 1])

        output:
        set val(name), val(type), file("${asm.simpleName}.sam") into alignedReads
        set val(name), val(type), file("*.txt") into alignedStats

        """
        # Figure out it the file is compressed or not.
        # NB. pigz still shows up as gzip in here.
        # dereference needed because of symlinks
        FWD_FTYPE=\$(file -b --dereference *R1.fastq)
        REV_FTYPE=\$(file -b --dereference *R2.fastq)
        # This keeps only stuff up to first whitespace.
        FWD_FTYPE=\${FWD_FTYPE%% *}
        REV_FTYPE=\${REV_FTYPE%% *}

        if [ "\${FWD_FTYPE}" = "gzip" ]; then
            EXT=".gz"
        elif [ "\${FWD_FTYPE}" = "bzip2" ]; then
            EXT=".bz2"
        elif [ "\${FWD_FTYPE}" = "bzip" ]; then
            EXT=".bz"
        elif [ "\${FWD_FTYPE}" = "XZ" ]; then
            EXT=".xz"
        else
            EXT=""
        fi

        # Combine the reads into a single file.
        cat *R1.fastq > forward.fastq\${EXT}
        cat *R2.fastq > reverse.fastq\${EXT}

        bbmap.sh \
          -Xmx${task.memory.toGiga()}g \
          threads=${task.cpus} \
          in1="forward.fastq\${EXT}" \
          in2="reverse.fastq\${EXT}" \
          ref="${asm}" \
          out="${asm.simpleName}.sam" \
          fast \
          local \
          nodisk \
          covstats="${asm.simpleName}_constats.txt" \
          covhist="${asm.simpleName}_covhist.txt" \
          basecov="${asm.simpleName}_basecov.txt" \
          bincov="${asm.simpleName}_bincov.txt" \
          bhist="${asm.simpleName}_bhist.txt" \
          qhist="${asm.simpleName}_qhist.txt" \
          aqhist="${asm.simpleName}_aqhist.txt" \
          lhist="${asm.simpleName}_lhist.txt" \
          ihist="${asm.simpleName}_ihist.txt" \
          ehist="${asm.simpleName}_ehist.txt" \
          qahist="${asm.simpleName}_qahist.txt" \
          indelhist="${asm.simpleName}_indelhist.txt" \
          mhist="${asm.simpleName}_mhist.txt" \
          gchist="${asm.simpleName}_gchist.txt" \
          idhist="${asm.simpleName}_idhist.txt" \
          scafstats="${asm.simpleName}_scafstats.txt" \
          gcbins=auto \
          idbins=auto
        """
    }


    /*
     * Get additional alignment statistics with samtools.
     */
    process alignmentStats {
        label "samtools"
        label "small_task"

        publishDir "${params.outdir}/asm_read_aligned/${type}"

        tag "${name} - ${type}"

        input:
        set val(name), val(type), file(sam) from alignedReads

        output:
        set val(name),
            val(type),
            file("${sam.simpleName}.idxstats"),
            file("${sam.simpleName}.flagstat"),
            file("${sam.simpleName}.stats") into samtoolsStats

        """
        samtools sort -O bam -o "${sam.simpleName}.bam" "${sam}"
        samtools index "${sam.simpleName}.bam"
        samtools idxstats "${sam.simpleName}.bam" > "${sam.simpleName}.idxstats"
        samtools flagstat "${sam.simpleName}.bam" > "${sam.simpleName}.flagstat"
        samtools stats "${sam.simpleName}.bam" > "${sam.simpleName}.stats"

        rm -f ${sam.simpleName}.{bam,bai}
        """
    }

    joined4AlignmentMultiQC = samtoolsStats
        .flatMap { n, t, i, f, s -> [[t, i], [t, f], [t, s]] }
        .concat(alignedStats.flatMap { n, t, s -> s.collect {f -> [t, f]} })
        .filter { f -> !f.name.endsWith("qhist.txt") }
        .filter { f -> !f.name.endsWith("qahist.txt") }

    /*
     * Produce a multiqc report per reference for the isolates.
     */
    process alignmentMultiQC {
        label "multiqc"
        label "small_task"

        publishDir "${params.outdir}/asm_read_aligned/${type}"

        input:
        set val(type), file("*") from joined4AlignmentMultiQC
            .groupTuple(by: 0)

        output:
        set val(type),
            file("multiqc.html"),
            file("multiqc_data") into alignmentMultiQCResults

        """
        multiqc . --filename "multiqc"
        """
    }
}


if ( !params.nosourmash && params.sourmashdb ) {

    /*
     * Classify the reads using Kraken to detect uncommon contamination.
     */
    process searchSourmash {
        label "sourmash"
        label "medium_task"

        publishDir "${params.outdir}/asm_contaminants"

        tag { name }

        input:
        file "genbank_lca.json.gz" from sourmashDB
        set val(name), file(asm) from assemblies4Sourmash

        output:
        file "${asm.simpleName}.csv" into sourmashResults

        """
        sourmash \
          compute \
          --scaled 1000 \
          -k 21,31,51 \
          --singleton \
          -o "${asm}.sig" \
          "${asm}"

        sourmash \
          lca \
          classify \
          --query "${asm}.sig" \
          --db genbank_lca.json.gz \
          --output "${asm.simpleName}.csv" \
        """
    }
}
