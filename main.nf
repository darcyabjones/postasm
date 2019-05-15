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


params.assemblies = false
params.table = false
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

if ( params.table ) {
    table = Channel.fromPath(params.table, checkIfExists: true)
        .splitCsv(by: 1, sep: '\t', header: true)
        .filter { it.assembly != null && it.read1_file != null && it.read2_file != null }
        .map {[
            it.assembly,
            file(it.read1_file, checkIfExists: true),
            file(it.read2_file, checkIfExists: true)
        ]}

    table.into {
        table4AssemblyChannel;
        table4Align;
        table4Kat;
        table4FilterLength;
    }
}


if ( params.assemblies ) {
    assemblies = Channel
        .fromPath(params.assemblies, checkIfExists: true, type: "file")
        .map { file -> [file.name, file] }
} else if ( params.table ) {
    assemblies = table4AssemblyChannel
        .unique { asm, r1, r2 -> asm}
        .map { asm, r1, r2 -> [asm, file(asm, checkIfExists: true)] }
} else {
    log.info "Hey I need some assemblies to assess."
    exit 1
}

assemblies.into {
    assemblies4Stats;
    assemblies4Busco;
    assemblies4Sourmash;
    assemblies4AlignContigs;
    assemblies4FilterCoveredContigs;
    assemblies4Align;
    assemblies4FilterLength;
}


process assemblyStats {

    label "bbmap"
    tag { name }
    publishDir "${params.outdir}/asm_stats"

    input:
    set val(name), file(fasta) from assemblies4Stats

    output:
    set val(name),
        file("${fasta.simpleName}_gc.txt"),
        file("${fasta.simpleName}_gchist.txt"),
        file("${fasta.simpleName}_shist.txt"),
        file("${fasta.simpleName}_stats.txt") into assemblyStatsResults

    // TODO add extra columns for sample and concatenation step.
    """
    stats.sh \
      in=${fasta} \
      gc=${fasta.simpleName}_gc.txt \
      gchist=${fasta.simpleName}_gchist.txt \
      shist=${fasta.simpleName}_shist.txt \
      extended=t \
      score=t \
      format=3 \
      gcformat=1 \
    > ${fasta.simpleName}_stats.txt
    """
}


process alignContigs {
    label "mummer"
    tag { name }
    publishDir "${params.outdir}/asm_self_aligned"

    input:
    set val(name), file(fasta) from assemblies4AlignContigs

    output:
    set val(name), file("${name}.delta"), file("${name}.coords") into alignedContigs

    """
    nucmer \
      --maxmatch \
      --nosimplify \
      --threads ${task.cpus} \
      --prefix ${name} \
      ${fasta} \
      ${fasta}

    show-coords -T -c -l -H ${name}.delta > ${name}.coords
    """
}


process filterCoveredContigs {
    label "python3"
    tag { name }
    publishDir "${params.outdir}/asm_self_aligned"

    input:
    set val(name), file(fasta) from assemblies4FilterCoveredContigs
    set val(name), file("${name}.delta"), file("${name}.coords") from alignedContigs

    output:
    set val(name), file("${name}.summary"), file("${name}.filtered") into filteredContigs


}


if ( params.table ) {

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
        set val(name), file(asm) from assemblies4FilterLength
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
        .combine(table4Kat, by: 0)

    // Analyse kmer content in different spectra
    //
    // We could also look within spectra.
    // https://kat.readthedocs.io/en/latest/walkthrough.html#in-assemblies
    process assemblySpectra {
        label "kat"
        tag { "${name} - ${thres}" }
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
        tag { name }
        publishDir "${params.outdir}/asm_busco"

        input:
        set val(name), file(fasta) from assemblies4Busco
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


if ( params.table && params.map) {

    /*
     * Align reads to assemblies and collect stats.
     */
    process genomeAlign {
        label "bbmap"
        label "biggish_task"

        publishDir "${params.outdir}/asm_self_aligned"

        tag { "${name}" }

        input:
        set val(name), file(asm), file("*R1.fastq"), file("*R2.fastq") from assemblies4Align
            .combine(table4Align, by: 0)
            .groupTuple(by: [0, 1])

        output:
        set val(name), file("${asm.simpleName}.sam") into alignedReads
        set val(name), file("*.txt") into alignedStats

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

        publishDir "${params.outdir}/asm_self_aligned"

        tag { "${name}" }

        input:
        set val(name), file(sam) from alignedReads

        output:
        set val(name),
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
        .flatMap { r, i, f, s -> [i, f, s] }
        .concat(alignedStats.flatMap { r, s -> s })
        .filter { f -> !f.name.endsWith("qhist.txt") }
        .filter { f -> !f.name.endsWith("qahist.txt") }

    /*
     * Produce a multiqc report per reference for the isolates.
     */
    process alignmentMultiQC {
        label "multiqc"
        label "small_task"

        publishDir "${params.outdir}/asm_self_aligned"

        input:
        file "*" from joined4AlignmentMultiQC.collect()

        output:
        set file("multiqc.html"), file("multiqc_data") into alignmentMultiQCResults

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
