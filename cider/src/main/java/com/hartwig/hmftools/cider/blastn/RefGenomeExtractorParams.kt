package com.hartwig.hmftools.cider.blastn

import com.beust.jcommander.Parameter

class RefGenomeExtractorParams
{
    @Parameter(
        names = ["-config"],
        required = true,
        description = "Path to the config file."
    )
    lateinit var config: String

    @Parameter(
        names = ["-input_ref_genome_regions"],
        description = "TSV file of ref genome regions found."
    )
    var inputRefGenomeRegions: String? = null

    @Parameter(
        names = ["-output_ref_genome_regions"],
        required = true,
        description = "TSV file of ref genome regions found."
    )
    lateinit var outputRefGenomeRegions: String

    @Parameter(
        names = ["-blast"],
        required = true,
        description = "Location of blast installation")
    lateinit var blast: String

    @Parameter(names = ["-blast_db"], required = true, description = "Location of blast database")
    lateinit var blastDb: String

    @Parameter(
        names = ["-temp_dir"],
        required = true,
        description = "temp directory to store blastn input / output."
    )
    lateinit var tempDir: String

    @Parameter(names = ["-threads"], description = "Number of threads")
    var threadCount = 1
}