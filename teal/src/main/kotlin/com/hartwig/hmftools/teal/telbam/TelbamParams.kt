package com.hartwig.hmftools.teal.telbam

import com.beust.jcommander.Parameter
import com.beust.jcommander.converters.IParameterSplitter
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource

class SemiColonSplitter : IParameterSplitter
{
    override fun split(value: String): List<String>
    {
        return value.split(";")
    }
}

data class TelbamParams
(
    @Parameter(names = ["-bam_file"], required = true, description = "Path to bam/cram file")
    var bamFile: String? = null,

    @Parameter(names = ["-" + RefGenomeSource.REF_GENOME],
                        description = "Path to reference genome fasta file if using CRAM files")
    var refGenomeFile: String? = null,

    @Parameter(names = ["-telbam_file"], required = true, description = "Path to write output telbam file")
    var telbamFile: String? = null,

    @Parameter(names = ["-tsv_file"], required = false, description = "Path to write output tsv file")
    var tsvFile: String? = null,

    @Parameter(names = ["-threads"], description = "Number of bam reader threads")
    var threadCount: Int = 1,

    @Parameter(names = ["-specific_chr"],
                splitter = SemiColonSplitter::class,
                description = "Optional: list of chromosomes separated by ;")
    var specificChromosomes: List<String> = emptyList()
)
