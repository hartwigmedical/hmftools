package com.hartwig.hmftools.cider.blastn

import com.hartwig.hmftools.cider.CiderConstants
import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.common.genome.region.Strand

data class BlastnMatch(
    val querySeqLen: Int,
    val subjectTitle: String,
    val percentageIdent: Double,
    val queryCoverage: Double,
    val alignmentLength: Int,
    val numMismatch: Int,
    val numGapOpenings: Int,
    val queryAlignStart: Int,
    val queryAlignEnd: Int,
    val subjectAlignStart: Int,
    val subjectAlignEnd: Int,
    val subjectFrame: Strand,
    val expectedValue: Double,
    val bitScore: Double,
    val alignedPartOfQuerySeq: String,
    val alignedPartOfSubjectSeq: String
)
{
    fun toGenomicLocation() : GenomicLocation?
    {
        // parse the chromosome / assembly
        val m = CHROMOSOME_ASSEMBLY_REGEX.matchEntire(subjectTitle)

        if (m == null)
        {
            // might not have a chromosome, or it is mitochondrion
            return null
        }

        val chromosome = CiderConstants.BLAST_REF_GENOME_VERSION.versionedChromosome(m.groupValues[1]).intern()
        val assembly = m.groupValues[2].intern()

        val start: Int
        val end: Int

        when (subjectFrame)
        {
            Strand.FORWARD -> { start = subjectAlignStart; end = subjectAlignEnd }
            Strand.REVERSE -> { start = subjectAlignEnd; end = subjectAlignStart }
        }

        return GenomicLocation(chromosome, start, end, subjectFrame, if (isPrimaryAssembly(assembly)) null else assembly)
    }

    companion object
    {
        // val REF_GENOME_NAME = "GRCh38.p13"
        val CHROMOSOME_ASSEMBLY_REGEX = Regex("^Homo sapiens chromosome (\\w+).*, GRCh38.p13 (.+)$")
        val PRIMARY_ASSEMBLY_NAME = "Primary Assembly".intern()

        fun isPrimaryAssembly(assemblyName: String) : Boolean
        {
            return assemblyName == PRIMARY_ASSEMBLY_NAME
        }
    }
}
