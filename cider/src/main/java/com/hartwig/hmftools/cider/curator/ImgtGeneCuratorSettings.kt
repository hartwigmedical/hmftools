package com.hartwig.hmftools.cider.curator

import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.cider.genes.VJGeneType
import com.hartwig.hmftools.common.bam.CigarUtils.cigarElementsFromStr
import com.hartwig.hmftools.common.genome.region.Strand

// some extra settings to do with the gene curator
object ImgtGeneCuratorSettings
{
    const val SPECIES = "Homo sapiens"
    const val IMGT_V_ANCHOR_INDEX = 282
    const val IMGT_ANCHOR_LENGTH = 30
    const val BLASTN_EVALUE_CUTOFF = 1000.0
    const val BLASTN_MAX_MISMATCH = 1

    const val ANCHOR_MISMATCH_MAX = 6

    const val FASTA_REF_CONTEXT = 100
    const val REF_CONTEXT_CHECK = 5
    const val REF_CONTEXT_CHECK_MISMATCH_MAX = 1

    // shift the V anchor index such that it starts at the first base
    val IGKINTR_SEQ = ".".repeat(IMGT_V_ANCHOR_INDEX) + "CACCGCGCTCTTGGGGCAGCCGCCTTGCCGCTAGTGGCCGTGGCCACCCTGTGTCTGCCCGATT"
    val IGKDEL_SEQ = "GGAGCCCTAGTGGCAGCCCAGGGCGACTCCTCATGAGTCTGCAGCTGCATTTTTGCCATATCCACTATTTGGAGTCTGACCTCCCTAGGAAGCCTCCCTGC"

    // extra override for this one gene that does not seem to map nicely
    val genomicLocationOverrides = mapOf(
        "IGHV3-54" to LocationInfo(GenomicLocation("chr14", 106601338, 106601641, Strand.REVERSE), cigarElementsFromStr("1X14=1X4=1X2=1X100=1X12=1X41=1X73=8D43="))
    )

    // following genes liftover from v38 to v37 produce incorrect genomic locations
    val liftOverBlacklist = setOf(
        "IGKV1/OR1-1"
    )

    fun getGenomicLocationOverrides(geneName: String): LocationInfo?
    {
        return genomicLocationOverrides[geneName]
    }

    fun jAnchorSignatures(geneName: String) : List<String>
    {
        return if (geneName.startsWith("IGHJ"))
            listOf("TGGGG")
        else if (geneName == VJGeneType.IGKDEL)
            listOf("GCCC")
        else
            listOf("TTTG", "TTCG")
    }
}