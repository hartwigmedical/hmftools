package com.hartwig.hmftools.cider.curator

import com.hartwig.hmftools.cider.VJGeneType
import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.common.genome.region.Strand

// some extra settings to do with the gene curator
object ImgtGeneCuratorSettings
{
    const val SPECIES = "Homo sapiens"
    const val IMGT_V_ANCHOR_INDEX = 282
    const val IMGT_ANCHOR_LENGTH = 30
    const val BLASTN_EVALUE_CUTOFF = 1000.0
    const val BLASTN_MAX_MISMATCH = 1

    // shift the V anchor index such that it starts at the first base
    val IGKINTR_SEQ = ".".repeat(IMGT_V_ANCHOR_INDEX) + "CACCGCGCTCTTGGGGCAGCCGCCTTGCCGCTAGTGGCCGTGGCCACCCTGTGTCTGCCCGATT"
    val IGKDEL_SEQ = "GGAGCCCTAGTGGCAGCCCAGGGCGACTCCTCATGAGTCTGCAGCTGCATTTTTGCCATATCCACTATTTGGAGTCTGACCTCCCTAGGAAGCCTCCCTGC"

    // extra override for this one gene that does not seem to map nicely
    val genomicLocationOverrides = mapOf(
        "IGHV3-54" to GenomicLocation("chr14", 106601338, 106601641, Strand.REVERSE)
    )

    // following genes liftover from v38 to v37 produce incorrect genomic locations
    val liftOverBlacklist = setOf(
        "IGKV1/OR1-1"
    )

    fun getGenomicLocationOverrides(geneName: String): GenomicLocation?
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