package com.hartwig.hmftools.cider.curator

import com.hartwig.hmftools.cider.VJGeneType
import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.region.ChrBaseRegion

// some extra settings to do with the gene curator
object ImgtGeneCuratorSettings
{
    const val SPECIES = "Homo sapiens"
    const val IMGT_V_ANCHOR_INDEX = 282
    const val IMGT_ANCHOR_LENGTH = 30
    const val ALIGNMENT_MAX_MISMATCH = 1

    // shift the V anchor index such that it starts at the first base
    val IGKINTR_SEQ = ".".repeat(IMGT_V_ANCHOR_INDEX) + "CACCGCGCTCTTGGGGCAGCCGCCTTGCCGCTAGTGGCCGTGGCCACCCTGTGTCTGCCCGATT"
    val IGKDEL_SEQ = "GGAGCCCTAGTGGCAGCCCAGGGCGACTCCTCATGAGTCTGCAGCTGCATTTTTGCCATATCCACTATTTGGAGTCTGACCTCCCTAGGAAGCCTCCCTGC"

    // extra override for this one gene that does not seem to map nicely
    val genomicLocationOverrides = mapOf(
        RefGenomeVersion.V38 to mapOf(
            "IGHV3-54" to GenomicLocation(ChrBaseRegion("chr14", 106601338, 106601641), Strand.REVERSE)
        ),
        RefGenomeVersion.V37 to mapOf(
        )
    )

    fun getGenomicLocationOverrides(geneName: String, genomeVersion: RefGenomeVersion): GenomicLocation?
    {
        return genomicLocationOverrides[genomeVersion]?.get(geneName)
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