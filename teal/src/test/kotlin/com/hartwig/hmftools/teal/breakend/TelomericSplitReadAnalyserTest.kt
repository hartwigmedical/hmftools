package com.hartwig.hmftools.teal.breakend

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import htsjdk.samtools.SAMRecord
import com.hartwig.hmftools.teal.TealUtils
import kotlin.test.*

/*
    splitTelomereMatchThreshold: Double,
    alignedTelomereMatchThreshold: Double,
    fuzzyMatchDistance: Int,
    maxAlignedPolyBaseThreshold: Double,
 */

class TelomericSplitReadAnalyserTest
{
    @Test
    fun testFindTelomericBreakEnd()
    {
        val config = BreakEndParams()
        config.sampleId = "s"
        config.tumorTelbamFile = ""
        config.germlineTelbamFile = ""
        config.markDuplicateDistance = 10
        config.telomereMatchThreshold = 0.9
        config.alignedSegmentTelomereMatchThreshold = 0.8
        config.maxAlignedPolyBaseThreshold = 0.9
        config.refGenomeVersion = RefGenomeVersion.V37
        val analyser = CandidateBreakEndFinder(config)
        analyser.setBoundaryZone(2)
        //analyser.setMinSoftClipLength(10)
        val record = SAMRecord(null)
        record.readName = "TestRead"
        record.referenceName = "1"
        record.alignmentStart = 1000
        record.cigarString = "20S30M"
        record.mateReferenceName = "2"
        record.mateAlignmentStart = 200
        record.readPairedFlag = true
        record.firstOfPairFlag = true

        // positive strand
        record.readNegativeStrandFlag = false

        // first 18 is telomeric
        record.readString = "TAACCCTAACCCTAACCC" + "AGCTAGCTAGCTACGTTGCAACCTGACGAA"
        var tbe = analyser.findCandidateBreakEnd(record)
        assertEquals(TelomericBreakEndType.LEFT_C_TELOMERIC, tbe!!.type)
        assertEquals("1", tbe.chromosome)
        assertEquals(1000, tbe.position)

        // right soft clip
        record.cigarString = "30M20S"
        // last 18 is telomeric
        record.readString = "AGCTAGCTAGCTACGTTGCAACCTGACGAA" + "TAACCCTAACCCTAACCC"
        record.readNegativeStrandFlag = false
        tbe = analyser.findCandidateBreakEnd(record)
        assertEquals(TelomericBreakEndType.RIGHT_C_TELOMERIC, tbe!!.type)
        assertEquals("1", tbe.chromosome)
        assertEquals(1029, tbe.position)

        // test G telomeric
        record.readNegativeStrandFlag = false
        record.readString = "TTAGGGTTAGGGTTAGGG" + "AGCTAGCTAGCTACGTTGCAACCTGACGAA"
        record.cigarString = "20S30M"
        tbe = analyser.findCandidateBreakEnd(record)
        assertEquals(TelomericBreakEndType.LEFT_G_TELOMERIC, tbe!!.type)
        assertEquals("1", tbe.chromosome)
        assertEquals(1000, tbe.position)
    }

    @Test
    fun testLikelyTelomeric()
    {
        val readBasesG = "AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGTGTTAGGG"
        val readBasesC = "CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCTAACCCAAACCCTAACCCAAACCCTAACCCTAACCCTAAC"
        assertTrue(TealUtils.isLikelyGTelomeric(readBasesG))
        assertFalse(TealUtils.isLikelyCTelomeric(readBasesG))
        assertFalse(TealUtils.isLikelyGTelomeric(readBasesC))
        assertTrue(TealUtils.isLikelyCTelomeric(readBasesC))

        //readBasesC = "CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCAAACCCAAACCCTAACCCAAACCCACACCCCCACACCAACCCCCACCCCCACCACAACACCCCCCCCCCCCCCCCCCCCACC";
    }
}