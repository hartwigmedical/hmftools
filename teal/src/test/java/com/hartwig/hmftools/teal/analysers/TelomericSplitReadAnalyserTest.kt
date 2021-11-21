package com.hartwig.hmftools.teal.analysers

import com.hartwig.hmftools.teal.breakend.SplitReadAnalyser
import htsjdk.samtools.SAMRecord
import com.hartwig.hmftools.teal.breakend.*
import com.hartwig.hmftools.teal.TeloUtils
import junit.framework.TestCase
import org.junit.Test

class TelomericSplitReadAnalyserTest
{
    @Test
    fun testFindTelomericBreakEnd()
    {
        val analyser = SplitReadAnalyser(0.9, 0.8, 10)
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
        var tbe = analyser.findTelomericBreakEnd(record)
        TestCase.assertEquals(TelomericBreakEndType.LEFT_C_TELOMERIC, tbe!!.type)
        TestCase.assertEquals("1", tbe.chromosome)
        TestCase.assertEquals(999, tbe.position)

        // try the negative strand
        record.readNegativeStrandFlag = true
        record.readString = TeloUtils.reverseComplementSequence(record.readString)
        tbe = analyser.findTelomericBreakEnd(record)
        TestCase.assertEquals(TelomericBreakEndType.LEFT_C_TELOMERIC, tbe!!.type)
        TestCase.assertEquals("1", tbe.chromosome)
        TestCase.assertEquals(999, tbe.position)

        // right soft clip
        record.cigarString = "30M20S"
        // last 18 is telomeric
        record.readString = "AGCTAGCTAGCTACGTTGCAACCTGACGAA" + "TAACCCTAACCCTAACCC"
        record.readNegativeStrandFlag = false
        tbe = analyser.findTelomericBreakEnd(record)
        TestCase.assertEquals(TelomericBreakEndType.RIGHT_C_TELOMERIC, tbe!!.type)
        TestCase.assertEquals("1", tbe.chromosome)
        TestCase.assertEquals(1030, tbe.position)
        record.readNegativeStrandFlag = true
        record.readString = TeloUtils.reverseComplementSequence(record.readString)
        tbe = analyser.findTelomericBreakEnd(record)
        TestCase.assertEquals(TelomericBreakEndType.RIGHT_C_TELOMERIC, tbe!!.type)
        TestCase.assertEquals("1", tbe.chromosome)
        TestCase.assertEquals(1030, tbe.position)

        // test G telomeric
        record.readNegativeStrandFlag = false
        record.readString = "TTAGGGTTAGGGTTAGGG" + "AGCTAGCTAGCTACGTTGCAACCTGACGAA"
        record.cigarString = "20S30M"
        tbe = analyser.findTelomericBreakEnd(record)
        TestCase.assertEquals(TelomericBreakEndType.LEFT_G_TELOMERIC, tbe!!.type)
        TestCase.assertEquals("1", tbe.chromosome)
        TestCase.assertEquals(999, tbe.position)
    }

    @Test
    fun testLikelyTelomeric()
    {
        val readBasesG = "AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGTGTTAGGG"
        val readBasesC = "CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCTAACCCAAACCCTAACCCAAACCCTAACCCTAACCCTAAC"
        TestCase.assertTrue(TeloUtils.isLikelyGTelomeric(readBasesG))
        TestCase.assertFalse(TeloUtils.isLikelyCTelomeric(readBasesG))
        TestCase.assertFalse(TeloUtils.isLikelyGTelomeric(readBasesC))
        TestCase.assertTrue(TeloUtils.isLikelyCTelomeric(readBasesC))

        //readBasesC = "CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCAAACCCAAACCCTAACCCAAACCCACACCCCCACACCAACCCCCACCCCCACCACAACACCCCCCCCCCCCCCCCCCCCACC";
    }
}