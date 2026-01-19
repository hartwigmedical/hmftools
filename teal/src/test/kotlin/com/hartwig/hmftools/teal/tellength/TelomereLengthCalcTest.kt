package com.hartwig.hmftools.teal.tellength

import com.hartwig.hmftools.common.sequencing.SequencingType
import com.hartwig.hmftools.teal.ReadGroup
import htsjdk.samtools.SAMRecord
import junit.framework.TestCase
import kotlin.test.Test

class TelomereLengthCalcTest
{
    @Test
    fun testTelomereLengthCalcZeroRead()
    {
        val telomereLengthCalc = TelomereLengthCalc(
            SequencingType.ILLUMINA, 1.0, 2.0, 0.0,
            100.0, 100.0, 1000.0)
        TestCase.assertEquals(0.0, telomereLengthCalc.calcTelomereLength())
        TestCase.assertEquals(0.0, telomereLengthCalc.calcTumorTelomereLength())
    }

    @Test
    fun testTelomereLengthCalc()
    {
        val telomereLengthCalc = TelomereLengthCalc(SequencingType.ILLUMINA, 1.0, 2.0, 0.0,
            2.0, 2.0, 100.0)

        // create a read group
        val rg = ReadGroup("ReadGroup")
        var record = SAMRecord(null)
        record.readName = rg.name
        record.referenceName = "1"
        record.alignmentStart = 100
        record.mateReferenceName = "2"
        record.mateAlignmentStart = 200
        record.readPairedFlag = true
        record.firstOfPairFlag = true

        // we need at least 6 repeats of TTAGGG like pattern
        record.readString = "TTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG"
        rg.mutableReads.add(record)

        record = SAMRecord(null)
        record.readName = rg.name
        record.referenceName = "2"
        record.alignmentStart = 200
        record.mateReferenceName = "1"
        record.mateAlignmentStart = 100
        record.readPairedFlag = true
        record.firstOfPairFlag = false

        // we need at least 6 repeats of CCCTAA like pattern
        record.readString = "CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA"
        rg.mutableReads.add(record)
        TestCase.assertTrue(rg.isComplete())
        TestCase.assertTrue(rg.invariant())

        telomereLengthCalc.onReadGroup(rg)

        TestCase.assertEquals(0.869143, telomereLengthCalc.calcTelomereLength(), 1e-5)
        TestCase.assertEquals(0.869143, telomereLengthCalc.calcTumorTelomereLength(), 1e-5)
    }
}