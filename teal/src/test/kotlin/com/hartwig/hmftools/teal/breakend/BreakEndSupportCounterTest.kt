package com.hartwig.hmftools.teal.breakend
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.teal.ReadGroup
import htsjdk.samtools.SAMRecord
import junit.framework.TestCase
import kotlin.test.*

class BreakEndSupportCounterTest
{
    @Test
    fun testBreakEndSupport1()
    {
        val rg = ReadGroup("ReadGroup")

        val breakEnd = TelomericBreakEnd(TelomericBreakEndType.LEFT_G_TELOMERIC, "X", 100)

        // left G telomeric
        // first read is break
        var record = SAMRecord(null)
        record.readName = rg.name
        record.referenceName = breakEnd.chromosome
        record.alignmentStart = breakEnd.position
        record.mateReferenceName = "2"
        record.mateAlignmentStart = 200
        record.readPairedFlag = true
        record.firstOfPairFlag = true
        record.readNegativeStrandFlag = true // must be facing left
        record.readString = "TTAGGGTTAGGGCTGAAGCTTGACCG"
        record.cigarString = "12S14M"
        rg.mutableReads.add(record)
        TestCase.assertFalse(rg.isComplete())
        record = SAMRecord(null)
        record.readName = rg.name
        record.referenceName = "2"
        record.alignmentStart = 200
        record.mateReferenceName = breakEnd.chromosome
        record.mateAlignmentStart = breakEnd.position
        record.readPairedFlag = true
        record.firstOfPairFlag = false
        record.readString = "TTAGGGTTAGGGTTA"
        rg.mutableReads.add(record)
        TestCase.assertTrue(rg.isComplete())
        TestCase.assertTrue(rg.invariant())

        val supportCounter = BreakEndSupportCounter(RefGenomeVersion.V37, 0.9)

        val breakEndSupport = supportCounter.countBreakEndSupports(listOf(breakEnd), listOf(rg)).first()

        TestCase.assertTrue(breakEndSupport.fragments.size == 1)

        val fragment = breakEndSupport.fragments.first()

        TestCase.assertSame(fragment.alignedRead, rg.firstOfPair)
        TestCase.assertSame(fragment.pairedRead, rg.secondOfPair)
        TestCase.assertEquals(fragment.alignedReadType, Fragment.AlignedReadType.SPLIT_READ_TELOMERIC)
        TestCase.assertEquals(fragment.pairedReadType, Fragment.PairedReadType.DISCORDANT_PAIR_TELOMERIC)
    }
}