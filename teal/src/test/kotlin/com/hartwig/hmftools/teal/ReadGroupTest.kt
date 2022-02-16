package com.hartwig.hmftools.teal

import com.hartwig.hmftools.teal.ReadGroup.Companion.suppAlignmentPositions
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMTag
import junit.framework.TestCase
import kotlin.test.*

class ReadGroupTest
{
    @Test
    fun testSuppAlignmentPositions()
    {
        var saString = "16,59998,+,64S38M49S,55,0;"
        var sas = suppAlignmentPositions(true, saString)
        TestCase.assertEquals(sas!!.size, 1)
        TestCase.assertEquals(sas[0].chromosome, "16")
        TestCase.assertEquals(sas[0].position, 59998)
        saString = "16,59998,+,64S38M49S,55,0;15,42243201,+,108S35M8S,9,0;,"
        sas = suppAlignmentPositions(true, saString)
        TestCase.assertEquals(sas!!.size, 2)
        TestCase.assertEquals(sas[0].chromosome, "16")
        TestCase.assertEquals(sas[0].position, 59998)
        TestCase.assertEquals(sas[1].chromosome, "15")
        TestCase.assertEquals(sas[1].position, 42243201)
    }

    @Test
    fun testReadGroupIsComplete()
    {
        val rg = ReadGroup("ReadGroup")
        var record = SAMRecord(null)
        record.readName = rg.name
        record.referenceName = "1"
        record.alignmentStart = 100
        record.mateReferenceName = "2"
        record.mateAlignmentStart = 200
        record.readPairedFlag = true
        record.firstOfPairFlag = true
        rg.mutableReads.add(record)
        TestCase.assertFalse(rg.isComplete())
        record = SAMRecord(null)
        record.readName = rg.name
        record.referenceName = "2"
        record.alignmentStart = 200
        record.mateReferenceName = "1"
        record.mateAlignmentStart = 100
        record.readPairedFlag = true
        record.firstOfPairFlag = false
        rg.mutableReads.add(record)
        TestCase.assertTrue(rg.isComplete())
        TestCase.assertTrue(rg.invariant())
    }

    @Test
    fun testReadGroupIsCompleteWithSupplementary()
    {
        val rg = ReadGroup("ReadGroup")
        var record = SAMRecord(null)
        val cigarRead1 = "64S38M49S"
        val cigarSupplRead1 = "13S43M26S"
        val cigarSupplRead2 = "66M43S"
        record.readName = rg.name
        record.referenceName = "1"
        record.alignmentStart = 100
        record.mateReferenceName = "2"
        record.mateAlignmentStart = 200
        record.readPairedFlag = true
        record.firstOfPairFlag = true
        record.readNegativeStrandFlag = true
        rg.mutableReads.add(record)
        TestCase.assertFalse(rg.isComplete())
        record = SAMRecord(null)
        record.readName = rg.name
        record.referenceName = "2"
        record.alignmentStart = 200
        record.mateReferenceName = "1"
        record.mateAlignmentStart = 100
        record.readPairedFlag = true
        record.firstOfPairFlag = false
        record.readNegativeStrandFlag = false
        record.cigarString = cigarRead1
        record.setAttribute(SAMTag.SA.name, "3,300,+,${cigarSupplRead1},55,0;4,400,-,${cigarSupplRead2},55,0;,")
        rg.mutableReads.add(record)
        TestCase.assertFalse(rg.isComplete())

        // add supplementary
        record = SAMRecord(null)
        record.readName = rg.name
        record.referenceName = "3"
        record.alignmentStart = 300
        record.mateReferenceName = "1"
        record.mateAlignmentStart = 100
        record.readPairedFlag = true
        record.firstOfPairFlag = false
        record.supplementaryAlignmentFlag = true
        record.readNegativeStrandFlag = false
        record.cigarString = cigarSupplRead1
        record.setAttribute(SAMTag.SA.name, "2,200,+,${cigarRead1},55,0;")
        rg.mutableSupplementaryReads.add(record)

        // we still missing one supplementary read
        TestCase.assertFalse(rg.isComplete())
        record = SAMRecord(null)
        record.readName = rg.name
        record.referenceName = "4"
        record.alignmentStart = 400
        record.mateReferenceName = "1"
        record.mateAlignmentStart = 100
        record.readPairedFlag = true
        record.firstOfPairFlag = false
        record.supplementaryAlignmentFlag = true
        record.readNegativeStrandFlag = true
        record.cigarString = cigarSupplRead2
        record.setAttribute(SAMTag.SA.name, "2,200,+,${cigarRead1},55,0;")
        rg.mutableSupplementaryReads.add(record)
        TestCase.assertTrue(rg.isComplete())
        TestCase.assertTrue(rg.invariant())
    }
}