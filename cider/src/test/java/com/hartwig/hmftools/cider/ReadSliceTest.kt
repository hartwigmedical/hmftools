package com.hartwig.hmftools.cider

import htsjdk.samtools.SAMRecord
import org.junit.Test
import kotlin.test.assertEquals

class ReadSliceTest
{
    @Test
    fun testReadPosToSlicePos()
    {
        // create a SAM record
        val record = SAMRecord(null)
        record.readString = "AGCTAGCTAG"
        record.readNegativeStrandFlag = false
        record.mappingQuality = 20
        record.duplicateReadFlag = false
        record.readUnmappedFlag = false
        record.properPairFlag = true
        record.readPairedFlag = true

        // make a slice
        val readSlice = ReadSlice(record, false, 2, 7)
        assertEquals(5, readSlice.readLength)
        assertEquals("CTAGC", readSlice.readString)

        assertEquals(1, readSlice.readPositionToSlicePosition(3))
        assertEquals(3, readSlice.slicePositionToReadPosition(1))
    }

    @Test
    fun testReadPosToSlicePosReverseComp()
    {
        // create a SAM record
        val record = SAMRecord(null)
        record.readString = "CTAGCTAGCT" // reverse comp = AGCTAGCTAG
        record.readNegativeStrandFlag = false
        record.mappingQuality = 20
        record.duplicateReadFlag = false
        record.readUnmappedFlag = false
        record.properPairFlag = true
        record.readPairedFlag = true

        // make a slice
        val readSlice = ReadSlice(record, true, 2, 7)
        assertEquals(5, readSlice.readLength)
        assertEquals("CTAGC", readSlice.readString)

        assertEquals(4, readSlice.readPositionToSlicePosition(3))
        assertEquals(3, readSlice.slicePositionToReadPosition(4))
        assertEquals('C', readSlice.baseAt(4))
        assertEquals(2, readSlice.readPositionToSlicePosition(5))
        assertEquals(5, readSlice.slicePositionToReadPosition(2))
        assertEquals('A', readSlice.baseAt(2))

    }
}