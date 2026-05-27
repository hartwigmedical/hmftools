package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadKeyTest
{
    @Test
    public void testSupplementarySecondaryIds()
    {
        String readId = "read1";
        String readBases = "ACGT";
        String readCigar = "4M";

        int readStart = 2000;
        int mateStart = 2200;

        SAMRecord read = SamRecordTestUtils.createSamRecord(
                readId, CHR_1, readStart, readBases, readCigar, CHR_1, mateStart,
                false, false, null);
        read.setMappingQuality(20);

        // this is not a supplementary, should have 0 as supplementary index
        ReadKey readKey = ReadKey.from(read);
        assertEquals(readId, readKey.ReadName);
        assertTrue(readKey.FirstInPair);
        assertNull(readKey.SuppSecondaryId);

        read.setSupplementaryAlignmentFlag(true);
        readKey = ReadKey.from(read);
        assertEquals(readId, readKey.ReadName);
        assertTrue(readKey.FirstInPair);
        assertEquals("supp:1:2000:4M:+", readKey.SuppSecondaryId);

        read.setSupplementaryAlignmentFlag(false);
        read.setReadNegativeStrandFlag(true);
        read.setSecondaryAlignment(true);
        readKey = ReadKey.from(read);
        assertEquals(readId, readKey.ReadName);
        assertTrue(readKey.FirstInPair);
        assertEquals("second:1:2000:4M:-", readKey.SuppSecondaryId);

    }
}