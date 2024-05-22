package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.bam.CigarUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadContextUtilsTest
{
    @Test
    public void testReadBasesCigar()
    {
        String readBases = REF_BASES_200.substring(0, 40);
        byte[] baseQuals = buildDefaultBaseQuals(readBases.length());
        String readCigar = "40M";
        SAMRecord read = buildSamRecord(100, readCigar, readBases, baseQuals);

        // test 1: basic core in aligned section
        int readFlankStart = 10;
        int readFlankEnd = 34;

        ReadCigarInfo readCigarInfo = ReadCigarInfo.buildReadCigar(read, readFlankStart, readFlankEnd);
        assertEquals("25M", CigarUtils.cigarStringFromElements(readCigarInfo.Cigar));
        assertEquals(10, readCigarInfo.FlankIndexStart);
        assertEquals(34, readCigarInfo.FlankIndexEnd);
        assertEquals(110, readCigarInfo.UnclippedStart);
        assertEquals(134, readCigarInfo.UnclippedEnd);

        // test 2: flanks ending in soft-clipped bases
        readCigar = "20S10M10S";
        read = buildSamRecord(100, readCigar, readBases, baseQuals);

        // test 2: flanks both going into the soft-clips
        readFlankStart = 10;
        readFlankEnd = 34;

        readCigarInfo = ReadCigarInfo.buildReadCigar(read, readFlankStart, readFlankEnd);
        assertEquals("10S10M5S", CigarUtils.cigarStringFromElements(readCigarInfo.Cigar));
        assertEquals(10, readCigarInfo.FlankIndexStart);
        assertEquals(34, readCigarInfo.FlankIndexEnd);
        assertEquals(90, readCigarInfo.UnclippedStart);
        assertEquals(114, readCigarInfo.UnclippedEnd);

        // test 3: flanks straddling deletes
        readCigar = "10M10D20M10D10M";
        read = buildSamRecord(100, readCigar, readBases, baseQuals);

        readFlankStart = 8;
        readFlankEnd = 32;

        readCigarInfo = ReadCigarInfo.buildReadCigar(read, readFlankStart, readFlankEnd);
        assertEquals("2M10D20M10D3M", CigarUtils.cigarStringFromElements(readCigarInfo.Cigar));
        assertEquals(108, readCigarInfo.UnclippedStart);
        assertEquals(150, readCigarInfo.UnclippedEnd);

        // test 4: flanks ending in inserts at both ends
        readCigar = "10M5I20M5I10M";
        read = buildSamRecord(100, readCigar, readBases, baseQuals);

        readFlankStart = 10;
        readFlankEnd = 35;

        readCigarInfo = ReadCigarInfo.buildReadCigar(read, readFlankStart, readFlankEnd);
        assertEquals("1M5I20M5I1M", CigarUtils.cigarStringFromElements(readCigarInfo.Cigar));
        assertEquals(9, readCigarInfo.FlankIndexStart);
        assertEquals(40, readCigarInfo.FlankIndexEnd);
        assertEquals(109, readCigarInfo.UnclippedStart);
        assertEquals(130, readCigarInfo.UnclippedEnd);
    }

    @Test
    public void testDetermineAltIndexRanges()
    {


    }
}
