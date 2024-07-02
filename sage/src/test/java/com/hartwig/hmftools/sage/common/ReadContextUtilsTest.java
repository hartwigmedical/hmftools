package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_SEQUENCE_200;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.bam.CigarUtils;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadContextUtilsTest
{
    private static ReadCigarInfo buildReadCigar(
            final SAMRecord read, int leftFlankIndex, int leftCoreIndex, int rightCoreIndex, int rightFlankIndex)
    {
        return ReadCigarInfo.buildReadCigar(
                read.getAlignmentStart(), read.getCigar().getCigarElements(), leftFlankIndex, leftCoreIndex, rightCoreIndex, rightFlankIndex);
    }

    @Test
    public void testReadBasesCigar()
    {
        String readBases = REF_BASES_200.substring(0, 40);
        byte[] baseQuals = buildDefaultBaseQuals(readBases.length());
        String readCigar = "40M";
        SAMRecord read = buildSamRecord(100, readCigar, readBases, baseQuals);

        int flankSize = DEFAULT_FLANK_LENGTH;

        // test 1: basic core in aligned section
        int readFlankStart = 10;
        int readFlankEnd = 34;
        int readCoreStart = readFlankStart + flankSize;
        int readCoreEnd = readFlankEnd - flankSize;

        ReadCigarInfo readCigarInfo = buildReadCigar(read, readFlankStart, readCoreStart, readCoreEnd, readFlankEnd);
        assertEquals("25M", CigarUtils.cigarStringFromElements(readCigarInfo.Cigar));
        assertEquals(110, readCigarInfo.FlankPositionStart);
        assertEquals(120, readCigarInfo.CorePositionStart);
        assertEquals(124, readCigarInfo.CorePositionEnd);
        assertEquals(134, readCigarInfo.FlankPositionEnd);
        assertEquals(10, readCigarInfo.FlankIndexStart);
        assertEquals(34, readCigarInfo.FlankIndexEnd);

        // test 2: flanks ending in soft-clipped bases
        readCigar = "20S10M10S";
        read = buildSamRecord(100, readCigar, readBases, baseQuals);

        // test 2: flanks both going into the soft-clips
        readFlankStart = 10;
        readFlankEnd = 34;
        readCoreStart = readFlankStart + flankSize;
        readCoreEnd = readFlankEnd - flankSize;

        readCigarInfo = buildReadCigar(read, readFlankStart, readCoreStart, readCoreEnd, readFlankEnd);
        assertEquals("10S10M5S", CigarUtils.cigarStringFromElements(readCigarInfo.Cigar));
        assertEquals(90, readCigarInfo.FlankPositionStart);
        assertEquals(100, readCigarInfo.CorePositionStart);
        assertEquals(104, readCigarInfo.CorePositionEnd);
        assertEquals(114, readCigarInfo.FlankPositionEnd);
        assertEquals(10, readCigarInfo.FlankIndexStart);
        assertEquals(34, readCigarInfo.FlankIndexEnd);

        // test 3: flanks straddling deletes
        readCigar = "10M10D20M10D10M";
        read = buildSamRecord(100, readCigar, readBases, baseQuals);

        // elements:
        // 100-109  10M   index 0-9
        // 110-119  10D   no change
        // 120-139  20M   index 10-29
        // 140-149  10D   no change
        // 150-159  10M   index 30-39

        readFlankStart = 8;
        readFlankEnd = 32;
        readCoreStart = readFlankStart + flankSize;
        readCoreEnd = readFlankEnd - flankSize;

        readCigarInfo = buildReadCigar(read, readFlankStart, readCoreStart, readCoreEnd, readFlankEnd);
        assertEquals("2M10D20M10D3M", CigarUtils.cigarStringFromElements(readCigarInfo.Cigar));
        assertEquals(108, readCigarInfo.FlankPositionStart);
        assertEquals(128, readCigarInfo.CorePositionStart);
        assertEquals(132, readCigarInfo.CorePositionEnd);
        assertEquals(152, readCigarInfo.FlankPositionEnd);
        assertEquals(8, readCigarInfo.FlankIndexStart);
        assertEquals(32, readCigarInfo.FlankIndexEnd);


        // test 4: cores end in deletes
        readCigar = "20M10D20M10D20M";
        read = buildSamRecord(100, readCigar, readBases, baseQuals);

        // elements:
        // 100-119  20M   index 0-19
        // 120-129  10D   no change
        // 130-149  20M   index 20-39
        // 150-159  10D   no change
        // 160-179  10M   index 40-59

        readFlankStart = 8;
        readCoreStart = 18;
        readCoreEnd = 42;
        readFlankEnd = 52;

        readCigarInfo = buildReadCigar(read, readFlankStart, readCoreStart, readCoreEnd, readFlankEnd);
        assertEquals("12M10D20M10D13M", CigarUtils.cigarStringFromElements(readCigarInfo.Cigar));
        assertEquals(108, readCigarInfo.FlankPositionStart);
        assertEquals(118, readCigarInfo.CorePositionStart);
        assertEquals(162, readCigarInfo.CorePositionEnd);
        assertEquals(172, readCigarInfo.FlankPositionEnd);
        assertEquals(8, readCigarInfo.FlankIndexStart);
        assertEquals(52, readCigarInfo.FlankIndexEnd);

        // test 4b: cores start in deletes
        // elements:
        // 100-120    21M     index 0-20
        // 120-121    1D      no change
        // 122-146    25M     index 21-45

        readBases = REF_BASES_200.substring(100, 121) + REF_BASES_200.substring(122, 124) + "A" + REF_BASES_200.substring(125, 147);
        readCigar = "21M1D25M";
        read = buildSamRecord(100, readCigar, readBases);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);
        SimpleVariant var = createSimpleVariant(124, "T", "A");
        VariantReadContext readContext = builder.createContext(var, read, 23, REF_SEQUENCE_200);
        assertEquals(122, readContext.CorePositionStart);

        // repeat just testing the cigar building
        readFlankStart = 11;
        readFlankEnd = 35;
        readCoreStart = readFlankStart + flankSize;
        readCoreEnd = readFlankEnd - flankSize;

        readCigarInfo = buildReadCigar(read, readFlankStart, readCoreStart, readCoreEnd, readFlankEnd);
        assertEquals(122, readCigarInfo.CorePositionStart);


        // test 5: flanks ending in inserts at both ends
        readCigar = "10M5I20M5I10M";
        read = buildSamRecord(100, readCigar, readBases, baseQuals);

        // elements:
        // 100-109  10M   index 0-9
        // unch     5I    index 10-14
        // 110-129  20M   index 15-34
        // unc      5I    index 35-39
        // 130-139  10M   index 40-49

        readFlankStart = 11;
        readFlankEnd = 36;
        readCoreStart = readFlankStart + flankSize;
        readCoreEnd = readFlankEnd - flankSize;

        readCigarInfo = buildReadCigar(read, readFlankStart, readCoreStart, readCoreEnd, readFlankEnd);
        assertEquals("1M5I20M5I1M", CigarUtils.cigarStringFromElements(readCigarInfo.Cigar));
        assertEquals(109, readCigarInfo.FlankPositionStart);
        assertEquals(116, readCigarInfo.CorePositionStart);
        assertEquals(121, readCigarInfo.CorePositionEnd);
        assertEquals(130, readCigarInfo.FlankPositionEnd);
        assertEquals(9, readCigarInfo.FlankIndexStart);
        assertEquals(40, readCigarInfo.FlankIndexEnd);


        // test 6: core positions ending in inserts
        readFlankStart = 3;
        readCoreStart = 13;
        readCoreEnd = 36;
        readFlankEnd = 46;

        readCigarInfo = buildReadCigar(read, readFlankStart, readCoreStart, readCoreEnd, readFlankEnd);
        assertEquals("7M5I20M5I7M", CigarUtils.cigarStringFromElements(readCigarInfo.Cigar));
        assertEquals(103, readCigarInfo.FlankPositionStart);
        assertEquals(109, readCigarInfo.CorePositionStart);
        assertEquals(130, readCigarInfo.CorePositionEnd);
        assertEquals(136, readCigarInfo.FlankPositionEnd);
        assertEquals(3, readCigarInfo.FlankIndexStart);
        assertEquals(46, readCigarInfo.FlankIndexEnd);
    }
}
