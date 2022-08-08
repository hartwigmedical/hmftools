package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.DEFAULT_BASE_QUAL;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.DEFAULT_MAP_QUAL;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.READ_FILTERS;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.CHR_1;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.buildFlags;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.createSamRecord;
import static com.hartwig.hmftools.svprep.SvPrepTestUtils.readIdStr;
import static com.hartwig.hmftools.svprep.reads.JunctionTracker.hasDiscordantJunctionSupport;
import static com.hartwig.hmftools.svprep.reads.JunctionTracker.hasExactJunctionSupport;
import static com.hartwig.hmftools.svprep.reads.ReadFilterType.INSERT_MAP_OVERLAP;
import static com.hartwig.hmftools.svprep.reads.ReadFilterType.SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.svprep.reads.ReadFilters.isRepetitiveSectionBreak;
import static com.hartwig.hmftools.svprep.reads.ReadRecord.hasPolyATSoftClip;

import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.svprep.reads.JunctionData;
import com.hartwig.hmftools.svprep.reads.ReadFilterType;
import com.hartwig.hmftools.svprep.reads.ReadRecord;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadSupportTest
{
    private static final String REF_BASES =
              "AAAAATTTTTGGGGGCCCCCAATTGGCCATGCAATTGGCCAAAAATTTTTGGGGGCCCCCAATTGGCCATGCAATTGGCC"
            + "AAAAATTTTTGGGGGCCCCCAATTGGCCATGCAATTGGCCAAAAATTTTTGGGGGCCCCCAATTGGCCATGCAATTGGCC"
            + "AAAAATTTTTGGGGGCCCCCAATTGGCCATGCAATTGGCCAAAAATTTTTGGGGGCCCCCAATTGGCCATGCAATTGGCC"
            + "AAAAATTTTTGGGGGCCCCCAATTGGCCATGCAATTGGCCAAAAATTTTTGGGGGCCCCCAATTGGCCATGCAATTGGCC" + generateRandomBases(200);
            // 01234567890123456789012345678901234567890123456789012345678901234567890123456789
            // 200       210       220       230       240       250       260       270
            // 280       290       300       310       320       330       340       350
            // 360       370       380       390       400       410       420       430
            // 440       450       460       470       480       490       500       510

    @Test
    public void testExactSupportingReads()
    {
        int readId = 1;
        ReadRecord junctionRead = ReadRecord.from(createSamRecord(
                readIdStr(readId++), CHR_1, 230, REF_BASES.substring(0, 100), "30S70M"));

        JunctionData junctionData = new JunctionData(230, NEG_ORIENT, junctionRead);

        // exact position match
        ReadRecord supportRead = ReadRecord.from(createSamRecord(
                readIdStr(readId++), CHR_1, 230, REF_BASES.substring(29, 100), "1S70M"));

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // offset to an earlier base
        supportRead = ReadRecord.from(createSamRecord(
                readIdStr(readId++), CHR_1, 225, REF_BASES.substring(22, 100), "3S72M"));

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // offset to an later base
        SAMRecord record = createSamRecord(readIdStr(readId++), CHR_1, 234, REF_BASES.substring(24, 100), "10S66M");
        supportRead = ReadRecord.from(record);

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // testing mismatches - 1 is allowed
        record.getReadBases()[0] = (byte)'A';

        supportRead = ReadRecord.from(record);
        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        record.getReadBases()[1] = (byte)'A';

        supportRead = ReadRecord.from(record);
        assertFalse(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        record.getBaseQualities()[0] = (byte)10;
        supportRead = ReadRecord.from(record);
        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // for positive orientation
        junctionRead = ReadRecord.from(createSamRecord(
                "01", CHR_1, 200, REF_BASES.substring(0, 80), "50M30S"));

        junctionData = new JunctionData(249, POS_ORIENT, junctionRead);

        // exact position match
        supportRead = ReadRecord.from(createSamRecord(
                "02", CHR_1, 200, REF_BASES.substring(0, 53), "50M3S"));

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        supportRead = ReadRecord.from(createSamRecord(
                "02", CHR_1, 200, REF_BASES.substring(0, 79), "49M30S"));

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        supportRead = ReadRecord.from(createSamRecord(
                "02", CHR_1, 200, REF_BASES.substring(0, 71), "51M20S"));

        // first soft-clip base in read is index 51, then testing index 51 - 70, actual base range 251 - 270, end pos = 250
        // in the junc read: end pos = 249, end pos index = 49
        // junc read index for read's first soft-clip base = 49 + (251 - 249)
        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // offset to an earlier base
        // read index 45, being the first SC base ought to match the junction read index 45
        supportRead = ReadRecord.from(createSamRecord(
                "02", CHR_1, 200, REF_BASES.substring(0, 60), "45M15S"));

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        supportRead = ReadRecord.from(createSamRecord(
                "02", CHR_1, 210, REF_BASES.substring(10, 60), "35M15S"));

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // offset to an later base
        record = createSamRecord("02", CHR_1, 200, REF_BASES.substring(0, 65), "55M10S");
        supportRead = ReadRecord.from(record);

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // testing mismatches
        record.getReadBases()[58] = (byte)'T';
        record.getReadBases()[62] = (byte)'T';

        supportRead = ReadRecord.from(record);
        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        record.getReadBases()[61] = (byte)'T';

        supportRead = ReadRecord.from(record);
        assertFalse(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        record.getBaseQualities()[58] = (byte)10;
        record.getBaseQualities()[62] = (byte)10;
        supportRead = ReadRecord.from(record);
        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));
    }

    @Test
    public void testDistantSupportingReads()
    {
        int readId = 1;

        // first negative orientation
        ReadRecord junctionRead = ReadRecord.from(SvPrepTestUtils.createSamRecord(
                readIdStr(readId++), CHR_1, 230, REF_BASES.substring(0, 100), "30S70M"));

        JunctionData junctionData = new JunctionData(230, NEG_ORIENT, junctionRead);

        // distant within range, correct orientation
        ReadRecord supportRead = createRead(readIdStr(readId++), 500, false, true);
        assertTrue(hasDiscordantJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // too far
        supportRead = createRead(readIdStr(readId++), 1500, false, true);
        assertFalse(hasDiscordantJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // wrong orientation
        supportRead = createRead(readIdStr(readId++), 500, false, false);
        assertFalse(hasDiscordantJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // wrong side of junction
        supportRead = createRead(readIdStr(readId++), 220, false, true);
        assertFalse(hasDiscordantJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // other junction side
        junctionRead = ReadRecord.from(SvPrepTestUtils.createSamRecord(
                readIdStr(readId++), CHR_1, 1500, REF_BASES.substring(0, 100), "70M30S"));

        junctionData = new JunctionData(1500, POS_ORIENT, junctionRead);

        supportRead = createRead(readIdStr(readId++), 1000, true, false);
        assertTrue(hasDiscordantJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // too far
        supportRead = createRead(readIdStr(readId++), 100, true, false);
        assertFalse(hasDiscordantJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // wrong orientation
        supportRead = createRead(readIdStr(readId++), 100, true, true);
        assertFalse(hasDiscordantJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // wrong side of junction
        supportRead = createRead(readIdStr(readId++), 1510, true, false);
        assertFalse(hasDiscordantJunctionSupport(supportRead, junctionData, READ_FILTERS));
    }

    @Test
    public void testReadFilterFlags()
    {
        int filters = 0;
        filters = ReadFilterType.set(filters, SOFT_CLIP_LENGTH);
        assertTrue(ReadFilterType.isSet(filters, SOFT_CLIP_LENGTH));
        assertFalse(ReadFilterType.isSet(filters, INSERT_MAP_OVERLAP));

        filters = ReadFilterType.set(filters, INSERT_MAP_OVERLAP);
        assertTrue(ReadFilterType.isSet(filters, SOFT_CLIP_LENGTH));
        assertTrue(ReadFilterType.isSet(filters, INSERT_MAP_OVERLAP));

        filters = ReadFilterType.unset(filters, SOFT_CLIP_LENGTH);
        assertFalse(ReadFilterType.isSet(filters, SOFT_CLIP_LENGTH));
        assertTrue(ReadFilterType.isSet(filters, INSERT_MAP_OVERLAP));

        filters = ReadFilterType.unset(filters, INSERT_MAP_OVERLAP);
        assertFalse(ReadFilterType.isSet(filters, SOFT_CLIP_LENGTH));
        assertFalse(ReadFilterType.isSet(filters, INSERT_MAP_OVERLAP));
    }

    @Test
    public void testRepetitiveBreaks()
    {
        String bases = generateRandomBases(30);

        assertFalse(isRepetitiveSectionBreak(bases.getBytes(), true, 10));
        assertFalse(isRepetitiveSectionBreak(bases.getBytes(), false, 10));

        bases = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

        assertTrue(isRepetitiveSectionBreak(bases.getBytes(), true, 10));
        assertTrue(isRepetitiveSectionBreak(bases.getBytes(), false, 10));

        // 2-base repeats
        bases = "ATATATATATATATATATATATATATATAT";

        assertTrue(isRepetitiveSectionBreak(bases.getBytes(), true, 10));
        assertTrue(isRepetitiveSectionBreak(bases.getBytes(), false, 10));

        // with an error
        bases = "ATATATAGATATATATATATAGATATATAT";

        assertFalse(isRepetitiveSectionBreak(bases.getBytes(), true, 10));
        assertFalse(isRepetitiveSectionBreak(bases.getBytes(), false, 10));

        // 3-base repeats
        bases = "ATCATCATCATCATCATCATCATCATCATCATC";

        assertTrue(isRepetitiveSectionBreak(bases.getBytes(), true, 10));
        assertTrue(isRepetitiveSectionBreak(bases.getBytes(), false, 10));

        // with an error
        bases = "ATCATCATGATCATCATCATCATCGTCATCATC";

        assertFalse(isRepetitiveSectionBreak(bases.getBytes(), true, 10));
        assertFalse(isRepetitiveSectionBreak(bases.getBytes(), false, 10));
    }

    @Test
    public void testPolyATReads()
    {
        String aRepeat = "AAAAAAAAAACAAAAAAA";
        String tRepeat = "TTTTTGTTTTTTTTTTTT";
        String bases = aRepeat + generateRandomBases(30) + tRepeat;

        ReadRecord read = ReadRecord.from(createSamRecord("01",  CHR_1, 100, bases, "18S30M18S"));
        assertTrue(hasPolyATSoftClip(read, true));
        assertTrue(hasPolyATSoftClip(read, false));

        aRepeat = "AAAAACGAAACAAAAAAA";
        tRepeat = "TTTTTGTTTTTAGTTTTT";
        bases = aRepeat + generateRandomBases(30) + tRepeat;
        read = ReadRecord.from(createSamRecord("01",  CHR_1, 100, bases, "18S30M18S"));
        assertFalse(hasPolyATSoftClip(read, true));
        assertFalse(hasPolyATSoftClip(read, false));
    }

    private static ReadRecord createRead(
            final String readId, int readStart, boolean firstInPair, boolean reversed)
    {
        return ReadRecord.from(createSamRecord(
                readId, CHR_1, readStart, "", "100M",
                buildFlags(firstInPair, reversed, false),
                DEFAULT_MAP_QUAL, DEFAULT_BASE_QUAL));
    }
}
