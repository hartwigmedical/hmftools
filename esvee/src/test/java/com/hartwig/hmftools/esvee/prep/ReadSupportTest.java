package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.DEFAULT_BASE_QUAL;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.DEFAULT_MAP_QUAL;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.buildFlags;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;
import static com.hartwig.hmftools.esvee.TestUtils.makeCigarString;
import static com.hartwig.hmftools.esvee.TestUtils.readIdStr;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_ALIGNMENT_BASES;
import static com.hartwig.hmftools.esvee.prep.TestUtils.READ_FILTERS;
import static com.hartwig.hmftools.esvee.prep.JunctionTracker.hasOtherJunctionSupport;
import static com.hartwig.hmftools.esvee.prep.JunctionTracker.hasExactJunctionSupport;
import static com.hartwig.hmftools.esvee.prep.types.ReadFilterType.INSERT_MAP_OVERLAP;
import static com.hartwig.hmftools.esvee.prep.types.ReadFilterType.MIN_ALIGN_MATCH;
import static com.hartwig.hmftools.esvee.prep.types.ReadFilterType.SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.prep.types.ReadFilters.aboveRepeatTrimmedAlignmentThreshold;
import static com.hartwig.hmftools.esvee.prep.types.ReadFilters.isRepetitiveSectionBreak;

import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.esvee.prep.types.JunctionData;
import com.hartwig.hmftools.esvee.prep.types.ReadFilterType;
import com.hartwig.hmftools.esvee.prep.types.PrepRead;

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
        PrepRead junctionRead = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 230, REF_BASES.substring(0, 100), "30S70M"));

        JunctionData junctionData = new JunctionData(230, REVERSE, junctionRead);

        // exact position match
        PrepRead supportRead = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 230, REF_BASES.substring(29, 100), "1S70M"));

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // offset to an earlier base
        supportRead = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 225, REF_BASES.substring(22, 100), "3S72M"));

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // offset to an later base
        SAMRecord record = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 234, REF_BASES.substring(24, 100), "10S66M");
        supportRead = PrepRead.from(record);

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // testing mismatches - 1 is allowed
        record.getReadBases()[0] = (byte)'A';

        supportRead = PrepRead.from(record);
        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        record.getReadBases()[1] = (byte)'A';

        supportRead = PrepRead.from(record);
        assertFalse(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        record.getBaseQualities()[0] = (byte)10;
        supportRead = PrepRead.from(record);
        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // for positive orientation
        junctionRead = PrepRead.from(createSamRecord(
                "01", CHR_1, 200, REF_BASES.substring(0, 80), "50M30S"));

        junctionData = new JunctionData(249, FORWARD, junctionRead);

        // exact position match
        supportRead = PrepRead.from(createSamRecord(
                "02", CHR_1, 200, REF_BASES.substring(0, 53), "50M3S"));

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        supportRead = PrepRead.from(createSamRecord(
                "02", CHR_1, 200, REF_BASES.substring(0, 79), "49M30S"));

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        supportRead = PrepRead.from(createSamRecord(
                "02", CHR_1, 200, REF_BASES.substring(0, 71), "51M20S"));

        // first soft-clip base in read is index 51, then testing index 51 - 70, actual base range 251 - 270, end pos = 250
        // in the junc read: end pos = 249, end pos index = 49
        // junc read index for read's first soft-clip base = 49 + (251 - 249)
        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // offset to an earlier base
        // read index 45, being the first SC base ought to match the junction read index 45
        supportRead = PrepRead.from(createSamRecord(
                "02", CHR_1, 200, REF_BASES.substring(0, 60), "45M15S"));

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        supportRead = PrepRead.from(createSamRecord(
                "02", CHR_1, 210, REF_BASES.substring(10, 60), "35M15S"));

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // offset to an later base
        record = createSamRecord("02", CHR_1, 200, REF_BASES.substring(0, 65), "55M10S");
        supportRead = PrepRead.from(record);

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // testing mismatches
        record.getReadBases()[58] = (byte)'T';
        record.getReadBases()[62] = (byte)'T';

        supportRead = PrepRead.from(record);
        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        record.getReadBases()[61] = (byte)'T';

        supportRead = PrepRead.from(record);
        assertFalse(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        record.getBaseQualities()[58] = (byte)10;
        record.getBaseQualities()[62] = (byte)10;
        supportRead = PrepRead.from(record);
        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));
    }

    @Test
    public void testHardClipSupport()
    {
        PrepRead junctionRead = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 230, REF_BASES.substring(0, 100), "30S70M"));

        JunctionData junctionData = new JunctionData(230, REVERSE, junctionRead);

        // exact position match
        PrepRead supportRead = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 230, REF_BASES.substring(30, 100), "1H70M"));

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // a few ref bases extend past the junction
        supportRead = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 227, REF_BASES.substring(27, 100), "1H73M"));

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // same again or pos orientation
        junctionRead = PrepRead.from(createSamRecord(
                "01", CHR_1, 200, REF_BASES.substring(0, 80), "50M30S"));

        junctionData = new JunctionData(249, FORWARD, junctionRead);

        // exact position match
        supportRead = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 200, REF_BASES.substring(0, 50), "50M3H"));

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));

        supportRead = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 200, REF_BASES.substring(0, 53), "53M30H"));

        assertTrue(hasExactJunctionSupport(supportRead, junctionData, READ_FILTERS));
    }

    @Test
    public void testReadRepeatTrimming()
    {
        String readBases = "ACGTACGTAA" + REF_BASES.substring(10, 40) + REF_BASES.substring(50, 70) + "GGGGGGGGGG" + REF_BASES.substring(70, 100);

        PrepRead read = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 10, readBases, "10S30M10D20M10I30M"));

        assertTrue(aboveRepeatTrimmedAlignmentThreshold(read, MIN_ALIGNMENT_BASES));

        String repeat1 = "ACGTT";
        String repeat2 = "AAGG";
        readBases = "CGT" + repeat1 + repeat1 + repeat1 + repeat1 + repeat1 + repeat1
                + "TGAT" + repeat2 + repeat2 + repeat2 + repeat2 + repeat2 + repeat2 + repeat2 + repeat2 + "AAT";

        read = PrepRead.from(createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 10, readBases, makeCigarString(readBases, 0, 0)));

        assertFalse(aboveRepeatTrimmedAlignmentThreshold(read, MIN_ALIGNMENT_BASES));
    }

    @Test
    public void testDistantSupportingReads()
    {
        // first negative orientation
        PrepRead junctionRead = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 230, REF_BASES.substring(0, 100), "30S70M"));

        JunctionData junctionData = new JunctionData(230, REVERSE, junctionRead);

        // distant within range, correct orientation
        PrepRead supportRead = createRead(READ_ID_GENERATOR.nextId(), 500, false, true);
        assertTrue(hasOtherJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // too far
        supportRead = createRead(READ_ID_GENERATOR.nextId(), 1500, false, true);
        assertFalse(hasOtherJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // wrong orientation
        supportRead = createRead(READ_ID_GENERATOR.nextId(), 500, false, false);
        assertFalse(hasOtherJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // wrong side of junction
        supportRead = createRead(READ_ID_GENERATOR.nextId(), 220, false, true);
        assertFalse(hasOtherJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // other junction side
        junctionRead = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 1500, REF_BASES.substring(0, 100), "70M30S"));

        junctionData = new JunctionData(1500, FORWARD, junctionRead);

        supportRead = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 1000, CHR_3, 100, true, false, null));
        assertTrue(hasOtherJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // too far
        supportRead = createRead(READ_ID_GENERATOR.nextId(), 100, true, false);
        assertFalse(hasOtherJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // wrong orientation
        supportRead = createRead(READ_ID_GENERATOR.nextId(), 100, true, true);
        assertFalse(hasOtherJunctionSupport(supportRead, junctionData, READ_FILTERS));

        // wrong side of junction
        supportRead = createRead(READ_ID_GENERATOR.nextId(), 1510, true, false);
        assertFalse(hasOtherJunctionSupport(supportRead, junctionData, READ_FILTERS));
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

    private static PrepRead createRead(
            final String readId, int readStart, boolean firstInPair, boolean reversed)
    {
        return PrepRead.from(createSamRecord(
                readId, CHR_1, readStart, "", "100M", buildFlags(firstInPair, reversed, false), DEFAULT_MAP_QUAL, DEFAULT_BASE_QUAL));
    }
}
