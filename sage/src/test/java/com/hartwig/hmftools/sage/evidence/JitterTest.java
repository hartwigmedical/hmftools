package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.basequal.jitter.JitterModelParams.REPEAT_UNIT_3_PLUS_LABEL;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_SEQUENCE_200;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_SAMPLE;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContext;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContextMatcher;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadCounter;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;
import static com.hartwig.hmftools.sage.evidence.JitterMatch.LENGTHENED;
import static com.hartwig.hmftools.sage.evidence.JitterMatch.SHORTENED;
import static com.hartwig.hmftools.sage.evidence.JitterMatch.checkJitter;
import static com.hartwig.hmftools.sage.evidence.JitterMatch.hasJitterMatchType;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.basequal.jitter.JitterModelParams;
import com.hartwig.hmftools.sage.common.ReadContextMatcher;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.RepeatInfo;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;
import com.hartwig.hmftools.sage.quality.MsiJitterCalcs;
import com.hartwig.hmftools.sage.quality.MsiModelParams;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class JitterTest
{
    @Test
    public void testJitterMatching()
    {
        //                           10        20        30
        //                 012345678901234567890123456789012
        String refBases = "GTCGTCGTCGTAAAACAAAAACTCGTCGTCGTC";

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);
        RefSequence refSequence = new RefSequence(100, refBases.getBytes());

        String variantReadBases = refBases.substring(0, 15) + "A" + refBases.substring(16);
        SimpleVariant variant = createSimpleVariant(115, "C", "A");
        SAMRecord read = buildSamRecord(100, buildCigarString(variantReadBases.length()), variantReadBases);

        VariantReadContext readContext = builder.createContext(variant, read, 15, refSequence);
        assertTrue(readContext.isValid());
        assertEquals(15, readContext.VarIndex);
        assertEquals(1, readContext.AllRepeats.size());
        RepeatInfo repeat = readContext.AllRepeats.get(0);
        assertEquals(11, repeat.Index);
        assertEquals(10, repeat.Count);

        byte[] readQuals = buildDefaultBaseQuals(variantReadBases.length());

        ReadContextMatcher matcher = new ReadContextMatcher(readContext, true, true);

        // test jitter without indication of a indel

        // test 1: shortened at start
        String readBases = variantReadBases.substring(0, 11) + "G" + variantReadBases.substring(12);

        assertTrue(hasJitterMatchType(
                repeat, readContext, matcher.altIndexLower(), matcher.altIndexUpper(),
                15, readBases.getBytes(), readQuals, SHORTENED, false, true, 1));

        // test 2: shortened at end
        readBases = variantReadBases.substring(0, 20) + "G" + variantReadBases.substring(21);

        assertTrue(hasJitterMatchType(
                repeat, readContext, matcher.altIndexLower(), matcher.altIndexUpper(),
                15, readBases.getBytes(), readQuals, SHORTENED, false, false, 1));

        // test 3: lengthened at start
        readBases = variantReadBases.substring(0, 10) + "A" + variantReadBases.substring(11);

        assertTrue(hasJitterMatchType(
                repeat, readContext, matcher.altIndexLower(), matcher.altIndexUpper(),
                15, readBases.getBytes(), readQuals, LENGTHENED, false, true, 1));

        // test 4: lengthened at end
        readBases = variantReadBases.substring(0, 21) + "A" + variantReadBases.substring(22);

        assertTrue(hasJitterMatchType(
                repeat, readContext, matcher.altIndexLower(), matcher.altIndexUpper(),
                15, readBases.getBytes(), readQuals, LENGTHENED, false, false, 1));


        // test 5: low-qual base mismatches are permitted outside the core and within the core in specific locations

        String readMismatches = "G" + readBases.substring(1, 22) + "G" + readBases.substring(23);
        readQuals[0] = 11;
        readQuals[22] = 11;

        assertTrue(hasJitterMatchType(
                repeat, readContext, 15, matcher.altIndexLower(), matcher.altIndexUpper(),
                readMismatches.getBytes(), readQuals, LENGTHENED, false, false, 1));

        // 1 mismatch is permitted in the core outside the critical range
        readQuals = buildDefaultBaseQuals(variantReadBases.length());
        readMismatches = readBases.substring(0, 12) + "T" + readBases.substring(13);
        readQuals[12] = 11;

        assertTrue(hasJitterMatchType(
                repeat, readContext, 15, matcher.altIndexLower(), matcher.altIndexUpper(),
                readMismatches.getBytes(), readQuals, LENGTHENED, false, false, 1));

        // 2 mismatches are not permitted
        readQuals = buildDefaultBaseQuals(variantReadBases.length());
        readMismatches = readBases.substring(0, 12) + "TT" + readBases.substring(14);
        readQuals[12] = 11;
        readQuals[13] = 11;

        assertFalse(hasJitterMatchType(
                repeat, readContext, 15, matcher.altIndexLower(), matcher.altIndexUpper(),
                readMismatches.getBytes(), readQuals, LENGTHENED, false, false, 1));

        // and not within the critical range
        readQuals = buildDefaultBaseQuals(variantReadBases.length());
        readMismatches = readBases.substring(0, 15) + "T" + readBases.substring(16);
        readQuals[15] = 11;

        assertFalse(hasJitterMatchType(
                repeat, readContext, 15, matcher.altIndexLower(), matcher.altIndexUpper(),
                readMismatches.getBytes(), readQuals, LENGTHENED, false, false, 1));
    }

    @Test
    public void testDoubleJitter()
    {
        //                           10        20        30
        //                 0123456789012345678901234567890123456789012345678
        String refBases = "GTCGTCGTCGTATAAAAAAAAAACTCGTAGTCGTCATTCGTCATTCATT";

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        int varIndex = 23;
        int reqPosStart = 100;
        int position = reqPosStart + varIndex;
        RefSequence refSequence = new RefSequence(reqPosStart, refBases.getBytes());

        String variantReadBases = refBases.substring(0, varIndex) + "G" + refBases.substring(varIndex + 1);
        SimpleVariant variant = createSimpleVariant(position, "C", "G");
        SAMRecord read = buildSamRecord(reqPosStart, buildCigarString(variantReadBases.length()), variantReadBases);

        VariantReadContext readContext = builder.createContext(variant, read, varIndex, refSequence);
        assertTrue(readContext.isValid());
        assertEquals(21, readContext.VarIndex);
        assertEquals(1, readContext.AllRepeats.size());
        RepeatInfo repeat = readContext.AllRepeats.get(0);
        assertEquals(11, repeat.Index);
        assertEquals(10, repeat.Count);

        ReadContextMatcher matcher = createReadContextMatcher(readContext);

        byte[] readQuals = buildDefaultBaseQuals(variantReadBases.length());

        // double jitter for variants with long repeats
        String readBases = variantReadBases.substring(0, 16) + "A" + variantReadBases.substring(16);

        assertTrue(hasJitterMatchType(
                repeat, readContext, matcher.altIndexLower(), matcher.altIndexUpper(),
                varIndex + 1, readBases.getBytes(), readQuals, LENGTHENED, true, true, 1));

        readBases = variantReadBases.substring(0, 16) + "AA" + variantReadBases.substring(16);

        assertTrue(hasJitterMatchType(
                repeat, readContext, matcher.altIndexLower(), matcher.altIndexUpper(),
                varIndex + 2, readBases.getBytes(), readQuals, LENGTHENED, true, true, 2));

        readBases = variantReadBases.substring(0, 16) + variantReadBases.substring(18);

        assertTrue(hasJitterMatchType(
                repeat, readContext, matcher.altIndexLower(), matcher.altIndexUpper(),
                varIndex - 2, readBases.getBytes(), readQuals, SHORTENED, true, true, 2));
    }

    @Test
    public void testReadJitter()
    {
        SimpleVariant variant = createSimpleVariant(30, "A", "T");

        //                           10        20        30        40        50        60        70        80        90
        //                 0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
        String refBases = "ATGTCATTTTAGCGCGCATTCCTTCCTTCCAAAAAAAAAATGCTGGCACACACACATTAGTCGTAGATTAGTCGTAG";
        RefSequence refSequence = new RefSequence(0, refBases.getBytes());

        String readBases = refBases.substring(1, 30) + variant.alt() + refBases.substring(31, 71);
        String readCigar = "70M";
        SAMRecord read = buildSamRecord(1, readCigar, readBases);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        VariantReadContext readContext = builder.createContext(variant, read, 29, refSequence);
        assertTrue(readContext.isValid());
        assertEquals(23, readContext.VarIndex);
        assertEquals(2, readContext.AllRepeats.size());
        assertTrue(readContext.AllRepeats.stream().anyMatch(x -> x.matches("A", 24, 9)));
        assertTrue(readContext.AllRepeats.stream().anyMatch(x -> x.matches("TTCC", 11, 3)));

        ReadContextMatcher matcher = new ReadContextMatcher(readContext, true, true);

        SAMRecord read1 = buildSamRecord(1, readCigar, readBases);

        JitterMatch jitterMatch = checkJitter(readContext, matcher, read1, 29);
        assertEquals(SHORTENED, jitterMatch);
        // assertEquals(JitterMatch.NONE, jitterMatch); // TODO: matches since skips shortened bases

        readBases = refBases.substring(1, 30) + "TTCC" + variant.alt() + refBases.substring(31, 71);
        readCigar = buildCigarString(readBases.length());
        read1 = buildSamRecord(1, readCigar, readBases);

        jitterMatch = checkJitter(readContext, matcher, read1, 33);
        assertEquals(JitterMatch.LENGTHENED, jitterMatch);

        readBases = refBases.substring(1, 26) + variant.alt() + refBases.substring(31, 71);
        read1 = buildSamRecord(1, readCigar, readBases);

        jitterMatch = checkJitter(readContext, matcher, read1, 25);
        assertEquals(JitterMatch.SHORTENED, jitterMatch);


        // now test reads where the jitter is after the variant read index
        // 28 bases then the variant at the index = 29 then an extra 'A'
        readBases = refBases.substring(1, 30) + variant.alt() + "A" + refBases.substring(31, 71);
        readCigar = buildCigarString(readBases.length());
        read1 = buildSamRecord(1, readCigar, readBases);

        jitterMatch = checkJitter(readContext, matcher, read1, 29);
        assertEquals(JitterMatch.LENGTHENED, jitterMatch);

        readBases = refBases.substring(1, 30) + variant.alt() + refBases.substring(32, 71);
        read1 = buildSamRecord(1, readCigar, readBases);

        jitterMatch = checkJitter(readContext, matcher, read1, 29);
        assertEquals(JitterMatch.SHORTENED, jitterMatch);

        // test a partial core on each side - still valid for jitter

        // first missing the left core
        readBases = refBases.substring(13, 30) + variant.alt() + "A" + refBases.substring(31, 71);
        readCigar = buildCigarString(readBases.length());
        read1 = buildSamRecord(13, readCigar, readBases);

        jitterMatch = checkJitter(readContext, matcher, read1, 17);
        assertEquals(JitterMatch.LENGTHENED, jitterMatch);

        // now missing the right core
        readBases = refBases.substring(1, 30) + variant.alt() + refBases.substring(32, 46);
        readCigar = buildCigarString(readBases.length());
        read1 = buildSamRecord(1, readCigar, readBases);

        jitterMatch = checkJitter(readContext, matcher, read1, 29);
        assertEquals(JitterMatch.SHORTENED, jitterMatch);

        // test read not covering the repeat so cannot test it correctly
        readBases = refBases.substring(1, 30) + variant.alt() + refBases.substring(32, 40);
        readCigar = buildCigarString(readBases.length());
        read1 = buildSamRecord(1, readCigar, readBases);

        jitterMatch = checkJitter(readContext, matcher, read1, 29);
        assertEquals(JitterMatch.NONE, jitterMatch);
    }

    @Test
    public void testJitterSpanningVariantPos()
    {
        SimpleVariant variant = createSimpleVariant(30, "G", "A");

        //                           10        20        30        40        50        60        70
        //                 01234567890123456789012345678901234567890123456789012345678901234567890
        String refBases = "ATGTCATTTTAGCGCGCATTCTACACACACGCACACACACACACACACATTAGTCGTAGATTAGTCGTAGT";
        RefSequence refSequence = new RefSequence(0, refBases.getBytes());

        String readBases = refBases.substring(1, 30) + variant.alt() + refBases.substring(31, 71);
        String readCigar = buildCigarString(readBases.length());
        SAMRecord read = buildSamRecord(1, readCigar, readBases);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);
        VariantReadContext readContext = builder.createContext(variant, read, 29, refSequence);

        ReadContextMatcher matcher = createReadContextMatcher(readContext);

        // test 1: insert (ie lengthened) with a matching CIGAR
        String read1Bases = refBases.substring(1, 30) + "AC" + variant.alt() + refBases.substring(31, 71);
        SAMRecord readWithIndel = buildSamRecord(1, "21M2I49M", read1Bases);

        JitterMatch jitterMatch = checkJitter(readContext, matcher, readWithIndel, 31);
        assertEquals(JitterMatch.LENGTHENED, jitterMatch);

        // test 2: insert but without a matching CIGAR
        SAMRecord readWithoutIndel = buildSamRecord(1, buildCigarString(read1Bases.length()), read1Bases);
        jitterMatch = checkJitter(readContext, matcher, readWithoutIndel, 29);
        assertEquals(JitterMatch.LENGTHENED, jitterMatch);

        // test 3: deleted repeat bases matched in the CIGAR
        String read2Bases = refBases.substring(1, 28) + variant.alt() + refBases.substring(31, 71);
        readWithIndel = buildSamRecord(1, "21M2D47M", read2Bases);
        jitterMatch = checkJitter(readContext, matcher, readWithIndel, 27);
        assertEquals(JitterMatch.SHORTENED, jitterMatch);

        // test 4: deleted but all bases appear aligned
        readWithoutIndel = buildSamRecord(1, buildCigarString(read2Bases.length()), read2Bases);
        jitterMatch = checkJitter(readContext, matcher, readWithoutIndel, 29);
        assertEquals(JitterMatch.SHORTENED, jitterMatch);
    }

    @Test
    public void testJitterEndingOnVariantPos()
    {
        SimpleVariant variant = createSimpleVariant(49, "T", "C");

        //                           10        20        30        40        50        60        70
        //                 01234567890123456789012345678901234567890123456789012345678901234567890
        String refBases = "ATGTCATTTTAGCGCGCATTCTACACACACACACACACACACACACACATTAGTCGTAGATTAGTCGTAGT";
        RefSequence refSequence = new RefSequence(0, refBases.getBytes());

        String readBases = refBases.substring(1, 49) + variant.alt() + refBases.substring(50, 71);
        String readCigar = buildCigarString(readBases.length());
        SAMRecord read = buildSamRecord(1, readCigar, readBases);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);
        VariantReadContext readContext = builder.createContext(variant, read, 48, refSequence);

        ReadContextMatcher matcher = createReadContextMatcher(readContext);

        String read1Bases = refBases.substring(1, 49) + "CA" + variant.alt() + refBases.substring(50, 71);
        SAMRecord read1_withCigarElem = buildSamRecord(1, "21M2I49M", read1Bases);
        JitterMatch jitterMatch_withCigarElem = checkJitter(readContext, matcher, read1_withCigarElem, 50);
        assertEquals(JitterMatch.LENGTHENED, jitterMatch_withCigarElem);

        SAMRecord read1_noCigarElem = buildSamRecord(1, buildCigarString(read1Bases.length()), read1Bases);
        JitterMatch jitterMatch_noCigarElem = checkJitter(readContext, matcher, read1_noCigarElem, 48);
        assertEquals(JitterMatch.LENGTHENED, jitterMatch_noCigarElem);
    }

    @Test
    public void testMsiJitterCalcs()
    {
        JitterModelParams jitterParams = new JitterModelParams(
                "A/T",	0.05,0.05, 0.11,	0.0444,
                -0.1835,	1.0164);

        MsiModelParams modelParams = new MsiModelParams(jitterParams);

        double calc = modelParams.calcErrorRate(8, -1, null);
        assertEquals(0.002986, calc, 0.0001);

        calc = modelParams.calcErrorRate(6, 1, jitterParams.OptimalScaleRepeat6);
        assertEquals(0.00011, calc, 0.0001);
    }

    @Test
    public void testMsiJitterErrorRate()
    {
        MsiJitterCalcs msiJitterCalcs = new MsiJitterCalcs();

        JitterModelParams jitterParams1 = new JitterModelParams(
                "A/T",	0.05,0.06, 0.07,	0.0444,
                -0.1835,	1.0164);

        JitterModelParams jitterParams2 = new JitterModelParams(
                "AT/TA",	0.02,0.03, 0.04,	0.0444,
                -0.1835,	1.0164);

        JitterModelParams jitterParams3 = new JitterModelParams(
                REPEAT_UNIT_3_PLUS_LABEL,	0.1, 0.11, 0.12, 0.0133,
                0.1115, 1.6087);

        String sampleId = TEST_SAMPLE;
        msiJitterCalcs.setSampleParams(sampleId, List.of(jitterParams1, jitterParams2, jitterParams3));

        SimpleVariant variant = createSimpleVariant(100, "TA", "T");
        VariantReadContext readContext = createReadContext(variant, "GT", "AAAAT");

        double errorRate = msiJitterCalcs.calcErrorRate(readContext, sampleId);
        assertEquals(5.87e-8, errorRate, 1e-9);

        variant = createSimpleVariant(100, "T", "TAAAA");
        readContext = createReadContext(variant, "GT", "AAAAT");

        errorRate = msiJitterCalcs.calcErrorRate(readContext, sampleId);
        assertEquals(1.69e-35, errorRate, 1e-35);

        // too many repeats changing
        variant = createSimpleVariant(100, "T", "TAAAAAAA");
        readContext = createReadContext(variant, "GT", "AAAAT");

        errorRate = msiJitterCalcs.calcErrorRate(readContext, sampleId);
        assertEquals(0, errorRate, 0.001);

        // 4-base repeat
        variant = createSimpleVariant(100, "ACGTACGTA", "A");
        readContext = createReadContext(variant, "GT", "CGTACGTACGTACGTACGTGG");

        errorRate = msiJitterCalcs.calcErrorRate(readContext, sampleId);
        assertEquals(1.80e-7, errorRate, 1e-8);
    }

    @Test
    public void testJitterFilterQualState()
    {
        MsiJitterCalcs msiJitterCalcs = new MsiJitterCalcs();

        JitterModelParams jitterParams1 = new JitterModelParams(
                "A/C/G/T",	0.05,0.06, 0.07,	0.05,
                -0.02,	0.05);

        msiJitterCalcs.setSampleParams(TEST_SAMPLE, List.of(jitterParams1));

        SimpleVariant var = createSimpleVariant(100, "T", "TA");

        String readBases = REF_BASES_200.substring(80, 100) + var.alt() + "AAAAAA" + REF_BASES_200.substring(101, 121);
        byte[] baseQuals = buildDefaultBaseQuals(readBases.length());
        String readCigar = buildCigarString(readBases.length());
        SAMRecord read = buildSamRecord(30, readCigar, readBases, baseQuals);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        VariantReadContext readContext = builder.createContext(var, read, 20, REF_SEQUENCE_200);
        ReadContextCounter readContextCounter = createReadCounter(0, readContext);

        JitterData jitterData = new JitterData();

        // lengthened is noise
        jitterData.setValues(10, 10);
        readContextCounter.readSupportCounts().Full = 50;

        jitterData.setJitterQualFilterState(msiJitterCalcs, readContextCounter);

        assertEquals(1.2, jitterData.qualBoost(), 0.01);
        assertFalse(jitterData.filterOnNoise());

        // shortened is noise
        jitterData.setValues(5, 35);
        readContextCounter.readSupportCounts().Full = 50;

        jitterData.setJitterQualFilterState(msiJitterCalcs, readContextCounter);

        assertEquals(1.1, jitterData.qualBoost(), 0.01);
        assertFalse(jitterData.filterOnNoise());

        // both noise
        jitterData.setValues(5, 2);
        readContextCounter.readSupportCounts().Full = 50;

        jitterData.setJitterQualFilterState(msiJitterCalcs, readContextCounter);

        assertEquals(1.14, jitterData.qualBoost(), 0.01);
        assertFalse(jitterData.filterOnNoise());

        // shortened too high
        jitterData.setValues(21, 1);
        readContextCounter.readSupportCounts().Full = 10;

        jitterData.setJitterQualFilterState(msiJitterCalcs, readContextCounter);

        assertTrue(jitterData.filterOnNoise());
        assertEquals(1, jitterData.qualBoost(), 0.01);

        // again for lengthened
        jitterData.setValues(1, 60);
        readContextCounter.readSupportCounts().Full = 2;

        jitterData.setJitterQualFilterState(msiJitterCalcs, readContextCounter);

        assertTrue(jitterData.filterOnNoise());
        assertEquals(1, jitterData.qualBoost(), 0.01);

        // combined
        jitterData.setValues(126, 126);
        readContextCounter.readSupportCounts().Full = 35;

        jitterData.setJitterQualFilterState(msiJitterCalcs, readContextCounter);

        assertTrue(jitterData.filterOnNoise());
        assertEquals(1, jitterData.qualBoost(), 0.01);
    }
}
