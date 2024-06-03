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
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadCounter;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.basequal.jitter.JitterModelParams;
import com.hartwig.hmftools.sage.common.RefSequence;
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
        assertEquals(3, readContext.AllRepeats.size());

        SAMRecord read1 = buildSamRecord(1, readCigar, readBases);

        JitterMatch jitterMatch = JitterData.checkJitter(readContext, read1, 29);
        assertEquals(JitterMatch.NONE, jitterMatch);

        readBases = refBases.substring(1, 30) + "TTCC" + variant.alt() + refBases.substring(31, 71);
        read1 = buildSamRecord(1, readCigar, readBases);

        jitterMatch = JitterData.checkJitter(readContext, read1, 33);
        assertEquals(JitterMatch.LENGTHENED, jitterMatch);

        readBases = refBases.substring(1, 26) + variant.alt() + refBases.substring(31, 71);
        read1 = buildSamRecord(1, readCigar, readBases);

        jitterMatch = JitterData.checkJitter(readContext, read1, 25);
        assertEquals(JitterMatch.SHORTENED, jitterMatch);

        readBases = refBases.substring(1, 30) + variant.alt() + "A" + refBases.substring(31, 71);
        read1 = buildSamRecord(1, readCigar, readBases);

        jitterMatch = JitterData.checkJitter(readContext, read1, 30);
        assertEquals(JitterMatch.LENGTHENED, jitterMatch);

        readBases = refBases.substring(1, 30) + variant.alt() + refBases.substring(32, 71);
        read1 = buildSamRecord(1, readCigar, readBases);

        jitterMatch = JitterData.checkJitter(readContext, read1, 28);
        assertEquals(JitterMatch.SHORTENED, jitterMatch);
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
                REPEAT_UNIT_3_PLUS_LABEL,	0.05,0.06, 0.07,	0.05,
                -0.02,	0.05);

        msiJitterCalcs.setSampleParams(TEST_SAMPLE, List.of(jitterParams1));

        SimpleVariant var = createSimpleVariant(100, "A", "AA");

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
        assertFalse(jitterData.isWithinNoise());

        // shortened is noise
        jitterData.setValues(5, 35);
        readContextCounter.readSupportCounts().Full = 50;

        jitterData.setJitterQualFilterState(msiJitterCalcs, readContextCounter);

        assertEquals(1.1, jitterData.qualBoost(), 0.01);
        assertFalse(jitterData.isWithinNoise());

        // both noise
        jitterData.setValues(5, 2);
        readContextCounter.readSupportCounts().Full = 50;

        jitterData.setJitterQualFilterState(msiJitterCalcs, readContextCounter);

        assertEquals(1.14, jitterData.qualBoost(), 0.01);
        assertFalse(jitterData.isWithinNoise());

        // shortened too high
        jitterData.setValues(21, 1);
        readContextCounter.readSupportCounts().Full = 10;

        jitterData.setJitterQualFilterState(msiJitterCalcs, readContextCounter);

        assertTrue(jitterData.isWithinNoise());
        assertEquals(0, jitterData.qualBoost(), 0.01);

        // again for lengthened
        jitterData.setValues(1, 60);
        readContextCounter.readSupportCounts().Full = 2;

        jitterData.setJitterQualFilterState(msiJitterCalcs, readContextCounter);

        assertTrue(jitterData.isWithinNoise());
        assertEquals(0, jitterData.qualBoost(), 0.01);

        // combined
        jitterData.setValues(125, 125);
        readContextCounter.readSupportCounts().Full = 50;

        jitterData.setJitterQualFilterState(msiJitterCalcs, readContextCounter);

        assertTrue(jitterData.isWithinNoise());
        assertEquals(0, jitterData.qualBoost(), 0.01);
    }
}
