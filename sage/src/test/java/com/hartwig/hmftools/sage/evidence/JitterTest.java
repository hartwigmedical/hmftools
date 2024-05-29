package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContext;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;

import static org.junit.Assert.assertEquals;
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

        double calc = modelParams.calcSkew(8, -1);
        assertEquals(0.002986, calc, 0.0001);
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
                "3-5bp repeat",	0.1, 0.11, 0.12, 0.0133,
                0.1115, 1.6087);

        String sampleId = "SAMPLE";
        msiJitterCalcs.setSampleParams("SAMPLE", List.of(jitterParams1, jitterParams2, jitterParams3));

        SimpleVariant variant = createSimpleVariant(100, "AA", "A");
        VariantReadContext readContext = createReadContext(variant, "GT", "AAAAT");

        double errorRate = msiJitterCalcs.calcErrorRate(readContext, sampleId);
        assertEquals(0.06, errorRate, 0.001);

        variant = createSimpleVariant(100, "A", "AAAA");
        readContext = createReadContext(variant, "GT", "AAAAT");

        errorRate = msiJitterCalcs.calcErrorRate(readContext, sampleId);
        assertEquals(0.05, errorRate, 0.001);

        // too many repeats changing
        variant = createSimpleVariant(100, "A", "AAAAAAA");
        readContext = createReadContext(variant, "GT", "AAAAT");

        errorRate = msiJitterCalcs.calcErrorRate(readContext, sampleId);
        assertEquals(0, errorRate, 0.001);

        // 4-base repeat
        variant = createSimpleVariant(100, "ACGTACGTA", "A");
        readContext = createReadContext(variant, "GT", "CGTACGTACGTACGTACGTGG");

        errorRate = msiJitterCalcs.calcErrorRate(readContext, sampleId);
        assertEquals(0.12, errorRate, 0.001);
    }
}
