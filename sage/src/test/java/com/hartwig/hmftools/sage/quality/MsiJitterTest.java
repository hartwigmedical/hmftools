package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContext;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.basequal.jitter.JitterModelParams;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.junit.Test;

public class MsiJitterTest
{
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
