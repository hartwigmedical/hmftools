package com.hartwig.hmftools.amber.purity;

import org.junit.Assert;
import org.junit.Test;

public class RangeStepPeakTest
{
    @Test
    public void exampleTest()
    {
        int refSupport = 1739;
        int altSupport = 7;
        int depth = refSupport + altSupport;
        double symmetricVaf = (double) altSupport / depth;
        CandidatePeak.RangeStepPeak homPeak = new CandidatePeak.RangeStepPeak(depth, altSupport, 0.006, 0.001103);
        Assert.assertTrue(homPeak.isCaptured(symmetricVaf));
        CandidatePeak.RangeStepPeak hetPeak = new CandidatePeak.RangeStepPeak(depth, altSupport, 0.003, 0.001103 / 2.0);
        Assert.assertFalse(hetPeak.isCaptured(symmetricVaf));
    }
}
