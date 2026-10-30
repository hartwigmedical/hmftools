package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.amber.purity.CandidatePeakResult.tooFewPointsCaptured;
import static com.hartwig.hmftools.amber.purity.CandidatePeakResult.tooFewPointsOfSufficientDepth;

import org.junit.Assert;
import org.junit.Test;

public class CandidatePeakEvaluationResultTest
{
    @Test
    public void tooFewPointsCapturedTest()
    {
        CandidatePeakResult result = tooFewPointsCaptured(new CandidatePeak(0.1), 14);
        Assert.assertEquals(0.0, result.score(), 0.0001);
    }

    @Test
    public void tooFewPointsOfSufficientDepthTest()
    {
        CandidatePeakResult result = tooFewPointsOfSufficientDepth(new CandidatePeak(0.1), 14);
        Assert.assertEquals(0.0, result.score(), 0.0001);
    }
}
