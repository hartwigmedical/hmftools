package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.amber.purity.CandidatePeakEvaluationResult.tooFewPointsCaptured;
import static com.hartwig.hmftools.amber.purity.CandidatePeakEvaluationResult.tooFewPointsOfSufficientDepth;

import org.junit.Assert;
import org.junit.Test;

public class CandidatePeakEvaluationResultTest
{
    @Test
    public void tooFewPointsCapturedTest()
    {
        CandidatePeakEvaluationResult result = tooFewPointsCaptured(new CandidatePeak(0.1), 14);
        Assert.assertEquals(0.0, result.score(), 0.0001);
    }

    @Test
    public void tooFewPointsOfSufficientDepthTest()
    {
        CandidatePeakEvaluationResult result = tooFewPointsOfSufficientDepth(new CandidatePeak(0.1), 14);
        Assert.assertEquals(0.0, result.score(), 0.0001);
    }
}
