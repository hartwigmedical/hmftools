package com.hartwig.hmftools.amber.purity;

import java.util.List;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.amber.PositionEvidence;

public class VafLevelEvaluation implements Callable<VafLevelEvaluationResult>
{
    private static final int MINIMUM_CAPTURED_POINTS = 15;
    private final VafLevel Level;
    private final List<PositionEvidence> Evidence;
    private VafLevelEvaluationResult Result;

    public VafLevelEvaluation(final VafLevel level, final List<PositionEvidence> evidence)
    {
        Level = level;
        Evidence = evidence;
    }

    public VafLevelEvaluationResult call()
    {
        List<PositionEvidence> testable = Evidence.stream().filter(Level::hasSufficientDepthForEventDetection).toList();
        if(testable.size() < MINIMUM_CAPTURED_POINTS)
        {
            return VafLevelEvaluationResult.tooFewPointsOfSufficientDepth(Level, Level.numberOfCapturedEvidencePoints());
        }
        testable.forEach(Level::test);
        if(Level.numberOfCapturedEvidencePoints() < MINIMUM_CAPTURED_POINTS)
        {
            Result = VafLevelEvaluationResult.tooFewPointsCaptured(Level, Level.numberOfCapturedEvidencePoints());
        }
        else
        {
            double score = (double) Level.numberOfCapturedEvidencePoints() / testable.size();
            Result = new VafLevelEvaluationResult(Level, score, "");
        }
        return Result;
    }

    public boolean hasScore()
    {
        return Result != null && Result.Score != null;
    }

    public VafLevelEvaluationResult result()
    {
        return Result;
    }
}
