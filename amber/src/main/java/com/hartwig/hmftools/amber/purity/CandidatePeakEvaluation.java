package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.amber.purity.CandidatePeakEvaluationResult.tooFewPointsOfSufficientDepth;

import java.util.List;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.amber.PositionEvidence;

public class CandidatePeakEvaluation implements Callable<Void>
{
    private static final int MINIMUM_CAPTURED_POINTS = 15;
    private final CandidatePeak CandidatePeak;
    private final List<PositionEvidence> Evidence;
    private CandidatePeakEvaluationResult Result;

    public CandidatePeakEvaluation(final CandidatePeak candidatePeak, final List<PositionEvidence> evidence)
    {
        CandidatePeak = candidatePeak;
        Evidence = evidence;
    }

    public Void call()
    {
        List<PositionEvidence> testable = Evidence.stream().filter(CandidatePeak::hasSufficientDepthForEventDetection).toList();
        if(testable.size() < MINIMUM_CAPTURED_POINTS)
        {
            Result = tooFewPointsOfSufficientDepth(CandidatePeak, CandidatePeak.numberOfCapturedEvidencePoints());
            return null;
        }
        testable.forEach(CandidatePeak::test);
        if(CandidatePeak.numberOfCapturedEvidencePoints() < MINIMUM_CAPTURED_POINTS)
        {
            Result = CandidatePeakEvaluationResult.tooFewPointsCaptured(CandidatePeak, CandidatePeak.numberOfCapturedEvidencePoints());
        }
        else
        {
            double score = (double) CandidatePeak.numberOfCapturedEvidencePoints() / testable.size();
            Result = new CandidatePeakEvaluationResult(CandidatePeak, score, "");
        }
        return null;
    }

    public CandidatePeakEvaluationResult result()
    {
        return Result;
    }
}
