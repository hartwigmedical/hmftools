package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.amber.AmberConstants.MINIMUM_CAPTURED_POINTS;
import static com.hartwig.hmftools.amber.purity.CandidatePeakEvaluationResult.tooFewPointsOfSufficientDepth;

import java.util.List;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.amber.PositionEvidence;

public class CandidatePeakEvaluation implements Callable<Void>
{
    private final CandidatePeak mCandidatePeak;
    private final List<PositionEvidence> mEvidence;
    private CandidatePeakEvaluationResult mResult;

    public CandidatePeakEvaluation(final CandidatePeak candidatePeak, final List<PositionEvidence> evidence)
    {
        mCandidatePeak = candidatePeak;
        mEvidence = evidence;
    }

    public Void call()
    {
        List<PositionEvidence> testable = mEvidence.stream().filter(mCandidatePeak::hasSufficientDepthForEventDetection).toList();
        if(testable.size() < MINIMUM_CAPTURED_POINTS)
        {
            mResult = tooFewPointsOfSufficientDepth(mCandidatePeak, mCandidatePeak.numberOfCapturedEvidencePoints());
            return null;
        }
        testable.forEach(mCandidatePeak::test);
        if(mCandidatePeak.numberOfCapturedEvidencePoints() < MINIMUM_CAPTURED_POINTS)
        {
            mResult = CandidatePeakEvaluationResult.tooFewPointsCaptured(mCandidatePeak, mCandidatePeak.numberOfCapturedEvidencePoints());
        }
        else
        {
            double score = (double) mCandidatePeak.numberOfCapturedEvidencePoints() / testable.size();
            mResult = new CandidatePeakEvaluationResult(mCandidatePeak, score, "");
        }
        return null;
    }

    public CandidatePeakEvaluationResult result()
    {
        return mResult;
    }
}
