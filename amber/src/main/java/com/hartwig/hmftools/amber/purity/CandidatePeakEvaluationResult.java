package com.hartwig.hmftools.amber.purity;

public record CandidatePeakEvaluationResult(CandidatePeak candidatePeak, double score, String message) implements Score
{
    public static CandidatePeakEvaluationResult tooFewPointsOfSufficientDepth(CandidatePeak vaf, int numberOfPoints)
    {
        return new CandidatePeakEvaluationResult(vaf, 0.0,
                "Too few points of sufficient depth to test contamination level: " + numberOfPoints);
    }

    public static CandidatePeakEvaluationResult tooFewPointsCaptured(CandidatePeak vaf, int numberOfPoints)
    {
        return new CandidatePeakEvaluationResult(vaf, 0.0, "Too few points in bands: " + numberOfPoints);
    }
}
