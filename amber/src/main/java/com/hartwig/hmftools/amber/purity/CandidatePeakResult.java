package com.hartwig.hmftools.amber.purity;

public record CandidatePeakResult(CandidatePeak candidatePeak, double score, String message) implements Score
{
    public static CandidatePeakResult tooFewPointsOfSufficientDepth(CandidatePeak vaf, int numberOfPoints)
    {
        return new CandidatePeakResult(vaf, 0.0,
                "Too few points of sufficient depth to test contamination level: " + numberOfPoints);
    }

    public static CandidatePeakResult tooFewPointsCaptured(CandidatePeak vaf, int numberOfPoints)
    {
        return new CandidatePeakResult(vaf, 0.0, "Too few points in bands: " + numberOfPoints);
    }
}
