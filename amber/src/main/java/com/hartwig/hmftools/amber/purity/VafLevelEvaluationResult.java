package com.hartwig.hmftools.amber.purity;

public record VafLevelEvaluationResult(VafLevel Vaf, Double Score, String Message) implements Score
{
    public static VafLevelEvaluationResult tooFewPointsOfSufficientDepth(VafLevel vaf, int numberOfPoints)
    {
        return new VafLevelEvaluationResult(vaf, 0.0, "Too few points of sufficient depth to test contamination level: " + numberOfPoints);
    }

    public static VafLevelEvaluationResult tooFewPointsCaptured(VafLevel vaf, int numberOfPoints)
    {
        return new VafLevelEvaluationResult(vaf, 0.0, "Too few points in bands: " + numberOfPoints);
    }

    @Override
    public double score()
    {
        return Score;
    }
}
