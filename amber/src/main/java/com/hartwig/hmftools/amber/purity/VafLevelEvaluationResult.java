package com.hartwig.hmftools.amber.purity;

public class VafLevelEvaluationResult implements Score
{
    public static VafLevelEvaluationResult tooFewPointsOfSufficientDepth(VafLevel vaf, int numberOfPoints)
    {
        return new VafLevelEvaluationResult(vaf, 0.0, "Too few points of sufficient depth to test contamination level: " + numberOfPoints);
    }

    public static VafLevelEvaluationResult tooFewPointsCaptured(VafLevel vaf, int numberOfPoints)
    {
        return new VafLevelEvaluationResult(vaf, 0.0, "Too few points in bands: " + numberOfPoints);
    }

    public final VafLevel Vaf;
    public final Double Score;
    public final String Message;

    public VafLevelEvaluationResult(final VafLevel vaf, final Double score, final String message)
    {
        Vaf = vaf;
        Score = score;
        Message = message;
    }

    @Override
    public double score()
    {
        return Score;
    }
}
