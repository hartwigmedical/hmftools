package com.hartwig.hmftools.common.segmentation;

public class FixedPenalty implements PenaltyCalculator
{
    private final double Penalty;

    public FixedPenalty(final double penalty)
    {
        Penalty = penalty;
    }

    @Override
    public double getPenalty(double[] y)
    {
        return Penalty;
    }
}
