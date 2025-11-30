package com.hartwig.hmftools.common.segmentation.copynumber;

public class GammaPenaltyCalculator implements PenaltyCalculator
{
    private final double UserGamma;
    private final boolean Normalise;

    public GammaPenaltyCalculator(final double userGamma, final boolean normalise)
    {
        UserGamma = userGamma;
        Normalise = normalise;
    }

    @Override
    public double getPenalty(double[] y)
    {
        return new Gamma(y, UserGamma, Normalise).getSegmentPenalty();
    }
}
