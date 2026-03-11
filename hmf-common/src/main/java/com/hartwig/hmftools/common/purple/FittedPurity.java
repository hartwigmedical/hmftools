package com.hartwig.hmftools.common.purple;

public record FittedPurity(double purity, double normFactor, double ploidy, double score,
                           double diploidProportion, double somaticPenalty) implements Comparable<FittedPurity>
{
    @Override
    public int compareTo(final FittedPurity o)
    {
        return Double.compare(score(), o.score());
    }
}
