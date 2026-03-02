package com.hartwig.hmftools.amber.contamination;

import java.util.List;

import com.hartwig.hmftools.amber.VafReading;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class CumulativeProbabilityContaminationCost extends CumulativeContaminationScore
{
    public CumulativeProbabilityContaminationCost(final List<VafReading> contaminationPoints)
    {
        super(contaminationPoints);
    }

    @Override
    double cost(VafReading point, double value)
    {
        int support = point.altSupport();
        final int depth = point.readDepth();
        if(support > depth / 2)
        {
            support = depth - support;
        }
        BinomialDistribution hetDist = new BinomialDistribution(depth, value * 0.5);
        double hetScore = Math.abs(0.5 - hetDist.cumulativeProbability(support));
        BinomialDistribution homDist = new BinomialDistribution(depth, value);
        double homScore = Math.abs(0.5 - homDist.cumulativeProbability(support));
        return Math.min(hetScore, homScore);
    }
}
