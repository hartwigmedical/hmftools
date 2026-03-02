package com.hartwig.hmftools.amber.contamination;

import java.util.List;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.amber.VafReading;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class LogLikelihoodContaminationCost extends CumulativeContaminationScore
{
    public LogLikelihoodContaminationCost(final List<VafReading> contaminationPoints)
    {
        super(contaminationPoints);
    }

    @Override
    double cost(VafReading point, double contaminationLevel)
    {
        Preconditions.checkArgument(contaminationLevel > 0 && contaminationLevel < 1, "Contamination level must be in (0,1)");
        int n = point.readDepth();
        int k = point.altSupport();
        if(k > n / 2)
        {
            k = n - k;
        }
        double f = point.minorAlleleFrequency();
        double minimumProbability = 0.0001;

        double probabilityWhenContaminationIsHomozygousRef = new BinomialDistribution(n, minimumProbability).probability(k);
        double probabilityWhenContaminationIsHeterozygous = new BinomialDistribution(n, contaminationLevel * 0.5).probability(k);
        double probabilityWhenContaminationIsHomozygousAlt = new BinomialDistribution(n, contaminationLevel).probability(k);

        double rawCost = (1 - f) * (1 - f) * probabilityWhenContaminationIsHomozygousRef
                + 2 * f * (1 - f) * probabilityWhenContaminationIsHeterozygous
                + f * f * probabilityWhenContaminationIsHomozygousAlt;

        return Math.log(Math.max(rawCost, Double.MIN_VALUE));
    }
}
