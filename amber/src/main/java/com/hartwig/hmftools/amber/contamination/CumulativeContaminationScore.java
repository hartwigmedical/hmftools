package com.hartwig.hmftools.amber.contamination;

import java.util.List;

import com.hartwig.hmftools.amber.VafReading;

import org.apache.commons.math3.distribution.BinomialDistribution;

abstract class CumulativeContaminationScore implements SearchGrid.Calculator
{
    private final List<VafReading> ContaminationPoints;

    public CumulativeContaminationScore(final List<VafReading> contaminationPoints)
    {
        ContaminationPoints = contaminationPoints;
    }

    @Override
    public SearchGrid.ValueScore calculate(final double value)
    {
        double result = 0;
        for(VafReading point : ContaminationPoints)
        {
            result += cost(point, value);
        }
        return new SearchGrid.ValueScore(value, result);
    }

    abstract double cost(VafReading point, double value);
}
