package com.hartwig.hmftools.neo.cohort;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;

public class NeoPredictionData
{
    public final int NeId;

    public int Peptides;
    public double MaxLikelihood;
    public double SumLikelihood;

    public NeoPredictionData(final int neId)
    {
        NeId = neId;

        Peptides = 0;
        MaxLikelihood = 0;
        SumLikelihood = 0;
    }

    public void processPredictionData(final BindingPredictionData predData)
    {
        ++Peptides;
        MaxLikelihood = max(predData.likelihood(), MaxLikelihood);
        SumLikelihood += predData.likelihood();
    }

}
