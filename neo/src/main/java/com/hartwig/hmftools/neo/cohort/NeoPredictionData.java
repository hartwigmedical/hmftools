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
    public double MinPresentationPerc;
    public double MinAffinity;
    public double SumPresentation;
    public double SumAffinity;

    public NeoPredictionData(final int neId)
    {
        NeId = neId;

        Peptides = 0;
        MaxLikelihood = 0;
        SumLikelihood = 0;
        MinPresentationPerc = -1;
        MinAffinity = -1;
        SumPresentation = 0;
        SumAffinity = 0;
    }

    public void processPredictionData(final BindingPredictionData predData, double mcfFactor)
    {
        ++Peptides;
        MaxLikelihood = max(predData.likelihood(), MaxLikelihood);
        SumLikelihood += predData.likelihood();

        if(MinAffinity < 0)
            MinAffinity = predData.affinity();
        else
            MinAffinity = min(MinAffinity, predData.affinity());

        if(MinPresentationPerc < 0)
            MinPresentationPerc = predData.presentationPerc();
        else
            MinPresentationPerc = min(MinPresentationPerc, predData.presentationPerc());

        SumPresentation += pow(predData.presentation(), mcfFactor);
        SumAffinity += 1.0 / pow(predData.affinity(), mcfFactor);
    }

}
