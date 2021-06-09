package com.hartwig.hmftools.neo.predict;

public class HlaPeptidePrediction
{
    public final String HlaType;
    public final String Peptide;
    public final PredictionValues Prediction;

    public HlaPeptidePrediction(final String hlaType, final String peptide, final PredictionValues prediction)
    {
        HlaType = hlaType;
        Peptide = peptide;
        Prediction = prediction;
    }
}
