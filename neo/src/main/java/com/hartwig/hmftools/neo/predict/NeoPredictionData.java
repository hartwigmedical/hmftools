package com.hartwig.hmftools.neo.predict;

import java.util.List;

import com.google.common.collect.Lists;

public class NeoPredictionData
{
    public final int NeId;

    public final List<HlaPeptidePrediction> Predictions;

    public NeoPredictionData(final int neId)
    {
        NeId = neId;
        Predictions = Lists.newArrayList();
    }

    // NeId,HlaAllele,Peptide,affinity,affinity_percentile,processing_score,presentation_score,presentation_percentile
    //0,HLA-A2410,VEGCPSGTF,21473.741935853723,4.006624999999998,0.06819962710142136,0.005919105212229773,46.224673913043475
    //0,HLA-A2410,LLAAVEGCPS,28958.924528141455,18.636375,0.0038219578564167023,0.003533626812025815,99.28660326086957
}
