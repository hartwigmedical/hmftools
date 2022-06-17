package com.hartwig.hmftools.neo.cohort;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.neo.bind.BindData;

public class NeoPredictionData
{
    public final int NeId;

    private final Map<String,List<BindData>> mPeptidePredictions; // keyed by allele

    public double MaxLikelihood;
    public double SumLikelihood;

    public NeoPredictionData(final int neId)
    {
        NeId = neId;

        mPeptidePredictions = Maps.newHashMap();

        MaxLikelihood = 0;
        SumLikelihood = 0;
    }

    public Map<String,List<BindData>> getPeptidePredictions() { return mPeptidePredictions; }

    public List<BindData> getPeptidePredictions(final String allele) { return mPeptidePredictions.get(allele); }
}
