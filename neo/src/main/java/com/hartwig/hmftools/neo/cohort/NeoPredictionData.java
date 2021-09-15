package com.hartwig.hmftools.neo.cohort;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.neo.bind.BindData;

import org.apache.commons.compress.utils.Lists;

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

    /*
    public void calcSummaryValues()
    {
        for(BindData predData : mPeptidePredictions)
        {
            MaxLikelihood = max(predData.likelihood(), MaxLikelihood);
            SumLikelihood += predData.likelihood();
        }
    }
    */

}
