package com.hartwig.hmftools.sig_analyser.cup;

import java.util.Map;

import com.google.common.collect.Maps;

public class CupSampleData
{
    public final String SampleId;
    public final String CancerType;
    public final String CancerSubtype;

    private int mIndex;
    private final Map<String,Double> mCancerCssTotals;

    public CupSampleData(final String sampleId, final String cancerType, final String cancerSubtype)
    {
        mIndex = -1;
        SampleId = sampleId;
        CancerType = cancerType;
        CancerSubtype = cancerSubtype;
        mCancerCssTotals = Maps.newHashMap();
    }

    public void setSampleIndex(final int sampleIndex) { mIndex = sampleIndex; }
    public int index() { return mIndex; }

    public void addSampleCss(final String cancerType, double cssWeight)
    {
        Double total = mCancerCssTotals.get(cancerType);

        if(total == null)
            mCancerCssTotals.put(cancerType, cssWeight);
        else
            mCancerCssTotals.put(cancerType, total + cssWeight);
    }

    public final Map<String,Double> getCancerCssTotals() { return mCancerCssTotals; }

    public double getTotalWeightedCss()
    {
        return mCancerCssTotals.values().stream().mapToDouble(x -> x).sum();
    }
}
