package com.hartwig.hmftools.cup.common;

import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_UNKNOWN;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.cup.drivers.SampleDriverData;

import org.apache.commons.compress.utils.Lists;

public class SampleData
{
    public final String Id;
    public final String CancerType;
    public final String CancerSubtype;

    // SNV signature data and results
    private int mSnvCountsIndex;
    private int mSnvCount;
    private final Map<String,Double> mCancerSnvCssTotals;
    private final Map<String,Double> mSnvSigAllocations;

    // drive genes and results
    private final List<SampleDriverData> mDrivers;
    private final Map<String,Double> mDriverCancerTypeProbs;

    public SampleData(final String id, final String cancerType, final String cancerSubtype)
    {
        Id = id;
        CancerType = cancerType;
        CancerSubtype = cancerSubtype;

        mSnvCountsIndex = -1;
        mSnvCount = 0;
        mCancerSnvCssTotals = Maps.newHashMap();
        mSnvSigAllocations = Maps.newHashMap();
        mDriverCancerTypeProbs = Maps.newHashMap();
        mDrivers = Lists.newArrayList();
    }

    public boolean isUnknownCancerType() { return CancerType.equalsIgnoreCase(CANCER_TYPE_UNKNOWN); }

    public void addSampleCss(final String cancerType, double cssWeight)
    {
        Double total = mCancerSnvCssTotals.get(cancerType);

        if(total == null)
            mCancerSnvCssTotals.put(cancerType, cssWeight);
        else
            mCancerSnvCssTotals.put(cancerType, total + cssWeight);
    }

    // signatures
    public void setSnvCountsIndex(int sampleIndex) { mSnvCountsIndex = sampleIndex; }
    public void setSnvCount(int snvCount) { mSnvCount = snvCount; }
    public int getSnvCount() { return mSnvCount; }
    public int getSnvIndex() { return mSnvCountsIndex; }
    public final Map<String,Double> getCancerCssTotals() { return mCancerSnvCssTotals; }
    public final Map<String,Double> getSnvSigAllocations() { return mSnvSigAllocations; }

    public double getTotalWeightedCss()
    {
        return mCancerSnvCssTotals.values().stream().mapToDouble(x -> x).sum();
    }

    // drivers
    public final List<SampleDriverData> getDrivers() { return mDrivers; }
    public final Map<String,Double> getDriverCancerTypeProbs() { return mDriverCancerTypeProbs; }
}
