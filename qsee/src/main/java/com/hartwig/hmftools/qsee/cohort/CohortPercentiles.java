package com.hartwig.hmftools.qsee.cohort;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureKey;

public class CohortPercentiles
{
    private final Map<SampleType, Map<FeatureKey, FeaturePercentiles>> mCohortData = new LinkedHashMap<>();

    public Map<SampleType, Map<FeatureKey, FeaturePercentiles>> getData() { return mCohortData; }

    public void add(FeaturePercentiles featurePercentiles)
    {
        mCohortData.computeIfAbsent(featurePercentiles.sampleType(), k -> new LinkedHashMap<>());

        mCohortData
                .get(featurePercentiles.sampleType())
                .put(featurePercentiles.featureKey(), featurePercentiles);
    }

    public void addAll(List<FeaturePercentiles> featurePercentilesList)
    {
        for(FeaturePercentiles featurePercentiles : featurePercentilesList)
        {
            add(featurePercentiles);
        }
    }

    public FeaturePercentiles getFeaturePercentiles(SampleType sampleType, FeatureKey featureKey)
    {
        return mCohortData.get(sampleType).get(featureKey);
    }
}
