package com.hartwig.hmftools.qsee.feature;

import org.jetbrains.annotations.Nullable;

public class Feature
{
    private final FeatureKey mKey;
    private final double mValue;
    @Nullable private final QcStatus mQcStatus;

    public Feature(String name, double value, FeatureType type, SourceTool sourceTool, @Nullable QcStatus qcStatus)
    {
        mKey = new FeatureKey(name, type, sourceTool);
        mValue = value;
        mQcStatus = qcStatus;
    }

    public Feature(FeatureKey key, double value, @Nullable QcStatus qcStatus)
    {
        mKey = key;
        mValue = value;
        mQcStatus = qcStatus;
    }

    public Feature(FeatureKey key, double value)
    {
        mKey = key;
        mValue = value;
        mQcStatus = null;
    }

    public FeatureKey key() { return mKey; }
    public double value() { return mValue; }
    public QcStatus qcStatus() { return mQcStatus; }
}
