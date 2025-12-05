package com.hartwig.hmftools.qsee.feature;

public class Feature
{
    private final FeatureKey mKey;
    private final double mValue;

    public Feature(String name, double value, FeatureType type, SourceTool sourceTool)
    {
        mKey = new FeatureKey(name, type, sourceTool);
        mValue = value;
    }

    public Feature(FeatureKey key, double value)
    {
        mKey = key;
        mValue = value;
    }

    public FeatureKey key() { return mKey; }

    public double value() { return mValue; }
}
