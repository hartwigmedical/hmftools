package com.hartwig.hmftools.qsee.feature;

public class Feature
{
    private final FeatureKey mKey;
    private final double mValue;
    private final SourceTool mSourceTool;

    public Feature(FeatureType type, String name, double value, SourceTool sourceTool)
    {
        mKey = new FeatureKey(type, name);
        mValue = value;
        mSourceTool = sourceTool;
    }

    public Feature(FeatureKey key, double value, SourceTool sourceTool)
    {
        mKey = key;
        mValue = value;
        mSourceTool = sourceTool;
    }

    public FeatureKey key() { return mKey; }

    public double value() { return mValue; }

    public SourceTool sourceTool() { return mSourceTool; }
}
