package com.hartwig.hmftools.qsee.feature;

public class Feature
{
    private final FeatureKey mKey;
    private final double mValue;
    private final FeatureMetadata mFeatureMetadata;

    public Feature(String name, double value, FeatureType type, SourceTool sourceTool, FeatureMetadata featureMetadata)
    {
        mKey = new FeatureKey(name, type, sourceTool);
        mValue = value;
        mFeatureMetadata = featureMetadata;
    }

    public Feature(String name, double value, FeatureType type, SourceTool sourceTool)
    {
        mKey = new FeatureKey(name, type, sourceTool);
        mValue = value;
        mFeatureMetadata = FeatureMetadata.createEmpty();
    }

    public Feature(FeatureKey key, double value, FeatureMetadata featureMetadata)
    {
        mKey = key;
        mValue = value;
        mFeatureMetadata = featureMetadata;
    }

    public Feature(FeatureKey key, double value)
    {
        mKey = key;
        mValue = value;
        mFeatureMetadata = FeatureMetadata.createEmpty();
    }

    public FeatureKey key() { return mKey; }
    public double value() { return mValue; }
    public FeatureMetadata metadata() { return mFeatureMetadata; }
}
