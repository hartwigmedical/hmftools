package com.hartwig.hmftools.qsee.feature;

public class Feature
{
    private final FeatureKey mKey;
    private final double mValue;
    private final String mQcStatus;

    public Feature(String name, double value, FeatureType type, SourceTool sourceTool, String qcStatus)
    {
        mKey = new FeatureKey(name, type, sourceTool);
        mValue = value;
        mQcStatus = qcStatus;
    }

    public Feature(FeatureKey key, double value, String qcStatus)
    {
        mKey = key;
        mValue = value;
        mQcStatus = qcStatus;
    }

    public Feature(FeatureKey key, double value)
    {
        mKey = key;
        mValue = value;
        mQcStatus = "";
    }

    public FeatureKey key() { return mKey; }
    public double value() { return mValue; }
    public String qcStatus() { return mQcStatus; }
}
