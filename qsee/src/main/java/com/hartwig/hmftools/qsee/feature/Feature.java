package com.hartwig.hmftools.qsee.feature;

public class Feature
{
    private final FeatureKey mKey;
    private final double mValue;
    private final PlotMetadata mPlotMetadata;

    public Feature(String name, double value, FeatureType type, SourceTool sourceTool, PlotMetadata plotMetadata)
    {
        mKey = new FeatureKey(name, type, sourceTool);
        mValue = value;
        mPlotMetadata = plotMetadata;
    }

    public Feature(FeatureKey key, double value, PlotMetadata plotMetadata)
    {
        mKey = key;
        mValue = value;
        mPlotMetadata = plotMetadata;
    }

    public Feature(FeatureKey key, double value)
    {
        mKey = key;
        mValue = value;
        mPlotMetadata = PlotMetadata.createEmpty();
    }

    public FeatureKey key() { return mKey; }
    public double value() { return mValue; }
    public PlotMetadata plotMetadata() { return mPlotMetadata; }
}
