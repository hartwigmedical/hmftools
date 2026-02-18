package com.hartwig.hmftools.qsee.feature;

import com.hartwig.hmftools.qsee.status.QcStatus;

import org.jetbrains.annotations.Nullable;

public class Feature
{
    private final FeatureKey mKey;
    private final double mValue;
    private final QcStatus mQcStatus;
    @Nullable private final String mPlotLabel;

    public Feature(String name, double value, FeatureType type, SourceTool sourceTool,
            @Nullable QcStatus qcStatus, @Nullable String plotLabel)
    {
        mKey = new FeatureKey(name, type, sourceTool);
        mValue = value;
        mQcStatus = qcStatus;
        mPlotLabel = plotLabel;
    }

    public Feature(FeatureKey key, double value, @Nullable QcStatus qcStatus, @Nullable String plotLabel)
    {
        mKey = key;
        mValue = value;
        mQcStatus = qcStatus;
        mPlotLabel = plotLabel;
    }

    public Feature(FeatureKey key, double value)
    {
        mKey = key;
        mValue = value;
        mQcStatus = QcStatus.createEmpty();
        mPlotLabel = null;
    }

    public FeatureKey key() { return mKey; }
    public double value() { return mValue; }
    public QcStatus qcStatus() { return mQcStatus; }
    public String plotLabel() { return mPlotLabel; }
}
