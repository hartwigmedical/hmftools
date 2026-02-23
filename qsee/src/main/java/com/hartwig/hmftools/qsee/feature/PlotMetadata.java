package com.hartwig.hmftools.qsee.feature;

import com.hartwig.hmftools.qsee.common.MultiFieldStringBuilder;
import com.hartwig.hmftools.qsee.status.QcStatus;

public class PlotMetadata
{
    private final String mFeatureGroup;
    private final String mPlotLabel;
    private final NumberFormat mNumberFormat;
    private final QcStatus mQcStatus;

    public static final String FLD_FEATURE_GROUP = "FeatureGroup";
    public static final String FLD_PLOT_LABEL = "PlotLabel";
    public static final String FLD_NUMBER_FORMAT = "NumberFormat";
    public static final String FLD_QC_STATUS = "QcStatus";

    public PlotMetadata(String featureGroup, String plotLabel, NumberFormat numberFormat, QcStatus qcStatus)
    {
        mFeatureGroup = featureGroup;
        mPlotLabel = plotLabel;
        mNumberFormat = numberFormat;
        mQcStatus = qcStatus;
    }

    public static Builder builder() { return new Builder(); }

    public static PlotMetadata createEmpty(){ return builder().build(); }

    public String featureGroup() { return mFeatureGroup; }
    public String plotLabel() { return mPlotLabel; }
    public NumberFormat numberFormat() { return mNumberFormat; }
    public QcStatus qcStatus() { return mQcStatus; }

    public String toString()
    {
        MultiFieldStringBuilder builder = new MultiFieldStringBuilder();

        builder.add(FLD_FEATURE_GROUP, mFeatureGroup);
        builder.add(FLD_PLOT_LABEL, mPlotLabel);
        builder.add(FLD_NUMBER_FORMAT, mNumberFormat.toString());
        builder.add(FLD_QC_STATUS, mQcStatus.toString());

        return builder.toString();
    }

    public static class Builder
    {
        private String mFeatureGroup = "";
        private String mPlotLabel = "";
        private NumberFormat mNumberFormat = NumberFormat.NUMBER;
        private QcStatus mQcStatus = QcStatus.createEmpty();

        public Builder featureGroup(String featureGroup)
        {
            mFeatureGroup = featureGroup;
            return this;
        }

        public Builder plotLabel(String plotLabel)
        {
            mPlotLabel = plotLabel;
            return this;
        }

        public Builder numberFormat(NumberFormat numberFormat)
        {
            mNumberFormat = numberFormat;
            return this;
        }

        public Builder qcStatus(QcStatus qcStatus)
        {
            mQcStatus = qcStatus;
            return this;
        }

        public PlotMetadata build()
        {
            return new PlotMetadata(mFeatureGroup, mPlotLabel, mNumberFormat, mQcStatus);
        }
    }
}
