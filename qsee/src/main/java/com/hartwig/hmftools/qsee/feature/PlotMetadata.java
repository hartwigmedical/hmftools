package com.hartwig.hmftools.qsee.feature;

import com.hartwig.hmftools.qsee.common.MultiFieldStringBuilder;
import com.hartwig.hmftools.qsee.status.QcStatus;
import com.hartwig.hmftools.qsee.status.QcStatusType;

public class PlotMetadata
{
    private final String mFeatureGroup;
    private final String mPlotLabel;
    private final NumberFormat mNumberFormat;
    private final QcStatus mQcStatus;

    public static final String FIELD_FEATURE_GROUP = "FeatureGroup";
    public static final String FIELD_PLOT_LABEL = "PlotLabel";
    public static final String FIELD_NUMBER_FORMAT = "NumberFormat";
    public static final String FIELD_QC_STATUS = "QcStatus";
    public static final String FIELD_QC_THRESHOLD = "QcThreshold";

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

        builder.add(FIELD_FEATURE_GROUP, mFeatureGroup);
        builder.add(FIELD_PLOT_LABEL, mPlotLabel);
        builder.add(FIELD_NUMBER_FORMAT, getNumberFormatString());
        builder.add(FIELD_QC_STATUS, getQcStatusString());
        builder.add(FIELD_QC_THRESHOLD, getQcThresholdString());

        return builder.toString();
    }

    private String getNumberFormatString()
    {
        return mNumberFormat == null ? "" : mNumberFormat.toString();
    }

    private String getQcStatusString()
    {
        return mQcStatus.type() == QcStatusType.NONE ? "" : mQcStatus.type().toString();
    }

    private String getQcThresholdString()
    {
        if(mQcStatus.type() == QcStatusType.NONE)
            return "";

        double thresholdValue = mQcStatus.threshold();
        boolean isPercent = mNumberFormat == NumberFormat.PERCENT;

        if(isPercent)
            thresholdValue = thresholdValue * 100;

        boolean isInteger = thresholdValue % 1 == 0;

        String thresholdString = isInteger
                ? String.valueOf((int) thresholdValue)
                : String.valueOf(thresholdValue);

        if(isPercent)
            thresholdString = thresholdString + "%";

        return mQcStatus.operator().operatorString() + thresholdString;
    }

    public static class Builder
    {
        private String mFeatureGroup = "";
        private String mPlotLabel = "";
        private NumberFormat mNumberFormat = null;
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
