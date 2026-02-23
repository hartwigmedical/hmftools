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
        builder.add(FIELD_NUMBER_FORMAT, mNumberFormat.toString());
        builder.add(FIELD_QC_STATUS, getFormattedQcStatus());

        return builder.toString();
    }

    private String getFormattedQcStatus()
    {
        QcStatusType qcStatusType = mQcStatus.type();
        if(qcStatusType == QcStatusType.NONE)
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

        return qcStatusType + " " + mQcStatus.operator().operatorString() + thresholdString;
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
