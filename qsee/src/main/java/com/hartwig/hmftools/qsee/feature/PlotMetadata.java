package com.hartwig.hmftools.qsee.feature;

import com.hartwig.hmftools.qsee.common.MultiFieldStringBuilder;
import com.hartwig.hmftools.qsee.status.QcStatus;

public class PlotMetadata
{
    private final String mFeatureGroup;
    private final String mDisplayName;
    private final NumberFormat mNumberFormat;
    private final QcStatus mQcStatus;

    public static final String FIELD_FEATURE_GROUP = "FeatureGroup";
    public static final String FIELD_DISPLAY_NAME = "DisplayName";
    public static final String FIELD_NUMBER_FORMAT = "NumberFormat";
    public static final String FIELD_QC_STATUS = "QcStatus";
    public static final String FIELD_QC_THRESHOLD = "QcThreshold";

    public PlotMetadata(String featureGroup, String displayName, NumberFormat numberFormat, QcStatus qcStatus)
    {
        mFeatureGroup = featureGroup;
        mDisplayName = displayName;
        mNumberFormat = numberFormat;
        mQcStatus = qcStatus;
    }

    public static Builder builder() { return new Builder(); }

    public static PlotMetadata createEmpty(){ return builder().build(); }

    public String featureGroup() { return mFeatureGroup; }
    public String displayName() { return mDisplayName; }
    public NumberFormat numberFormat() { return mNumberFormat; }
    public QcStatus qcStatus() { return mQcStatus; }

    public String displayString()
    {
        MultiFieldStringBuilder builder = new MultiFieldStringBuilder();

        builder.add(FIELD_FEATURE_GROUP, mFeatureGroup);
        builder.add(FIELD_DISPLAY_NAME, mDisplayName);
        builder.add(FIELD_NUMBER_FORMAT, mNumberFormat.displayString());
        builder.add(FIELD_QC_STATUS, mQcStatus.type().displayString());
        builder.add(FIELD_QC_THRESHOLD, mQcStatus.displayString(mNumberFormat));

        return builder.toString();
    }

    @Override
    public String toString(){ return displayString(); }

    public static class Builder
    {
        private String mFeatureGroup = "";
        private String mDisplayName = "";
        private NumberFormat mNumberFormat = NumberFormat.NONE;
        private QcStatus mQcStatus = QcStatus.createEmpty();

        public Builder featureGroup(String featureGroup)
        {
            mFeatureGroup = featureGroup;
            return this;
        }

        public Builder displayName(String displayName)
        {
            mDisplayName = displayName;
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
            return new PlotMetadata(mFeatureGroup, mDisplayName, mNumberFormat, mQcStatus);
        }
    }
}
