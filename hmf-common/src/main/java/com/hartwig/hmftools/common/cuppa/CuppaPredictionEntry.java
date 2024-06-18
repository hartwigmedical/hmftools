package com.hartwig.hmftools.common.cuppa;

import java.util.List;

public class CuppaPredictionEntry
{
    protected static final String FLD_SAMPLE_ID = "sample_id";
    public static final String FLD_DATA_TYPE = "data_type";
    public static final String FLD_CLASSIFIER_GROUP = "clf_group"; // Potentially unused
    public static final String FLD_CLASSIFIER_NAME = "clf_name";
    public static final String FLD_FEATURE_NAME = "feat_name";
    public static final String FLD_FEATURE_VALUE = "feat_value";
    protected static final String FLD_CANCER_TYPE = "cancer_type";
    public static final String FLD_DATA_VALUE = "data_value";
    public static final String FLD_RANK = "rank";
    public static final String FLD_RANK_GROUP = "rank_group";

    public final String SampleId;
    public final DataType DataType;
    public final ClassifierGroup ClassifierGroup;
    public final ClassifierName ClassifierName;
    public final String FeatureName;
    public final double FeatureValue;
    public final String CancerType;
    public final double DataValue;
    public final int Rank;
    public final int RankGroup;

    public CuppaPredictionEntry(
            final String sampleId,
            final DataType dataType,
            final ClassifierGroup classifierGroup,
            final ClassifierName classifierName,
            final String featureName,
            final double featureValue,
            final String cancerType,
            final double dataValue,
            final int rank,
            final int rankGroup
    )
    {
        SampleId = sampleId;
        DataType = dataType;
        ClassifierGroup = classifierGroup;
        ClassifierName = classifierName;
        FeatureName = featureName;
        FeatureValue = featureValue;
        CancerType = cancerType;
        DataValue = dataValue;
        Rank = rank;
        RankGroup = rankGroup;
    }

    public List<String> getFieldNames()
    {
        return List.of(
                FLD_SAMPLE_ID,
                FLD_DATA_TYPE,
                FLD_CLASSIFIER_GROUP,
                FLD_CLASSIFIER_NAME,
                FLD_FEATURE_NAME,
                FLD_FEATURE_VALUE,
                FLD_CANCER_TYPE,
                FLD_DATA_VALUE,
                FLD_RANK,
                FLD_RANK_GROUP
        );
    }

    public List<Object> getValues()
    {
        return List.of(
                SampleId,
                DataType,
                ClassifierGroup,
                ClassifierName,
                FeatureName,
                FeatureValue,
                CancerType,
                DataValue,
                Rank,
                RankGroup
        );
    }

    public String toString()
    {
        StringBuilder stringBuilder = new StringBuilder();

        List<String> fields = getFieldNames();
        List<Object> values = getValues();

        for(int i = 0; i < fields.size(); i++)
        {
            stringBuilder.append(fields.get(i));
            stringBuilder.append("=");
            stringBuilder.append(values.get(i));
            stringBuilder.append(", ");
        }

        return stringBuilder.toString();
    }
}
