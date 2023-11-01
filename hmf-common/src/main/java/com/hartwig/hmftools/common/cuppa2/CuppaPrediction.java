package com.hartwig.hmftools.common.cuppa2;

import java.util.ArrayList;
import java.util.List;

public class CuppaPrediction
{
    public static final String FLD_SAMPLE_ID = "sample_id";
    public static final String FLD_DATA_TYPE = "data_type";
    public static final String FLD_CLF_GROUP = "clf_group"; // Potentially unused
    public static final String FLD_CLF_NAME = "clf_name";
    public static final String FLD_FEAT_NAME = "feat_name";
    public static final String FLD_FEAT_VALUE = "feat_value";
    public static final String FLD_CANCER_TYPE = "cancer_type";
    public static final String FLD_DATA_VALUE = "data_value";

    public final String SampleId;
    public final Categories.DataType DataType;
    public final Categories.ClfGroup ClfGroup;
    public final Categories.ClfName ClfName;
    public final String FeatName;
    public final double FeatValue;
    public final String CancerType;
    public final double DataValue;

    CuppaPrediction(
            final String sampleId,
            final Categories.DataType dataType,
            final Categories.ClfGroup clfGroup,
            final Categories.ClfName clfName,
            final String featName,
            final double featValue,
            final String cancerType,
            final double dataValue
    )
    {
        SampleId = sampleId;
        DataType = dataType;
        ClfGroup = clfGroup;
        ClfName = clfName;
        FeatName = featName;
        FeatValue = featValue;
        CancerType = cancerType;
        DataValue = dataValue;
    }

    public List<Object> getFieldNames()
    {

        List<Object> fields = new ArrayList<>();

        fields.add(FLD_SAMPLE_ID);
        fields.add(FLD_DATA_TYPE);
        fields.add(FLD_CLF_GROUP);
        fields.add(FLD_CLF_NAME);
        fields.add(FLD_FEAT_NAME);
        fields.add(FLD_FEAT_VALUE);
        fields.add(FLD_CANCER_TYPE);
        fields.add(FLD_DATA_VALUE);

        return fields;
    }

    public List<Object> getValues()
    {

        List<Object> values = new ArrayList<>();

        values.add(SampleId);
        values.add(DataType);
        values.add(ClfGroup);
        values.add(ClfName);
        values.add(FeatName);
        values.add(FeatValue);
        values.add(CancerType);
        values.add(DataValue);

        return values;
    }

    public String toString()
    {
        StringBuilder stringBuilder = new StringBuilder();

        List<Object> fields = getFieldNames();
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
