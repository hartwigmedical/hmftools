package com.hartwig.hmftools.cup.common;

import java.util.Map;

import com.google.common.collect.Maps;

public class SampleResult
{
    public final String SampleId;
    public final CategoryType Category;
    public final ResultType ResultType;
    public final String DataType;
    public final Object Value;

    public final Map<String,Double> CancerTypeValues;

    public SampleResult(
            final String sampleId, final CategoryType category, final ResultType resultType, final String dataType, final Object value,
            final Map<String,Double> cancerTypeValues)
    {
        SampleId = sampleId;
        Category = category;
        DataType = dataType;
        ResultType = resultType;
        Value = value;
        CancerTypeValues = cancerTypeValues;
    }

    public static boolean checkIsValidCancerType(final SampleData sample, final String refCancerType, final Map<String,Double> cancerDataMap)
    {
        if(sample.isCandidateCancerType(refCancerType))
            return true;

        cancerDataMap.put(refCancerType, 0.0);
        return false;
    }

    public String toString()
    {
        return String.format("sample(%s) cat(%s) resultType(%s) type(%s) value(%s)",
                SampleId, Category, ResultType, DataType, Value.toString());
    }
}
