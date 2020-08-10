package com.hartwig.hmftools.cup.common;

import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CANCER_SUBTYPE_OTHER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;

import java.util.Map;

public class SampleData
{
    public final String Id;
    public final String CancerType;
    public final String CancerSubtype;
    public final String OriginalCancerType;

    public SampleData(final String id, final String cancerType, final String cancerSubtype, final String origType)
    {
        Id = id;
        CancerType = cancerType;
        CancerSubtype = cancerSubtype;
        OriginalCancerType = origType;
    }

    public static SampleData from(final Map<String,Integer> fieldsIndexMap, final String data)
    {
        final String[] items = data.split(DATA_DELIM, -1);
        final String sampleId = items[fieldsIndexMap.get("SampleId")];
        final String cancerType = items[fieldsIndexMap.get("CancerType")];

        final String cancerSubtype = fieldsIndexMap.containsKey("CancerSubtype") ?
                items[fieldsIndexMap.get("CancerSubtype")] : CANCER_SUBTYPE_OTHER;

        final String originalCancerType = fieldsIndexMap.containsKey("OrigCancerType") ?
                items[fieldsIndexMap.get("OrigCancerType")] : cancerType;

        return new SampleData(sampleId, cancerType, cancerSubtype, originalCancerType);

    }
}
