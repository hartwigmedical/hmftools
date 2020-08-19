package com.hartwig.hmftools.cup.common;

import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CANCER_SUBTYPE_OTHER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_OVARY;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_PROSTATE;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_UTERUS;

import java.util.Map;

import com.hartwig.hmftools.common.purple.gender.Gender;

public class SampleData
{
    public final String Id;
    public final String CancerType;
    public final String CancerSubtype;
    public final String OriginalCancerType;

    private boolean mIsRefSample;
    private Gender mGenderType;

    public SampleData(final String id, final String cancerType, final String cancerSubtype, final String origType)
    {
        Id = id;
        CancerType = cancerType;
        CancerSubtype = cancerSubtype;
        OriginalCancerType = origType;

        mGenderType = null;
        mIsRefSample = false;
    }

    public Gender gender() { return mGenderType; }
    public void setGender(final Gender gender) { mGenderType = gender; }

    public boolean isRefSample() { return mIsRefSample; }
    public void setRefSample() { mIsRefSample = true; }

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

    public boolean isCandidateCancerType(final String cancerType)
    {
        if(mGenderType == null)
            return true;

        if(cancerType.contains(CANCER_TYPE_UTERUS) || cancerType.contains(CANCER_TYPE_OVARY))
        {
            if(mGenderType != Gender.FEMALE)
                return false;
        }

        if(cancerType.contains(CANCER_TYPE_PROSTATE))
        {
            if(mGenderType == Gender.FEMALE)
                return false;
        }

        return true;
    }
}
