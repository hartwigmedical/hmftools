package com.hartwig.hmftools.cup.common;

import static com.hartwig.hmftools.cup.CuppaConfig.CANCER_SUBTYPE_OTHER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
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

    private boolean mIsRefSample;
    private Gender mGenderType;
    private boolean mCheckSubType;

    public SampleData(final String id, final String cancerType, final String cancerSubtype)
    {
        Id = id;
        CancerType = cancerType;
        CancerSubtype = cancerSubtype;

        mGenderType = null;
        mIsRefSample = false;
        mCheckSubType = false;
    }

    public Gender gender() { return mGenderType; }
    public void setGender(final Gender gender) { mGenderType = gender; }

    public boolean isRefSample() { return mIsRefSample; }
    public void setRefSample() { mIsRefSample = true; }

    public void setCheckSubType() { mCheckSubType = true; }
    public boolean checkSubType() { return mCheckSubType; }

    public static SampleData from(final Map<String,Integer> fieldsIndexMap, final String data)
    {
        final String[] items = data.split(DATA_DELIM, -1);
        final String sampleId = items[fieldsIndexMap.get("SampleId")];

        String cancerType = items[fieldsIndexMap.get("CancerType")];

        String cancerSubtype = fieldsIndexMap.containsKey("CancerSubtype") ?
                items[fieldsIndexMap.get("CancerSubtype")] : CANCER_SUBTYPE_OTHER;

        return new SampleData(sampleId, cancerType, cancerSubtype);
    }

    public static String CANCER_SUBTYPE_DELIM = ";";
    public static int CANCER_ITEM = 0;
    public static int CANCER_SUBTYPE_ITEM = 1;

    public static String parseCancerType(final String cancerSubTypeCombo)
    {
        final String[] items = cancerSubTypeCombo.split(CANCER_SUBTYPE_DELIM);
        return items.length == 2 ? items[CANCER_ITEM] : "";
    }

    public static String parseCancersSubType(final String cancerSubTypeCombo)
    {
        final String[] items = cancerSubTypeCombo.split(CANCER_SUBTYPE_DELIM);
        return items.length == 2 ? items[CANCER_SUBTYPE_ITEM] : "";
    }

    public boolean isCandidateCancerType(final String cancerType)
    {
        if(mCheckSubType)
        {
            String refCancerType = cancerType.contains(CANCER_SUBTYPE_DELIM) ? parseCancerType(cancerType) : cancerType;
            return parseCancerType(CancerType).equals(refCancerType);
        }

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
