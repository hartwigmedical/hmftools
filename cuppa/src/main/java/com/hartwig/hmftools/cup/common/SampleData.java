package com.hartwig.hmftools.cup.common;

import static com.hartwig.hmftools.cup.CuppaConfig.CANCER_SUBTYPE_OTHER;
import static com.hartwig.hmftools.cup.CuppaConfig.FLD_CANCER_SUBTYPE;
import static com.hartwig.hmftools.cup.CuppaConfig.FLD_CANCER_TYPE;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.CuppaConfig.FLD_RNA_READ_LENGTH;
import static com.hartwig.hmftools.cup.CuppaConfig.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_OTHER;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_UNKNOWN;

import java.util.Map;

import com.hartwig.hmftools.common.purple.Gender;

public class SampleData
{
    public final String Id;
    public final String CancerSubtype;
    public final String PrimaryLocation;
    public final String PrimarySubLocation;
    public final String PrimaryType;
    public final String PrimarySubtype;

    private String mCancerType; // can be overridden if matches a reference sample
    private boolean mIsRefSample;
    private Gender mGenderType;
    private int mRnaReadLength;

    public static final int RNA_READ_LENGTH_NONE = 0;

    public SampleData(final String id, final String cancerType, final String cancerSubtype)
    {
        this(id, cancerType, cancerSubtype, "", "", "", "");
    }

    public SampleData(
            final String id, final String cancerType, final String cancerSubtype,
            final String location, final String locationSubtype, final String primaryType, final String primarySubtype)
    {
        Id = id;
        mCancerType = cancerType;
        CancerSubtype = cancerSubtype;
        PrimaryLocation = location;
        PrimarySubLocation = locationSubtype;
        PrimaryType = primaryType;
        PrimarySubtype = primarySubtype;

        mGenderType = null;
        mIsRefSample = false;
        mRnaReadLength = RNA_READ_LENGTH_NONE;
    }

    public String cancerType()
    {
        return !CancerSubtype.isEmpty() && !CancerSubtype.equals(CANCER_TYPE_UNKNOWN) ? CancerSubtype : mCancerType ;
    }

    public String cancerMainType() { return mCancerType; }

    public void setCancerType(final String cancerType) { mCancerType = cancerType; }

    public Gender gender() { return mGenderType; }
    public void setGender(final Gender gender) { mGenderType = gender; }

    public void setRnaReadLength(int readLength) { mRnaReadLength = readLength; }
    public int rnaReadLength() { return mRnaReadLength; }
    public boolean hasRna() { return mRnaReadLength != RNA_READ_LENGTH_NONE; }

    public boolean isRefSample() { return mIsRefSample; }
    public void setRefSample() { mIsRefSample = true; }

    public boolean isKnownCancerType() { return isKnownCancerType(mCancerType); }

    public static boolean isKnownCancerType(final String cancerType)
    {
        return !cancerType.equals(CANCER_TYPE_OTHER) && !cancerType.equals(CANCER_TYPE_UNKNOWN);
    }

    public static SampleData from(final Map<String,Integer> fieldsIndexMap, final String data)
    {
        final String[] items = data.split(DATA_DELIM, -1);
        final String sampleId = items[fieldsIndexMap.get(FLD_SAMPLE_ID)];

        String cancerType = extractOptionalField(fieldsIndexMap, items, FLD_CANCER_TYPE, CANCER_TYPE_UNKNOWN);
        String cancerSubtype = extractOptionalField(fieldsIndexMap, items, FLD_CANCER_SUBTYPE, CANCER_TYPE_UNKNOWN);
        String primaryLocation = extractOptionalField(fieldsIndexMap, items, "PrimaryTumorLocation", CANCER_SUBTYPE_OTHER);
        String primarySubLocation = extractOptionalField(fieldsIndexMap, items, "PrimaryTumorSubLocation", "");
        String primaryType = extractOptionalField(fieldsIndexMap, items, "PrimaryTumorType", "");
        String primarySubtype = extractOptionalField(fieldsIndexMap, items, "PrimaryTumorSubType", "");

        SampleData sample = new SampleData(sampleId, cancerType, cancerSubtype, primaryLocation, primarySubLocation, primaryType, primarySubtype);

        if(fieldsIndexMap.containsKey(FLD_RNA_READ_LENGTH))
        {
            sample.setRnaReadLength(Integer.parseInt(items[fieldsIndexMap.get(FLD_RNA_READ_LENGTH)]));
        }

        return sample;
    }

    private static String extractOptionalField(
            final Map<String,Integer> fieldsIndexMap, final String[] items, final String fieldName, final String defaultValue)
    {
        return fieldsIndexMap.containsKey(fieldName) ?items[fieldsIndexMap.get(fieldName)] : defaultValue;
    }

    public boolean isCandidateCancerType(final String cancerType) { return CupConstants.isCandidateCancerType(mGenderType, cancerType); }
}
