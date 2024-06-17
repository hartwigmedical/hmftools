package com.hartwig.hmftools.id;

import static java.lang.String.format;

import static com.hartwig.hmftools.id.HmfIdConfig.DATA_DELIM;

import java.util.Comparator;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.amber.AmberAnonymous;

public class SampleData
{
    public final int PatientId;
    public final int SampleIndex;
    public final boolean Deleted;

    public String SampleId; // original non-anonymised ID
    private String mSampleHash;

    private static final String PREFIX = "HMF";

    // file columns
    public static final String PATIENT_ID = "patientId";
    public static final String SAMPLE_INDEX = "sampleId";
    public static final String HMF_SAMPLE_ID = "hmfSampleId";
    public static final String SAMPLE_HASH = "sampleHash";
    public static final String DELETED = "deleted";

    public SampleData(final int patientId, final int sampleIndex, final String sampleId, final String sampleHash, final boolean deleted)
    {
        PatientId = patientId;
        SampleIndex = sampleIndex;
        SampleId = sampleId;
        Deleted = deleted;
        mSampleHash = sampleHash;
    }

    public void setSampleHash(final String hash) { mSampleHash = hash; }
    public String sampleHash() { return mSampleHash; }

    public static SampleData fromAmberAnonymous(final AmberAnonymous amberAnonymous)
    {
        String hmfSampleId = amberAnonymous.HmfSampleId;

        int sampleIndex = Integer.valueOf(hmfSampleId.charAt(hmfSampleId.length() - 1)) - 64;
        String patientIdStr = hmfSampleId.substring(PREFIX.length(), hmfSampleId.length() - 1);
        int patientId = Integer.parseInt(patientIdStr);

        return new SampleData(patientId, sampleIndex, amberAnonymous.SampleId, "", amberAnonymous.Deleted);
    }

    public String toString()
    {
        return format("index(patient=%d sample=%d) sampleId(%s) hmfId(%s) deleted(%s) hash(%s)",
                PatientId, SampleIndex, SampleId, hmfSampleId(), Deleted, mSampleHash);
    }

    public String hmfSampleId()
    {
        // allows for up to 1 million distinct patient IDs
        String patientString = format("%s%06d", PREFIX, PatientId);
        char sampleString = (char)(64 + SampleIndex);
        return patientString + sampleString;
    }

    public boolean matches(final SampleData other)
    {
        return PatientId == other.PatientId && SampleIndex == other.SampleIndex;
    }

    public static class SampleComparator implements Comparator<SampleData>
    {
        public int compare(final SampleData first, final SampleData second)
        {
            if(first.PatientId != second.PatientId)
                return first.PatientId < second.PatientId ? -1 : 1;

            if(first.SampleIndex != second.SampleIndex)
                return first.SampleIndex < second.SampleIndex ? -1 : 1;

            return 0;
        }
    }

    public AmberAnonymous toAmberAnonymous()
    {
        return new AmberAnonymous(hmfSampleId(), SampleId, Deleted);
    }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);
        sj.add(PATIENT_ID);
        sj.add(SAMPLE_INDEX);
        sj.add(HMF_SAMPLE_ID);
        sj.add(SAMPLE_HASH);
        sj.add(DELETED);
        return sj.toString();
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);
        sj.add(String.valueOf(PatientId));
        sj.add(String.valueOf(SampleIndex));
        sj.add(hmfSampleId());
        sj.add(String.valueOf(Deleted));
        sj.add(mSampleHash);
        return sj.toString();
    }


}
