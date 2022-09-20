package com.hartwig.hmftools.id;

import static java.lang.String.format;

import static com.hartwig.hmftools.id.HmfIdConfig.DATA_DELIM;

import java.util.Comparator;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.amber.AmberAnonymous;
import com.hartwig.hmftools.common.amber.ImmutableAmberAnonymous;

public class HmfSample
{
    public final int PatientId;
    public final int SampleId;
    public final String SampleHash;
    public final boolean Deleted;

    /*
    public static final String PATIENT_ID = "PatientId";
    public static final String SAMPLE_ID = "SampleId";
    public static final String HMF_SAMPLE_ID = "HmfSampleId";
    public static final String SAMPLE_HASH = "SampleHash";
    public static final String DELETED = "Deleted";
    */

    private static final String PREFIX = "HMF";

    public HmfSample(final int patientId, final int sampleId, final String sampleHash, final boolean deleted)
    {
        PatientId = patientId;
        SampleId = sampleId;
        SampleHash = sampleHash;
        Deleted = deleted;
    }

    public static HmfSample fromAmberAnonymous(final AmberAnonymous amberAnonymous)
    {
        String hmfSampleId = amberAnonymous.hmfSampleId();

        // int sampleId = hmfSampleId[hmfSampleId.length() - 1].code - 64;
        int sampleId = Integer.valueOf(hmfSampleId.charAt(hmfSampleId.length() - 1)) - 64;
        String patientIdStr = hmfSampleId.substring(PREFIX.length(), hmfSampleId.length() - 1);
        int patientId = Integer.parseInt(patientIdStr);

        return new HmfSample(patientId, sampleId, amberAnonymous.sampleId(), amberAnonymous.deleted());
    }

    public String toString() { return format("patient(%d) sample(%d) hash(%s) deleted(%s)", PatientId, SampleId, SampleHash, Deleted); }

    public String hmfSample()
    {
        String patientString = format("%s%06d", PREFIX, PatientId);
        char sampleString = (char)(64 + SampleId);
        return patientString + sampleString;
    }

    public static class SampleComparator implements Comparator<HmfSample>
    {
        public int compare(final HmfSample first, final HmfSample second)
        {
            if(first.PatientId != second.PatientId)
                return first.PatientId < second.PatientId ? -1 : 1;

            if(first.SampleId != second.SampleId)
                return first.SampleId < second.SampleId ? -1 : 1;

            return 0;
        }
    }

    public AmberAnonymous toAmberAnonymous()
    {
        return ImmutableAmberAnonymous.builder()
                .hmfSampleId(hmfSample())
                .sampleId(SampleHash)
                .deleted(Deleted)
                .build();
    }

    public static final String PATIENT_ID = "PatientId";
    public static final String SAMPLE_ID = "SampleId";
    public static final String HMF_SAMPLE_ID = "HmfSampleId";
    public static final String SAMPLE_HASH = "SampleHash";
    public static final String DELETED = "Deleted";

    public static String header()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);
        sj.add(PATIENT_ID);
        sj.add(SAMPLE_ID);
        sj.add(HMF_SAMPLE_ID);
        sj.add(SAMPLE_HASH);
        sj.add(DELETED);
        return sj.toString();
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);
        sj.add(PATIENT_ID);
        sj.add(SAMPLE_ID);
        sj.add(HMF_SAMPLE_ID);
        sj.add(SAMPLE_HASH);
        sj.add(DELETED);
        return sj.toString();
    }


}
