package com.hartwig.hmftools.common.amber;

public class AmberAnonymous
{
    public final String SampleId;
    public final String HmfSampleId;
    public final boolean Deleted;

    public AmberAnonymous(final String sampleId, final String hmfSampleId, final boolean deleted)
    {
        SampleId = sampleId;
        HmfSampleId = hmfSampleId;
        Deleted = deleted;
    }
}
