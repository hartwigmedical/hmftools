package com.hartwig.hmftools.svtools.sequence;

public class InsertSeqData
{
    public final String SampleId;
    public final int SvId;
    public final String InsertSeq;

    public InsertSeqData(final String sampleId, final int svId, final String insertSeq)
    {
        SampleId = sampleId;
        InsertSeq = insertSeq;
        SvId = svId;
    }

}
