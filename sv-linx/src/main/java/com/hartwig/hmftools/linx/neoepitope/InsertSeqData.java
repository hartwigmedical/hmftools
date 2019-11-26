package com.hartwig.hmftools.linx.neoepitope;

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
