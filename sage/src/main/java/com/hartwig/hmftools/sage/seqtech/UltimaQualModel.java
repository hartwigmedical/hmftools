package com.hartwig.hmftools.sage.seqtech;

import htsjdk.samtools.SAMRecord;

public abstract class UltimaQualModel
{
    private final UltimaModelType mType;

    public UltimaQualModel(final UltimaModelType type)
    {
        mType = type;
    }

    public UltimaModelType type() { return mType; }

    public abstract byte calculateQual(final SAMRecord record, int varReadIndex);

    public boolean canCompute() { return true; }
}
