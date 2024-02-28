package com.hartwig.hmftools.sage.quality;

import htsjdk.samtools.SAMRecord;

public abstract class UltimaQualModel
{
    private final UltimaVariantModel mType;

    public UltimaQualModel(final UltimaVariantModel type)
    {
        mType = type;
    }

    public UltimaVariantModel type() { return mType; }

    public abstract byte calculateQual(final SAMRecord record, int varReadIndex);
}
