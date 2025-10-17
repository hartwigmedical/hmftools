package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.sequencing.SbxBamUtils.isHomopolymerLowBaseQualAtStart;

import htsjdk.samtools.SAMRecord;

public class SbxBamSampler
{
    private boolean mHpLowQualBaseDirectionSet;

    public SbxBamSampler()
    {
        mHpLowQualBaseDirectionSet = false;
    }

    public void processSampleRecord(final SAMRecord record)
    {
        if(mHpLowQualBaseDirectionSet)
            return;

        Boolean hpLowBaseQualAtStart = isHomopolymerLowBaseQualAtStart(record);

        if(hpLowBaseQualAtStart == null)
            return;

        mHpLowQualBaseDirectionSet = true;

        if(hpLowBaseQualAtStart == !record.getReadNegativeStrandFlag())
            SbxRoutines.SBX_HOMOPOLYMER_5_PRIME_LOW_BASE_QUAL = true;
        else
            SbxRoutines.SBX_HOMOPOLYMER_5_PRIME_LOW_BASE_QUAL = false;
    }
}
