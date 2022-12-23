package com.hartwig.hmftools.sage.evidence;

import htsjdk.samtools.SAMRecord;

public class SyncFragmentOutcome
{
    public final SAMRecord CombinedRecord;
    public final SyncFragmentType SyncType;

    public SyncFragmentOutcome(final SyncFragmentType syncType)
    {
        this(null, syncType);
    }

    public SyncFragmentOutcome(final SAMRecord combinedRecord, final SyncFragmentType syncType)
    {
        CombinedRecord = combinedRecord;
        SyncType = syncType;
    }
}
