package com.hartwig.hmftools.sage.sync;

import htsjdk.samtools.SAMRecord;

public class FragmentSyncOutcome
{
    public final SAMRecord CombinedRecord;
    public final FragmentSyncType SyncType;

    public FragmentSyncOutcome(final FragmentSyncType syncType)
    {
        this(null, syncType);
    }

    public FragmentSyncOutcome(final SAMRecord combinedRecord, final FragmentSyncType syncType)
    {
        CombinedRecord = combinedRecord;
        SyncType = syncType;
    }
}
