package com.hartwig.hmftools.sage.evidence;

import htsjdk.samtools.SAMRecord;

public class FragmentSyncOutcome
{
    public final SAMRecord CombinedRecord;
    public final FragmentSyncType SyncType;

    public static final String ORIG_READ_COORDS = "OC";

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
