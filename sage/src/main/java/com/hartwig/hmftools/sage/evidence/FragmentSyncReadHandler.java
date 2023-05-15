package com.hartwig.hmftools.sage.evidence;

import htsjdk.samtools.SAMRecord;

public interface FragmentSyncReadHandler
{
    void processReadRecord(final SAMRecord record, boolean checkSync);
}
