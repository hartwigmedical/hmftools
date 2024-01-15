package com.hartwig.hmftools.sage.sync;

import htsjdk.samtools.SAMRecord;

public interface FragmentSyncReadHandler
{
    void processReadRecord(final SAMRecord record, boolean checkSync, final FragmentData fragmentData);

    void processReadRecord(final SAMRecord record, boolean checkSync);
}
