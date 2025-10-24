package com.hartwig.hmftools.redux.duplicate;

import htsjdk.samtools.SAMRecord;

public interface IReadCache
{
    void processRead(final SAMRecord read);
    int currentReadMinPosition();
    FragmentCoordReads popReads();
    FragmentCoordReads evictAll();
    int minCachedReadStart();
    int cachedReadCount();
    int cachedFragCoordGroups();
}
