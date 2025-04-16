package com.hartwig.hmftools.redux;

import com.hartwig.hmftools.redux.common.FragmentCoordReads;
import com.hartwig.hmftools.redux.common.FragmentCoords;

import htsjdk.samtools.SAMRecord;

public interface IReadCache
{
    void processRead(final SAMRecord read);
    boolean isEmpty();
    int currentReadMinPosition();
    FragmentCoordReads popReads();
    FragmentCoordReads evictAll();
    int minCachedReadStart();
    int cachedReadCount();
    int cachedFragCoordGroups();
}
