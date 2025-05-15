package com.hartwig.hmftools.redux.common;

import com.google.common.collect.Multiset;

// TODO: remove
public interface IAuxReadCache
{
    int minCachedReadStart();
    int cachedReadCount();
    int cachedFragCoordGroups();
    Multiset<String> readNames();
    void pushSingleRead(final ReadInfo readInfo);
    void pushDuplicateGroup(final DuplicateGroup duplicateGroup);
    FragmentCoordReads popReads(final int lastReadCacheBoundary);
    FragmentCoordReads evictAll();
}
