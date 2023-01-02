package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

import java.util.List;
import java.util.Map;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import com.google.common.collect.Maps;

public class PartitionCache
{
    // TODO: make private onve GroupCombiner is deleted
    public final String mChrPartition;

    // fragment status from resolved fragments, keyed by readId
    public final Map<String,FragmentStatus> FragmentStatus;

    // supplmentary and candidate duplicate reads, keyed by readId
    public final Map<String,Fragment> IncompleteFragments;

    // positions with candidate duplicate fragments, keyed by initial fragment coordinate position
    public final Map<Integer, List<Fragment>> CandidateDuplicatesMap;

    public PartitionCache(final String chrPartition)
    {
        mChrPartition = chrPartition;
        FragmentStatus = Maps.newHashMap();
        IncompleteFragments = Maps.newHashMap();
        CandidateDuplicatesMap = Maps.newHashMap();
    }

    public int incompleteFragments()
    {
        return IncompleteFragments.size();
    }
    public int resolvedFragments()
    {
        return FragmentStatus.size();
    }

    public void clear()
    {
        FragmentStatus.clear();
        IncompleteFragments.clear();
        CandidateDuplicatesMap.clear();
    }

    public String toString() { return format("%s status(%d) frags(%d) positions(%d)",
            mChrPartition, FragmentStatus.size(), IncompleteFragments.size(), CandidateDuplicatesMap.size()); }
}
