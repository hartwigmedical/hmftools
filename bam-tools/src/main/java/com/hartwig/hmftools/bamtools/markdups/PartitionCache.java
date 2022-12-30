package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public class PartitionCache
{
    public final String ChrPartition;

    // fragment status from resolved fragments, keyed by readId
    public final Map<String,FragmentStatus> FragmentStatus;

    // supplmentary reads, keyed by readId
    public final Map<String,Fragment> Supplementaries;

    public PartitionCache(final String chrPartition)
    {
        ChrPartition = chrPartition;
        FragmentStatus = Maps.newHashMap();
        Supplementaries = Maps.newHashMap();
    }

    public int incompleteFragments()
    {
        return Supplementaries.size();
    }

    public int resolvedFragments()
    {
        return FragmentStatus.size();
    }

    public void clear()
    {
        FragmentStatus.clear();
        Supplementaries.clear();
    }

    public String toString() { return format("%s status(%d) frags(%d)",
            ChrPartition, FragmentStatus.size(), Supplementaries.size()); }
}
