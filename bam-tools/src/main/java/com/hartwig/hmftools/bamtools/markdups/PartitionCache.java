package com.hartwig.hmftools.bamtools.markdups;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public class PartitionCache
{
    public final String ChrPartition;

    // fragment status from resolved fragments, keyed by (remote) chromosome-partition then readId
    public final Map<String, FragmentStatus> FragmentStatus;

    // incomplete fragments (unclear or supplmentaries), keyed by chromosome-partition then readId
    public final Map<String, Fragment> IncompleteFragments;

    // positions with candidate duplicate fragments, keyed by chromosome-partition then initial fragment coordinate position
    public final Map<Integer, List<Fragment>> IncompleteFragmentPositions;

    public PartitionCache(final String chrPartition)
    {
        ChrPartition = chrPartition;
        FragmentStatus = Maps.newHashMap();
        IncompleteFragments = Maps.newHashMap();
        IncompleteFragmentPositions = Maps.newHashMap();
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
        IncompleteFragmentPositions.clear();
    }
}
