package com.hartwig.hmftools.isofox.fusion;

import java.util.Map;
import java.util.concurrent.ConcurrentMap;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.Nullable;

public class RacFragmentCache
{
    private ConcurrentMap<String, Map<Integer, JunctionRacFragments>> mGroups;
    private int mTotalFragmentCount;

    public RacFragmentCache()
    {
        mGroups = Maps.newConcurrentMap();
        mTotalFragmentCount = 0;
    }

    public void addRacFragments(final String chromosome, int geneCollectionId, final JunctionRacFragments group)
    {
        if(group.fragmentCount() == 0)
            return;

        Map<Integer,JunctionRacFragments> chrGroups = mGroups.get(chromosome);

        if(chrGroups == null)
        {
            chrGroups = Maps.newHashMap();
            mGroups.put(chromosome, chrGroups);
        }

        chrGroups.put(geneCollectionId, group);
        mTotalFragmentCount += group.fragmentCount();
    }

    @Nullable
    public JunctionRacFragments getRacFragments(final String chromosome, int geneCollectionId)
    {
        Map<Integer,JunctionRacFragments> chrGroups = mGroups.get(chromosome);

        if(chrGroups == null)
            return null;

        return chrGroups.get(geneCollectionId);
    }

    public int totalFragmentCount() { return mTotalFragmentCount; }

    public int assignedFragmentCount()
    {
        return mGroups.values().stream().mapToInt(x -> x.values().stream().mapToInt(y -> y.assignedFragmentCount()).sum()).sum();
    }

    public int totalGroupCount() { return mGroups.values().stream().mapToInt(x -> x.size()).sum(); }
}