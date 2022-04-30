package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGN_CANDIDATE;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentMap;

import com.google.common.collect.Maps;

import org.apache.commons.compress.utils.Lists;
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

    public ConcurrentMap<String,Map<Integer,JunctionRacFragments>> getAllGroups() { return mGroups; }

    public List<FusionFragment> getUnassignedFragments()
    {
        List<FusionFragment> fragments = Lists.newArrayList();

        for(Map<Integer,JunctionRacFragments> gcEntry : mGroups.values())
        {
            for(JunctionRacFragments juncRacFragments : gcEntry.values())
            {
                juncRacFragments.getRacFragments().stream().filter(x -> x.type() == REALIGN_CANDIDATE).forEach(x -> fragments.add(x));
            }
        }

        return fragments;
    }


}