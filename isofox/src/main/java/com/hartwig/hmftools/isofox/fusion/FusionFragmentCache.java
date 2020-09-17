package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.isofox.fusion.ReadGroup.mergeChimericReadMaps;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.isofox.IsofoxConfig;

public class FusionFragmentCache
{
    private final Map<String,Map<Integer,List<FusionFragment>>> mChrRealignCandidates;
    private final Map<String,ReadGroup> mIncompleteReadGroups;

    // gene filters
    private final IsofoxConfig mConfig;

    public FusionFragmentCache(final IsofoxConfig config)
    {
        mChrRealignCandidates = Maps.newHashMap();
        mIncompleteReadGroups = Maps.newHashMap();

        mConfig = config;
    }

    public synchronized List<ReadGroup> addIncompleteGroups(final Map<String,ReadGroup> incompleteGroups)
    {
        final List<ReadGroup> completeGroups = Lists.newArrayList();
        mergeChimericReadMaps(mIncompleteReadGroups, completeGroups, incompleteGroups);
        return completeGroups;
    }

}
