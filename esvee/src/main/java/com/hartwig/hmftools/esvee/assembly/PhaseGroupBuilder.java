package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.region.BaseRegion.binarySearch;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.JunctionGroup;
import com.hartwig.hmftools.esvee.common.PrimaryPhaseGroup;
import com.hartwig.hmftools.esvee.common.RemoteRegion;

public class PhaseGroupBuilder
{
    private final SvConfig mConfig;
    private final Map<String, List<JunctionGroup>> mJunctionGroupMap;
    private final List<PrimaryPhaseGroup> mPrimaryPhaseGroups;
    private int mMissingRemoteGroups;

    public PhaseGroupBuilder(final SvConfig config, final Map<String, List<JunctionGroup>> junctionGroupMap)
    {
        mConfig = config;
        mJunctionGroupMap = junctionGroupMap;
        mPrimaryPhaseGroups = Lists.newArrayList();
        mMissingRemoteGroups = 0;
    }

    private static final int LOG_COUNT = 10000;

    public void buildGroups()
    {
        int processed = 0;
        for(List<JunctionGroup> junctionGroups : mJunctionGroupMap.values())
        {
            for(JunctionGroup junctionGroup : junctionGroups)
            {
                for(JunctionAssembly assembly : junctionGroup.junctionAssemblies())
                {
                    findLinkedAssemblies(junctionGroup, assembly);

                    ++processed;

                    if((processed % LOG_COUNT) == 0)
                    {
                        SV_LOGGER.debug("primary phasing processed {} assemblies, groups({})",
                                processed, mPrimaryPhaseGroups.size());
                    }
                }
            }
        }

        for(int i = 0; i < mPrimaryPhaseGroups.size(); ++i)
        {
            mPrimaryPhaseGroups.get(i).setId(i);
        }
    }

    public List<PrimaryPhaseGroup> primaryPhaseGroups() { return mPrimaryPhaseGroups; }
    public int missingRemoteGroups() { return mMissingRemoteGroups; }

    private void findLinkedAssemblies(final JunctionGroup assemblyJunctionGroup, final JunctionAssembly assembly)
    {
        if(assembly.remoteRegions().isEmpty())
            return;

        PrimaryPhaseGroup primaryPhaseGroup = assembly.primaryPhaseGroup(); // may have been set from an earlier assembly link
        boolean linksWithExisting = primaryPhaseGroup != null;

        Set<JunctionGroup> processedGroups = Sets.newHashSet();

        for(RemoteRegion region : assembly.remoteRegions())
        {
            if(isFiltered(region))
                continue;

            List<JunctionGroup> overlappingJunctions = findOverlappingJunctionGroups(assemblyJunctionGroup, region);

            if(overlappingJunctions == null)
            {
                ++mMissingRemoteGroups;
                continue;
            }

            boolean matched = false;

            for(JunctionGroup junctionGroup : overlappingJunctions)
            {
                if(processedGroups.contains(junctionGroup))
                    continue;

                processedGroups.add(junctionGroup);

                for(JunctionAssembly otherAssembly : junctionGroup.junctionAssemblies())
                {
                    if(assembly == otherAssembly)
                        continue;

                    if(!assembliesShareReads(assembly, otherAssembly))
                        continue;

                    matched = true;

                    if(primaryPhaseGroup == null)
                    {
                        if(otherAssembly.primaryPhaseGroup() != null)
                        {
                            primaryPhaseGroup = otherAssembly.primaryPhaseGroup();
                            primaryPhaseGroup.addAssembly(assembly);
                            linksWithExisting = true;
                        }
                        else
                        {
                            primaryPhaseGroup = new PrimaryPhaseGroup(assembly, otherAssembly);
                        }
                    }
                    else
                    {
                        if(otherAssembly.primaryPhaseGroup() != null)
                        {
                            if(otherAssembly.primaryPhaseGroup() != primaryPhaseGroup)
                            {
                                // transfer to the other one
                                primaryPhaseGroup.assemblies().forEach(x -> otherAssembly.primaryPhaseGroup().addAssembly(x));
                                primaryPhaseGroup = otherAssembly.primaryPhaseGroup();
                                linksWithExisting = true;
                            }
                        }
                        else
                        {
                            primaryPhaseGroup.addAssembly(otherAssembly);
                        }
                    }

                }
            }

            if(!matched)
            {
                ++mMissingRemoteGroups;
            }
        }

        if(primaryPhaseGroup != null && !linksWithExisting)
        {
            mPrimaryPhaseGroups.add(primaryPhaseGroup);
        }
    }

    private List<JunctionGroup> findOverlappingJunctionGroups(final JunctionGroup assemblyJunctionGroup, final RemoteRegion region)
    {
        if(assemblyJunctionGroup.containsRemoteRegion(region))
            return List.of(assemblyJunctionGroup);

        List<JunctionGroup> junctionGroups = mJunctionGroupMap.get(region.Chromosome);

        if(junctionGroups == null)
            return Collections.emptyList();

        // make use of the binary search for junction groups to find the closest index
        int startIndex = JunctionGroup.binarySearch(region.start(), junctionGroups);

        List<JunctionGroup> overlapGroups = Lists.newArrayList();

        for(int d = 0; d <= 1; ++d)
        {
            boolean searchUp = (d == 0);

            int index = searchUp ? startIndex + 1 : startIndex;

            while(index >= 0 && index < junctionGroups.size())
            {
                if(!junctionGroups.get(index).overlapsRemoteRegion(region))
                    break;

                overlapGroups.add(junctionGroups.get(index));

                if(searchUp)
                    ++index;
                else
                    --index;
            }
        }

        return overlapGroups;
    }


    private static boolean assembliesShareReads(final JunctionAssembly first, final JunctionAssembly second)
    {
        for(AssemblySupport support : first.support())
        {
            if(second.hasMatchingFragmentSupport(support.read()))
                return true;
        }

        if(first.refBaseAssembly() != null)
        {
            for(AssemblySupport support : first.refBaseAssembly().support())
            {
                if(second.hasMatchingFragmentSupport(support.read()))
                    return true;
            }
        }

        return false;
    }

    private boolean isFiltered(final RemoteRegion region)
    {
        // ignore any remote region filtered out by config
        if(!mConfig.SpecificChrRegions.hasFilters())
            return false;

        if(mConfig.SpecificChrRegions.excludeChromosome(region.Chromosome))
            return true;

        if(mConfig.SpecificChrRegions.Regions.isEmpty())
            return false;

        return mConfig.SpecificChrRegions.Regions.stream().noneMatch(x -> x.overlaps(region));
    }
}
