package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.assembliesShareReads;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.DiscordantGroup;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.JunctionGroup;
import com.hartwig.hmftools.esvee.common.PhaseGroup;
import com.hartwig.hmftools.esvee.common.RefSideSoftClip;
import com.hartwig.hmftools.esvee.common.RemoteRegion;

public class PhaseGroupBuilder
{
    private final SvConfig mConfig;
    private final Map<String, List<JunctionGroup>> mJunctionGroupMap;
    private final List<PhaseGroup> mPhaseGroups;

    public PhaseGroupBuilder(final SvConfig config, final Map<String, List<JunctionGroup>> junctionGroupMap)
    {
        mConfig = config;
        mJunctionGroupMap = junctionGroupMap;
        mPhaseGroups = Lists.newArrayList();

        // setting an index for each junction group allows easy access to an assembly's own group during look-ups
        for(List<JunctionGroup> junctionGroups : mJunctionGroupMap.values())
        {
            for(int index = 0; index < junctionGroups.size(); ++index)
            {
                junctionGroups.get(index).setIndex(index);
            }
        }
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
                        SV_LOGGER.debug("processed {} assemblies into {} phase groups",
                                processed, mPhaseGroups.size());
                    }
                }
            }
        }

        for(int i = 0; i < mPhaseGroups.size(); ++i)
        {
            mPhaseGroups.get(i).setId(i);
        }
    }

    public List<PhaseGroup> phaseGroups() { return mPhaseGroups; }

    private void findLinkedAssemblies(final JunctionGroup assemblyJunctionGroup, final JunctionAssembly assembly)
    {
        // for the given assembly, looks in all overlapping other junction groups (based on remote regions, ref-side soft-clips, and
        // the local junction group for indels) for other assemblies with shared reads
        if(assembly.remoteRegions().isEmpty() && assembly.refSideSoftClips().isEmpty() && !assembly.indel())
            return;

        PhaseGroup phaseGroup = assembly.phaseGroup(); // may have been set from an earlier assembly link
        boolean linksWithExisting = phaseGroup != null;

        Set<JunctionGroup> linkedJunctionGroups = Sets.newHashSet();

        for(RemoteRegion region : assembly.remoteRegions())
        {
            if(isFiltered(region))
                continue;

            List<JunctionGroup> overlappingJunctions = findOverlappingJunctionGroups(assemblyJunctionGroup, region);

            if(overlappingJunctions == null)
                continue;

            linkedJunctionGroups.addAll(overlappingJunctions);
        }

        if(!assembly.refSideSoftClips().isEmpty())
        {
            // if the assembly has candidate facing TI links, then add its own junction group - they will often share reads which are
            // discordant in one and a junction read in the other
            RefSideSoftClip refSideSoftClip = assembly.refSideSoftClips().get(0);

            if(positionWithin(refSideSoftClip.Position, assemblyJunctionGroup.minPosition(), assemblyJunctionGroup.maxPosition()))
            {
                linkedJunctionGroups.add(assemblyJunctionGroup);
            }
        }

        if(assembly.indel())
            linkedJunctionGroups.add(assemblyJunctionGroup);

        for(JunctionGroup junctionGroup : linkedJunctionGroups)
        {
            // the matching to other assemblies for each of this assembly's remote groups is purely for informational purposes
            // but in time may be replaced by actual linking to all expected remote mate reads

            for(JunctionAssembly otherAssembly : junctionGroup.junctionAssemblies())
            {
                if(assembly == otherAssembly)
                    continue;

                // for assemblies sharing read, the phase group scenarios are:
                // - assemblies are already in the same phase group
                // - no phase group exists for either so make a new one
                // - the other assembly already has an existing phase group and this one doesn't so just add it
                // - the other assembly already has an existing phase group, so does this one and so transfer this to the other
                if(phaseGroup != null && otherAssembly.phaseGroup() == phaseGroup)
                {
                    // already linked
                    continue;
                }

                // FIXME: could first of all check the regions overlap, since for large junction groups there's a high chance they won't
                // and then all reads need to be checked in each

                if(!assembliesShareReads(assembly, otherAssembly))
                    continue;

                if(phaseGroup == null)
                {
                    if(otherAssembly.phaseGroup() != null)
                    {
                        phaseGroup = otherAssembly.phaseGroup();
                        phaseGroup.addAssembly(assembly);
                        linksWithExisting = true;
                    }
                    else
                    {
                        phaseGroup = new PhaseGroup(assembly, otherAssembly);
                    }
                }
                else
                {
                    if(otherAssembly.phaseGroup() != null)
                    {
                        // transfer to the other one, and do this for assemblies in this group
                        otherAssembly.phaseGroup().transferAssemblies(phaseGroup);
                        mPhaseGroups.remove(phaseGroup);
                        phaseGroup = otherAssembly.phaseGroup();
                        linksWithExisting = true;
                    }
                    else
                    {
                        phaseGroup.addAssembly(otherAssembly);
                    }
                }
            }
        }

        for(JunctionGroup junctionGroup : linkedJunctionGroups)
        {
            for(DiscordantGroup discordantGroup : junctionGroup.discordantGroups())
            {
                // as above, check that remote regions overlap before checking reads
                if(assembly.remoteRegions().stream().noneMatch(x -> x.overlaps(discordantGroup.remoteRegion())))
                    continue;

                if(assembly.support().stream().noneMatch(x -> discordantGroup.hasRead(x.read())))
                    continue;

                if(phaseGroup == null)
                    phaseGroup = new PhaseGroup(assembly, null);

                phaseGroup.addDiscordantGroup(discordantGroup);
            }
        }

        if(phaseGroup != null && !linksWithExisting)
        {
            mPhaseGroups.add(phaseGroup);
        }
    }

    private List<JunctionGroup> findOverlappingJunctionGroups(final JunctionGroup assemblyJunctionGroup, final RemoteRegion region)
    {
        List<JunctionGroup> junctionGroups = mJunctionGroupMap.get(region.Chromosome);

        if(junctionGroups == null)
            return Collections.emptyList();

        int startIndex;

        if(assemblyJunctionGroup.chromosome().equals(region.Chromosome) && assemblyJunctionGroup.overlapsRemoteRegion(region))
        {
            startIndex = assemblyJunctionGroup.index();
        }
        else
        {
            // make use of the binary search for junction groups to find the closest index
            startIndex = JunctionGroup.binarySearch(region.start(), junctionGroups);
        }

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
