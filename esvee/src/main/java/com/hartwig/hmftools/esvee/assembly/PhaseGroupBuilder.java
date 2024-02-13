package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.AssemblySupport.hasMatchingFragment;
import static com.hartwig.hmftools.esvee.common.SupportType.JUNCTION_MATE;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
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

        buildPhasedAssemblies();

        for(int i = 0; i < mPhaseGroups.size(); ++i)
        {
            mPhaseGroups.get(i).setId(i);
        }
    }

    public List<PhaseGroup> phaseGroups() { return mPhaseGroups; }

    private void findLinkedAssemblies(final JunctionGroup assemblyJunctionGroup, final JunctionAssembly assembly)
    {
        if(assembly.remoteRegions().isEmpty() && assembly.refSideSoftClips().isEmpty())
            return;

        PhaseGroup phaseGroup = assembly.phaseGroup(); // may have been set from an earlier assembly link
        boolean linksWithExisting = phaseGroup != null;
        boolean hasBranchedAssemblies = !assembly.branchedAssemblies().isEmpty();

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
            RefSideSoftClip refSideSoftClip = assembly.refSideSoftClips().get(0);

            if(positionWithin(refSideSoftClip.Position, assemblyJunctionGroup.minPosition(), assemblyJunctionGroup.maxPosition()))
            {
                linkedJunctionGroups.add(assemblyJunctionGroup);
            }
        }

        if(hasBranchedAssemblies)
            linkedJunctionGroups.add(assemblyJunctionGroup);

        for(JunctionGroup junctionGroup : linkedJunctionGroups)
        {
            // the matching to other assemblies for each of this assembly's remote groups is purely for informational purposes
            // but in time may be replaced by actual linking to all expected remote mate reads

            for(JunctionAssembly otherAssembly : junctionGroup.junctionAssemblies())
            {
                if(assembly == otherAssembly)
                    continue;

                if(phaseGroup != null && otherAssembly.phaseGroup() == phaseGroup)
                {
                    // already linked
                    continue;
                }

                boolean branchedAssemblies = hasBranchedAssemblies && assembly.hasBranchedAssembly(otherAssembly);

                if(!branchedAssemblies && !assembliesShareReads(assembly, otherAssembly))
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
                        if(otherAssembly.phaseGroup() != phaseGroup)
                        {
                            // transfer to the other one
                            phaseGroup.assemblies().forEach(x -> otherAssembly.phaseGroup().addAssembly(x));
                            phaseGroup = otherAssembly.phaseGroup();
                            linksWithExisting = true;
                        }
                    }
                    else
                    {
                        phaseGroup.addAssembly(otherAssembly);
                    }
                }
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

    private static boolean assembliesShareReads(final JunctionAssembly first, final JunctionAssembly second)
    {
        // tests matching reads in both the junction reads and any extension reads (ie discordant)
        for(AssemblySupport support : first.support())
        {
            if(support.type() == JUNCTION_MATE)
                continue;

            if(hasMatchingFragment(second.support(), support.read()))
                return true;
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

    private static final int PHASE_GROUP_LOG_COUNT = 1000;

    private void buildPhasedAssemblies()
    {
        SV_LOGGER.info("building phase sets from {} phase groups", mPhaseGroups.size());

        int processed = 0;
        for(PhaseGroup phaseGroup : mPhaseGroups)
        {
            // where there are more than 2 assemblies, start with the ones with the most support and overlapping junction reads
            PhaseSetBuilder phaseSetBuilder = new PhaseSetBuilder(phaseGroup);

            phaseSetBuilder.buildPhaseSets();

            ++processed;

            if((processed % PHASE_GROUP_LOG_COUNT) == 0)
            {
                SV_LOGGER.debug("processed {} phase groups", processed);
            }
        }
    }
}
