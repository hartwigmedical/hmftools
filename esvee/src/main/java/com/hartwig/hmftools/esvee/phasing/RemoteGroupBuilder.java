package com.hartwig.hmftools.esvee.phasing;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.REMOTE_PHASING_MIN_READS;
import static com.hartwig.hmftools.esvee.types.AssemblySupport.hasMatchingFragment;
import static com.hartwig.hmftools.esvee.types.SupportType.JUNCTION_MATE;
import static com.hartwig.hmftools.esvee.phasing.PhaseGroupBuilder.linkToPhaseGroups;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.types.AssemblySupport;
import com.hartwig.hmftools.esvee.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.types.JunctionGroup;
import com.hartwig.hmftools.esvee.types.PhaseGroup;
import com.hartwig.hmftools.esvee.types.RefSideSoftClip;
import com.hartwig.hmftools.esvee.types.RemoteRegion;
import com.hartwig.hmftools.esvee.types.ThreadTask;
import com.hartwig.hmftools.esvee.output.PhaseGroupBuildWriter;
import com.hartwig.hmftools.esvee.read.Read;

class RemoteGroupBuilder extends ThreadTask
{
    private final AssemblyConfig mConfig;
    private final PhaseGroupBuildWriter mWriter;
    private final Queue<JunctionGroup> mJunctionGroups;
    private final Map<String, List<JunctionGroup>> mJunctionGroupMap;

    private final Set<PhaseGroup> mPhaseGroupsSets;
    private final List<PhaseGroup> mRemovedPhaseGroups;
    private final int mJunctionGroupCount;

    public RemoteGroupBuilder(
            final AssemblyConfig config, final Queue<JunctionGroup> junctionGroups,
            final Map<String,List<JunctionGroup>> junctionGroupMap, final PhaseGroupBuildWriter writer)
    {
        super("RemotePhaseGroups");

        mConfig = config;
        mWriter = writer;
        mJunctionGroups = junctionGroups;
        mJunctionGroupMap = junctionGroupMap;

        mJunctionGroupCount = junctionGroups.size();
        mPhaseGroupsSets = Sets.newHashSet();
        mRemovedPhaseGroups = Lists.newArrayList();
    }

    public Set<PhaseGroup> phaseGroups()
    {
        return mPhaseGroupsSets;
    }

    public List<PhaseGroup> removedPhaseGroups()
    {
        return mRemovedPhaseGroups;
    }

    private static final int LOG_COUNT = 10000;

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                int remainingCount = mJunctionGroups.size();
                int processedCount = mJunctionGroupCount - remainingCount;

                mPerfCounter.start();

                ++processedCount;

                JunctionGroup junctionGroup = mJunctionGroups.remove();

                if((processedCount % LOG_COUNT) == 0)
                {
                    SV_LOGGER.debug("processed {} junction groups into phase groups", processedCount, mPhaseGroupsSets.size());
                }

                for(JunctionAssembly assembly : junctionGroup.junctionAssemblies())
                {
                    findRemotePhasedAssemblies(junctionGroup, assembly);
                }

                stopCheckLog(junctionGroup.toString(), mConfig.PerfLogTime);
            }
            catch(NoSuchElementException e)
            {
                SV_LOGGER.trace("all phase tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    private void findRemotePhasedAssemblies(final JunctionGroup junctionGroup, final JunctionAssembly assembly)
    {
        // for the given assembly, looks in all overlapping other junction groups (based on remote regions, ref-side soft-clips, and
        // the local junction group for indels) for other assemblies with shared reads
        if(assembly.remoteRegions().isEmpty() && assembly.refSideSoftClips().isEmpty())
            return;

        Set<JunctionGroup> linkedJunctionGroups = Sets.newHashSet();

        for(RemoteRegion region : assembly.remoteRegions())
        {
            if(isFiltered(region))
                continue;

            // IDEA: first of all check the regions overlap, since for large junction groups there's a high chance they won't
            // and then all reads need to be checked in each
            // finds too many examples where they do share reads, maybe the 1K distance buffer is too small or they are sharing
            // supplementaries whose remote regions were purged

            List<JunctionGroup> overlappingJunctions = findOverlappingJunctionGroups(junctionGroup, region);

            if(overlappingJunctions == null)
                continue;

            linkedJunctionGroups.addAll(overlappingJunctions);
        }

        if(!assembly.refSideSoftClips().isEmpty())
        {
            // if the assembly has candidate facing TI links, then add its own junction group - they will often share reads which are
            // discordant in one and a junction read in the other
            RefSideSoftClip refSideSoftClip = assembly.refSideSoftClips().get(0);

            if(positionWithin(refSideSoftClip.Position, junctionGroup.minPosition(), junctionGroup.maxPosition()))
            {
                linkedJunctionGroups.add(junctionGroup);
            }
        }

        for(JunctionGroup otherJunctionGroup : linkedJunctionGroups)
        {
            // the matching to other assemblies for each of this assembly's remote groups is purely for informational purposes
            // but in time may be replaced by actual linking to all expected remote mate reads

            for(JunctionAssembly otherAssembly : otherJunctionGroup.junctionAssemblies())
            {
                if(assembly == otherAssembly)
                    continue;

                checkAssemblyPhasing(assembly, otherAssembly, mPhaseGroupsSets, mRemovedPhaseGroups, mWriter);
            }
        }
    }

    private synchronized static void checkAssemblyPhasing(
            final JunctionAssembly assembly, final JunctionAssembly otherAssembly,
            final Set<PhaseGroup> phaseGroups, final List<PhaseGroup> removedPhaseGroups, final PhaseGroupBuildWriter writer)
    {
        // this method is synchronised so that only task at a time can attempt set and update an assembly's phase group
        // and the phase group itself

        // for assemblies sharing read, the phase group scenarios are:
        // - assemblies are already in the same phase group
        // - no phase group exists for either so make a new one
        // - the other assembly already has an existing phase group and this one doesn't so just add it
        // - the other assembly already has an existing phase group, so does this one and so transfer this to the other

        PhaseGroup phaseGroup = assembly.phaseGroup();

        if(phaseGroup != null && otherAssembly.phaseGroup() == phaseGroup) // already linked
            return;

        if(!assembliesShareReads(assembly, otherAssembly, REMOTE_PHASING_MIN_READS))
            return;

        linkToPhaseGroups(phaseGroup, assembly, otherAssembly, phaseGroups, removedPhaseGroups, writer, "Remote");
    }

    private static boolean assembliesShareReads(final JunctionAssembly first, final JunctionAssembly second, int minSharedReads)
    {
        // tests matching reads in both the junction reads and any extension reads (ie discordant)
        Set<String> readsMatched = Sets.newHashSet();

        for(AssemblySupport support : first.support())
        {
            if(support.type() == JUNCTION_MATE)
                break;

            Read read = support.read();

            if(hasMatchingFragment(second.support(), read) || hasMatchingFragment(second.candidateSupport(), read))
            {
                readsMatched.add(read.getName());

                if(readsMatched.size() >= minSharedReads)
                    return true;
            }
        }

        // search amongst candidates for a remote juncton read match
        for(AssemblySupport support : first.candidateSupport())
        {
            if(hasMatchingFragment(second.support(), support.read()))
            {
                readsMatched.add(support.read().getName());

                if(readsMatched.size() >= minSharedReads)
                    return true;
            }
        }

        return false;
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
