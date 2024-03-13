package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.BAM_READ_JUNCTION_BUFFER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.isLocalAssemblyCandidate;
import static com.hartwig.hmftools.esvee.common.AssemblySupport.hasMatchingFragment;
import static com.hartwig.hmftools.esvee.common.SupportType.JUNCTION_MATE;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
import com.hartwig.hmftools.esvee.common.DiscordantGroup;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.JunctionGroup;
import com.hartwig.hmftools.esvee.common.PhaseGroup;
import com.hartwig.hmftools.esvee.common.RefSideSoftClip;
import com.hartwig.hmftools.esvee.common.RemoteRegion;
import com.hartwig.hmftools.esvee.common.ThreadTask;

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

    public List<PhaseGroup> phaseGroups() { return mPhaseGroups; }

    public void buildGroups()
    {
        if(mConfig.PhaseGroupMultiThread)
        {
            List<Thread> threadTasks = new ArrayList<>();

            List<PhaseGroupBuilderTask> groupBuilderTasks = Lists.newArrayList();

            Queue<JunctionGroup> junctionGroupQueue = new ConcurrentLinkedQueue<>();

            for(List<JunctionGroup> junctionGroups : mJunctionGroupMap.values())
            {
                junctionGroupQueue.addAll(junctionGroups);
            }

            int junctionGroupCount = junctionGroupQueue.size();

            int taskCount = min(mConfig.Threads, junctionGroupCount);

            for(int i = 0; i < taskCount; ++i)
            {
                PhaseGroupBuilderTask groupBuilderTask = new PhaseGroupBuilderTask(junctionGroupQueue);
                groupBuilderTasks.add(groupBuilderTask);
                threadTasks.add(groupBuilderTask);
            }

            if(taskCount > 1)
            {
                SV_LOGGER.debug("phase-grouping {} junction groups across {} threads", junctionGroupCount, taskCount);
            }

            if(!runThreadTasks(threadTasks))
                System.exit(1);

            groupBuilderTasks.forEach(x -> mPhaseGroups.addAll(x.phaseGroups()));

            // clean-up phase groups which were transferred into another group
            groupBuilderTasks.forEach(x -> x.removedPhaseGroups().forEach(y -> mPhaseGroups.remove(y)));
        }
        else
        {
            int processed = 0;

            for(List<JunctionGroup> junctionGroups : mJunctionGroupMap.values())
            {
                for(JunctionGroup junctionGroup : junctionGroups)
                {
                    for(JunctionAssembly assembly : junctionGroup.junctionAssemblies())
                    {
                        findLinkedAssembliesNonThreaded(junctionGroup, assembly);

                        ++processed;

                        if((processed % PhaseGroupBuilderTask.LOG_COUNT) == 0)
                        {
                            SV_LOGGER.debug("processed {} assemblies into {} phase groups",
                                    processed, mPhaseGroups.size());
                        }
                    }
                }
            }
        }

        for(int i = 0; i < mPhaseGroups.size(); ++i)
        {
            mPhaseGroups.get(i).setId(i);
        }

        // run validation
        if(mConfig.PerfDebug)
        {
            // check if an assembly is in 2 phase groups
            List<JunctionAssembly> assemblies = Lists.newArrayList();

            for(List<JunctionGroup> junctionGroups : mJunctionGroupMap.values())
            {
                junctionGroups.forEach(x -> assemblies.addAll(x.junctionAssemblies()));
            }

            for(JunctionAssembly assembly : assemblies)
            {
                List<PhaseGroup> phaseGroups = mPhaseGroups.stream().filter(x -> x.assemblies().contains(assembly)).collect(Collectors.toList());

                if(phaseGroups.size() > 1)
                {
                    SV_LOGGER.error("asm({}) in {} phase groups", assembly, phaseGroups.size());
                }
            }
        }
    }

    private class PhaseGroupBuilderTask extends ThreadTask
    {
        private final Queue<JunctionGroup> mJunctionGroups;

        private final Set<PhaseGroup> mPhaseGroupsSets;
        private final List<PhaseGroup> mRemovedPhaseGroups;
        private final int mJunctionGroupCount;

        public PhaseGroupBuilderTask(final Queue<JunctionGroup> junctionGroups)
        {
            super("PhaseGroups");

            mJunctionGroups = junctionGroups;
            mJunctionGroupCount = junctionGroups.size();
            mPhaseGroupsSets = Sets.newHashSet();
            mRemovedPhaseGroups = Lists.newArrayList();
        }

        public Set<PhaseGroup> phaseGroups() { return mPhaseGroupsSets; }
        public List<PhaseGroup> removedPhaseGroups() { return mRemovedPhaseGroups; }

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
                        findLinkedAssemblies(junctionGroup, assembly);
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

        private void findLinkedAssemblies(final JunctionGroup junctionGroup, final JunctionAssembly assembly)
        {
            // for the given assembly, looks in all overlapping other junction groups (based on remote regions, ref-side soft-clips, and
            // the local junction group for indels) for other assemblies with shared reads
            boolean hasLocalAssemblyCandidates = junctionGroup.junctionAssemblies().stream().anyMatch(x -> isLocalAssemblyCandidate(assembly, x));

            if(!hasLocalAssemblyCandidates && assembly.remoteRegions().isEmpty() && assembly.refSideSoftClips().isEmpty() && !assembly.indel())
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

            if(assembly.indel() || hasLocalAssemblyCandidates)
                linkedJunctionGroups.add(junctionGroup);

            for(JunctionGroup otherJunctionGroup : linkedJunctionGroups)
            {
                // the matching to other assemblies for each of this assembly's remote groups is purely for informational purposes
                // but in time may be replaced by actual linking to all expected remote mate reads

                for(JunctionAssembly otherAssembly : otherJunctionGroup.junctionAssemblies())
                {
                    if(assembly == otherAssembly)
                        continue;

                    boolean sameJunctionGroup = otherJunctionGroup == junctionGroup;

                    checkAssemblyPhasing(
                            assembly, otherAssembly, sameJunctionGroup, mPhaseGroupsSets, mRemovedPhaseGroups, mConfig.LogPhaseGroupLinks);
                }
            }

            if(!mConfig.SkipDiscordant)
            {
                for(JunctionGroup otherJunctionGroup : linkedJunctionGroups)
                {
                    // the linked junction group may have discordant read groups which are required by this assembly
                    for(DiscordantGroup discordantGroup : otherJunctionGroup.discordantGroups())
                    {
                        // ignore any discordant group which actually shares the same reads
                        if(discordantGroup.reads().stream().anyMatch(x -> assembly.hasReadSupport(x)))
                            continue;

                        // as above, check that remote regions overlap before checking reads
                        boolean matched = false;

                        for(RemoteRegion remoteRegion : assembly.remoteRegions())
                        {
                            if(remoteRegion.isSuppOnlyRegion())
                                continue;

                            if(!remoteRegion.overlaps(
                                    discordantGroup.chromosome(), discordantGroup.minAlignedPosition(), discordantGroup.maxAlignedPosition(),
                                    discordantGroup.orientation()))
                            {
                                continue;
                            }

                            if(remoteRegion.readIds().stream().anyMatch(x -> discordantGroup.hasFragment(x)))
                            {
                                matched = true;
                                break;
                            }
                        }

                        if(!matched)
                            continue;

                        // FIXME: not thread-safe
                        /*
                        if(phaseGroup == null)
                            phaseGroup = new PhaseGroup(assembly, null);

                        phaseGroup.addDiscordantGroup(discordantGroup);
                        */
                    }
                }
            }
        }
    }

    private synchronized static void checkAssemblyPhasing(
            final JunctionAssembly assembly, final JunctionAssembly otherAssembly, final boolean sameJunctionGroup,
            final Set<PhaseGroup> phaseGroups, final List<PhaseGroup> removedPhaseGroups, final boolean recordLinks)
    {
        // for assemblies sharing read, the phase group scenarios are:
        // - assemblies are already in the same phase group
        // - no phase group exists for either so make a new one
        // - the other assembly already has an existing phase group and this one doesn't so just add it
        // - the other assembly already has an existing phase group, so does this one and so transfer this to the other

        PhaseGroup phaseGroup = assembly.phaseGroup();

        if(phaseGroup != null && otherAssembly.phaseGroup() == phaseGroup) // already linked
            return;

        boolean isLocalCandidateLink = sameJunctionGroup && isLocalAssemblyCandidate(assembly, otherAssembly);

        if(!isLocalCandidateLink && !assembliesShareReads(assembly, otherAssembly))
            return;

        if(phaseGroup == null)
        {
            if(otherAssembly.phaseGroup() != null)
            {
                otherAssembly.phaseGroup().addAssembly(assembly, otherAssembly);
            }
            else
            {
                phaseGroups.add(new PhaseGroup(assembly, otherAssembly, recordLinks));
            }
        }
        else
        {
            if(otherAssembly.phaseGroup() != null)
            {
                // transfer to the other one, and do this for assemblies in this group
                otherAssembly.phaseGroup().transferAssemblies(phaseGroup);
                removedPhaseGroups.add(phaseGroup);
            }
            else
            {
                phaseGroup.addAssembly(otherAssembly, recordLinks ? assembly : null);
            }
        }
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

        for(AssemblySupport support : first.candidateSupport())
        {
            if(hasMatchingFragment(second.support(), support.read()))
                return true;
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

    private static boolean assembliesOverlap(final JunctionAssembly first, final JunctionAssembly second)
    {
        for(RemoteRegion region : first.remoteRegions())
        {
            int regionStart = region.start() - BAM_READ_JUNCTION_BUFFER;
            int regionEnd = region.start() - BAM_READ_JUNCTION_BUFFER;

            if(positionWithin(second.junction().Position, regionStart, regionEnd))
                return true;
        }

        for(RemoteRegion region : second.remoteRegions())
        {
            int regionStart = region.start() - BAM_READ_JUNCTION_BUFFER;
            int regionEnd = region.start() - BAM_READ_JUNCTION_BUFFER;

            if(positionWithin(first.junction().Position, regionStart, regionEnd))
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

    private void findLinkedAssembliesNonThreaded(final JunctionGroup assemblyGroup, final JunctionAssembly assembly)
    {
        // for the given assembly, looks in all overlapping other junction groups (based on remote regions, ref-side soft-clips, and
        // the local junction group for indels) for other assemblies with shared reads
        boolean hasLocalAssemblyCandidates = assemblyGroup.junctionAssemblies().stream().anyMatch(x -> isLocalAssemblyCandidate(assembly, x));

        if(!hasLocalAssemblyCandidates && assembly.remoteRegions().isEmpty() && assembly.refSideSoftClips().isEmpty() && !assembly.indel())
            return;

        PhaseGroup phaseGroup = assembly.phaseGroup(); // may have been set from an earlier assembly link
        boolean linksWithExisting = phaseGroup != null;

        Set<JunctionGroup> linkedJunctionGroups = Sets.newHashSet();

        for(RemoteRegion region : assembly.remoteRegions())
        {
            if(isFiltered(region))
                continue;

            List<JunctionGroup> overlappingJunctions = findOverlappingJunctionGroups(assemblyGroup, region);

            if(overlappingJunctions == null)
                continue;

            linkedJunctionGroups.addAll(overlappingJunctions);
        }

        if(!assembly.refSideSoftClips().isEmpty())
        {
            // if the assembly has candidate facing TI links, then add its own junction group - they will often share reads which are
            // discordant in one and a junction read in the other
            RefSideSoftClip refSideSoftClip = assembly.refSideSoftClips().get(0);

            if(positionWithin(refSideSoftClip.Position, assemblyGroup.minPosition(), assemblyGroup.maxPosition()))
            {
                linkedJunctionGroups.add(assemblyGroup);
            }
        }

        if(assembly.indel() || hasLocalAssemblyCandidates)
            linkedJunctionGroups.add(assemblyGroup);

        for(JunctionGroup junctionGroup : linkedJunctionGroups)
        {
            // the matching to other assemblies for each of this assembly's remote groups is purely for informational purposes
            // but in time may be replaced by actual linking to all expected remote mate reads

            for(JunctionAssembly otherAssembly : junctionGroup.junctionAssemblies())
            {
                if(assembly == otherAssembly)
                    continue;

                boolean sameJunctionGroup = otherAssembly.phaseGroup() == phaseGroup;

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

                boolean isLocalCandidateLink = (assemblyGroup == junctionGroup && isLocalAssemblyCandidate(assembly, otherAssembly));

                if(!isLocalCandidateLink && !assembliesShareReads(assembly, otherAssembly))
                    continue;

                // IDEA: first of all check the regions overlap, since for large junction groups there's a high chance they won't
                // and then all reads need to be checked in each
                // finds too many examples where they do share reads, maybe the 1K distance buffer is too small or they are sharing
                // supplementaries whose remote regions were purged

                if(phaseGroup == null)
                {
                    if(otherAssembly.phaseGroup() != null)
                    {
                        phaseGroup = otherAssembly.phaseGroup();
                        phaseGroup.addAssembly(assembly, otherAssembly);
                        linksWithExisting = true;
                    }
                    else
                    {
                        phaseGroup = new PhaseGroup(assembly, otherAssembly, mConfig.LogPhaseGroupLinks);
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
                        phaseGroup.addAssembly(otherAssembly, assembly);
                    }
                }
            }
        }

        for(JunctionGroup junctionGroup : linkedJunctionGroups)
        {
            // the linked junction group may have discordant read groups which are required by this assembly
            for(DiscordantGroup discordantGroup : junctionGroup.discordantGroups())
            {
                // ignore any discordant group which actually shares the same reads
                if(discordantGroup.reads().stream().anyMatch(x -> assembly.hasReadSupport(x)))
                    continue;

                // as above, check that remote regions overlap before checking reads
                boolean matched = false;

                for(RemoteRegion remoteRegion : assembly.remoteRegions())
                {
                    if(remoteRegion.isSuppOnlyRegion())
                        continue;

                    if(!remoteRegion.overlaps(
                            discordantGroup.chromosome(), discordantGroup.minAlignedPosition(), discordantGroup.maxAlignedPosition(),
                            discordantGroup.orientation()))
                    {
                        continue;
                    }

                    if(remoteRegion.readIds().stream().anyMatch(x -> discordantGroup.hasFragment(x)))
                    {
                        matched = true;
                        break;
                    }
                }

                if(!matched)
                    continue;

                if(phaseGroup == null)
                    phaseGroup = new PhaseGroup(assembly, null, false);

                phaseGroup.addDiscordantGroup(discordantGroup);
            }
        }

        if(phaseGroup != null && !linksWithExisting)
        {
            mPhaseGroups.add(phaseGroup);
        }
    }
}
