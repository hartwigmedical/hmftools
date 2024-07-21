package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.REMOTE_PHASING_MIN_READS;
import static com.hartwig.hmftools.esvee.assembly.phase.PhaseGroupBuilder.linkToPhaseGroups;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.JunctionGroup;
import com.hartwig.hmftools.esvee.assembly.types.PhaseGroup;
import com.hartwig.hmftools.esvee.assembly.types.RefSideSoftClip;
import com.hartwig.hmftools.esvee.assembly.types.RemoteRegion;
import com.hartwig.hmftools.esvee.assembly.types.ThreadTask;
import com.hartwig.hmftools.esvee.assembly.output.PhaseGroupBuildWriter;

public class RemoteGroupBuilder extends ThreadTask
{
    private final AssemblyConfig mConfig;
    private final PhaseGroupBuildWriter mWriter;
    private final Queue<JunctionGroup> mJunctionGroups;
    private final Map<String, List<JunctionGroup>> mJunctionGroupMap;

    private final Set<PhaseGroup> mPhaseGroupsSets;
    private final List<PhaseGroup> mRemovedPhaseGroups;
    private final int mJunctionGroupCount;

    private final RemoteBuildStats mBuildStats;

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
        mBuildStats = new RemoteBuildStats();
    }

    public Set<PhaseGroup> phaseGroups()
    {
        return mPhaseGroupsSets;
    }

    public List<PhaseGroup> removedPhaseGroups()
    {
        return mRemovedPhaseGroups;
    }

    public void logStats()
    {
        if(mConfig.PerfDebug || mConfig.PerfLogTime > 0)
        {
            // now appears inconsequential
            // SV_LOGGER.debug("remote phase group building stats: {}", mBuildStats);
        }
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

            List<JunctionGroup> overlappingJunctions = findOverlappingJunctionGroups(junctionGroup, region);

            if(overlappingJunctions == null)
                continue;

            linkedJunctionGroups.addAll(overlappingJunctions);
        }

        // CHECK: really necessary if already checked by local phase group (same junction) building?
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
            for(JunctionAssembly otherAssembly : otherJunctionGroup.junctionAssemblies())
            {
                if(assembly == otherAssembly)
                    continue;

                if(doCheckAssemblyPhasing(assembly, otherAssembly, false, null, null, null))
                {
                    ++mBuildStats.AssemblyChecks;
                    ++mBuildStats.AssemblyPreLinked;
                    return;
                }

                if(!canPhaseAssemblies(assembly, otherAssembly))
                    continue;

                doCheckAssemblyPhasing(assembly, otherAssembly, true, mPhaseGroupsSets, mRemovedPhaseGroups, mWriter);
            }
        }
    }

    private boolean canPhaseAssemblies(final JunctionAssembly assembly, final JunctionAssembly otherAssembly)
    {
        ++mBuildStats.AssemblyChecks;

        RemoteRegion overlappingRegion = assembly.remoteRegions().stream()
                .filter(x -> x.overlaps(
                        otherAssembly.junction().Chromosome, otherAssembly.minAlignedPosition(), otherAssembly.maxAlignedPosition()))
                .findFirst().orElse(null);

        if(overlappingRegion == null)
            return false;

        if(!assembliesShareReads(overlappingRegion, otherAssembly, REMOTE_PHASING_MIN_READS, mConfig.ApplyRemotePhasingReadCheckThreshold))
        {
            ++mBuildStats.AssemblyNonMatches;
            return false;
        }

        ++mBuildStats.AssemblyMatches;

        return true;
    }

    private synchronized static boolean doCheckAssemblyPhasing(
            final JunctionAssembly assembly, final JunctionAssembly otherAssembly, boolean formLink,
            final Set<PhaseGroup> phaseGroups, final List<PhaseGroup> removedPhaseGroups, final PhaseGroupBuildWriter writer)
    {
        if(formLink)
        {
            // check again if assemblies have been linked while this thread was checking if they needed to be
            if(!assembliesAlreadyPhased(assembly, otherAssembly))
            {
                linkToPhaseGroups(assembly.phaseGroup(), assembly, otherAssembly, phaseGroups, removedPhaseGroups, writer, "Remote");
            }

            return true;
        }
        else
        {
            return assembliesAlreadyPhased(assembly, otherAssembly);
        }
    }

    private static boolean assembliesAlreadyPhased(final JunctionAssembly first, final JunctionAssembly second)
    {
        // this and the linking method may need to be managed with a lock
        return first.phaseGroup() != null && first.phaseGroup() == second.phaseGroup();
    }

    private boolean assembliesShareReads(
            final RemoteRegion firstRegion, final JunctionAssembly second, int minSharedReads, boolean applyMatchThreshold)
    {
        // tests matching reads in both the junction reads and any extension reads (ie discordant)

        int firstReadCount = firstRegion.readIds().size();
        int secondReadCount = second.supportCount() + second.candidateSupport().size();

        if(firstReadCount < minSharedReads || secondReadCount < minSharedReads)
            return false;

        int maxMatchChecks = applyMatchThreshold ?
                calculateMaxFragmentCheckThreshold(firstReadCount, secondReadCount) : NO_FRAG_CHECK_THRESHOLD;

        int currentChecks = 0;
        int matchedCount = 0;

        for(String readId : firstRegion.readIds())
        {
            if(hasMatchingSupportRead(second.support(), readId, mBuildStats))
            {
                ++matchedCount;
                ++mBuildStats.ReadMatches;

                if(matchedCount >= minSharedReads)
                    return true;

                continue;
            }

            currentChecks += second.supportCount();

            if(maxMatchChecks != NO_FRAG_CHECK_THRESHOLD && currentChecks > maxMatchChecks)
            {
                ++mBuildStats.ReadCheckLimited;
                return false;
            }

            if(hasMatchingRead(second.candidateSupport(), readId, mBuildStats))
            {
                ++matchedCount;
                ++mBuildStats.ReadMatches;

                if(matchedCount >= minSharedReads)
                    return true;

                continue;
            }

            currentChecks += second.candidateSupport().size();

            if(maxMatchChecks != NO_FRAG_CHECK_THRESHOLD && currentChecks > maxMatchChecks)
            {
                ++mBuildStats.ReadCheckLimited;
                return false;
            }
        }

        return false;
    }

    private static boolean hasMatchingSupportRead(final List<SupportRead> support, final String readId, final RemoteBuildStats stats)
    {
        for(SupportRead supportRead : support)
        {
            ++stats.ReadChecks;

            if(supportRead.id().equals(readId))
            {
                ++stats.ReadMatches;
                return true;
            }
        }

        stats.ReadNonMatches += support.size();
        return false;
    }

    private static boolean hasMatchingRead(final List<Read> support, final String readId, final RemoteBuildStats stats)
    {
        for(Read supportRead : support)
        {
            ++stats.ReadChecks;

            if(supportRead.id().equals(readId))
            {
                ++stats.ReadMatches;
                return true;
            }
        }

        stats.ReadNonMatches += support.size();
        return false;
    }

    private static final int NO_FRAG_CHECK_THRESHOLD = -1;

    private static int calculateMaxFragmentCheckThreshold(int firstReadCount, int secondReadCount)
    {
        long maxChecks = firstReadCount * secondReadCount;

        if(maxChecks < 500)
            return NO_FRAG_CHECK_THRESHOLD;

        if(maxChecks <= 1000)
            return (int)round(0.6 * maxChecks);

        if(maxChecks <= 2000)
            return (int)round(0.5 * maxChecks);

        if(maxChecks <= 5000)
            return (int)round(0.4 * maxChecks);

        if(maxChecks <= 10000)
            return (int)round(0.3 * maxChecks);

        if(maxChecks <= 250000)
            return (int)round(0.1 * maxChecks);

        return 50000;
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

    private static class RemoteBuildStats
    {
        public long AssemblyChecks;
        public long AssemblyPreLinked;
        public long AssemblyMatches;
        public long AssemblyNonMatches;
        public long ReadChecks;
        public long ReadMatches;
        public long ReadNonMatches;
        public long ReadCheckLimited;

        public RemoteBuildStats()
        {
            AssemblyChecks = 0;
            AssemblyPreLinked = 0;
            AssemblyMatches = 0;
            AssemblyNonMatches = 0;
            ReadChecks = 0;
            ReadMatches = 0;
            ReadNonMatches = 0;
            ReadCheckLimited = 0;
        }

        public String toString()
        {
            return format("assembly(check=%d preLinked=%d matched=%d nonMatch=%d) reads(%d matched=%d nonMatch=%d limited=%d)",
                    AssemblyChecks, AssemblyPreLinked, AssemblyMatches, AssemblyNonMatches,
                    ReadChecks, ReadMatches, ReadNonMatches, ReadCheckLimited);
        }
    }
}
