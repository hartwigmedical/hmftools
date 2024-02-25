package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.BAM_READ_JUNCTION_BUFFER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyDeduper.dedupProximateAssemblies;

import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.common.JunctionGroup;
import com.hartwig.hmftools.esvee.common.ThreadTask;
import com.hartwig.hmftools.esvee.output.ResultsWriter;
import com.hartwig.hmftools.esvee.read.BamReader;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.read.ReadAdjustments;
import com.hartwig.hmftools.esvee.read.ReadFilters;
import com.hartwig.hmftools.esvee.read.ReadStats;

import htsjdk.samtools.SAMRecord;

public class JunctionGroupAssembler extends ThreadTask
{
    private final SvConfig mConfig;

    private final Queue<JunctionGroup> mJunctionGroups;
    private final int mJunctionCount;

    private JunctionGroup mCurrentJunctionGroup;
    private final BamReader mBamReader;

    private final Map<String,ReadGroup> mReadGroupMap;
    private final ReadStats mReadStats;

    public JunctionGroupAssembler(
            final SvConfig config, final BamReader bamReader, final Queue<JunctionGroup> junctionGroups)
    {
        super("PrimaryAssembly");
        mConfig = config;
        mBamReader = bamReader;
        mJunctionGroups = junctionGroups;
        mJunctionCount = junctionGroups.size();

        mReadGroupMap = Maps.newHashMap();
        mCurrentJunctionGroup = null;
        mReadStats = new ReadStats();
    }

    public static List<JunctionGroupAssembler> createThreadTasks(
            final List<JunctionGroup> junctionGroups, final List<BamReader> bamReaders, final SvConfig config,
            final int taskCount, final List<Thread> threadTasks)
    {
        List<JunctionGroupAssembler> primaryAssemblyTasks = Lists.newArrayList();

        Queue<JunctionGroup> junctionGroupQueue = new ConcurrentLinkedQueue<>();
        junctionGroupQueue.addAll(junctionGroups);

        int junctionGroupCount = junctionGroups.size();

        for(int i = 0; i < taskCount; ++i)
        {
            BamReader bamReader = bamReaders.get(i);

            JunctionGroupAssembler junctionGroupAssembler = new JunctionGroupAssembler(config, bamReader, junctionGroupQueue);
            primaryAssemblyTasks.add(junctionGroupAssembler);
            threadTasks.add(junctionGroupAssembler);
        }

        SV_LOGGER.debug("splitting {} junction groups across {} threads", junctionGroupCount, taskCount);

        return primaryAssemblyTasks;
    }

    private static final int TASK_LOG_COUNT = 1000;

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                int remainingCount = mJunctionGroups.size();
                int processedCount = mJunctionCount - remainingCount;

                JunctionGroup junctionGroup = mJunctionGroups.remove();

                mPerfCounter.start();
                processJunctionGroup(junctionGroup);

                stopCheckLog(junctionGroup.toString(), mConfig.PerfLogTime);

                if(processedCount > 0 && (processedCount % TASK_LOG_COUNT) == 0)
                {
                    SV_LOGGER.debug("processed {} junction groups, remaining({})", processedCount, remainingCount);
                }
            }
            catch(NoSuchElementException e)
            {
                SV_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    public ReadStats readStats() { return mReadStats; }

    private void processJunctionGroup(final JunctionGroup junctionGroup)
    {
        mCurrentJunctionGroup = junctionGroup;

        int sliceStart = junctionGroup.readRangeStart();
        int sliceEnd = junctionGroup.readRangeEnd();

        SV_LOGGER.trace("junctionGroup({}:{}-{} count={}) slice",
                junctionGroup.chromosome(), sliceStart, sliceEnd, junctionGroup.count());

        mBamReader.sliceBam(junctionGroup.chromosome(), sliceStart, sliceEnd, this::processRecord);

        SV_LOGGER.trace("junctionGroup({}:{}-{} count={}) slice complete, readCount({}) readGroups({})",
                junctionGroup.chromosome(), sliceStart, sliceEnd, junctionGroup.count(), junctionGroup.candidateReadCount(),
                mReadGroupMap.size());

        mReadGroupMap.values().forEach(x -> x.formReadLinks());
        mReadGroupMap.clear();

        List<JunctionAssembly> junctionGroupAssemblies = Lists.newArrayList();

        // now pass applicable reads to each junction assembler - any read overlapping the junction
        // due to SvPrep filtering, most reads crossing the junction will have met soft-clip criteria
        for(int i = 0; i < junctionGroup.junctions().size(); ++i)
        {
            Junction junction = junctionGroup.junctions().get(i);

            // FIXME: decide how to model these - for now their reads are captured but not used explicitly used
            if(junction.DiscordantOnly)
                continue;

            JunctionAssembler junctionAssembler = new JunctionAssembler(mConfig, junction);

            // FIXME: doesn't seem to be making a big difference, but this is in efficient for long-range junction groups
            // since both the junctions and reads are ordered. Could consider re-ordering by unclipped start and comparing to junction position

            int junctionBoundaryStart = junction.isForward() ? junction.Position - BAM_READ_JUNCTION_BUFFER : junction.Position;
            int junctionBoundaryEnd = junction.isForward() ? junction.Position : junction.Position + BAM_READ_JUNCTION_BUFFER;

            List<Read> junctionCandidateReads = junctionGroup.candidateReads().stream()
                    .filter(x -> positionsOverlap(junctionBoundaryStart, junctionBoundaryEnd, x.unclippedStart(), x.unclippedEnd()))
                    .collect(Collectors.toList());

            if(junctionCandidateReads.isEmpty())
                continue;

            List<JunctionAssembly> candidateAssemblies = junctionAssembler.processJunction(junctionCandidateReads);

            // dedup assemblies with close junction positions, same orientation
            dedupProximateAssemblies(junctionGroupAssemblies, candidateAssemblies);

            // junctionGroupAssemblies.addAll(candidateAssemblies);

            // extend assemblies with non-junction and discordant reads
            for(JunctionAssembly assembly : candidateAssemblies)
            {
                AssemblyExtender assemblyExtender = new AssemblyExtender(assembly);
                assemblyExtender.findAssemblyExtensions(junctionAssembler.nonJunctionReads());

                junctionGroupAssemblies.addAll(assemblyExtender.assemblies());
            }
        }

        junctionGroup.addJunctionAssemblies(junctionGroupAssemblies);

        mCurrentJunctionGroup = null;
    }

    private void processRecord(final SAMRecord record)
    {
        mConfig.logReadId(record, "JunctionGroupAssembler:processRecord");

        Read read = new Read(record);

        ++mReadStats.TotalReads;

        // CHECK: do in SvPrep if worthwhile
        if(!mConfig.NoReadFilters)
        {
            if(!ReadFilters.isAboveBaseQualAvgThreshold(record.getBaseQualities()))
            {
                ++mReadStats.FilteredBaseQual;
                return;
            }
            else if(!ReadFilters.isAboveBaseQualAvgThreshold(record.getBaseQualities()) || !ReadFilters.isAboveMapQualThreshold(read))
            {
                ++mReadStats.FilteredMapQual;
                return;
            }
        }

        if(mBamReader.currentIsReferenceSample())
            read.markReference();

        // CHECK: could track for stats
        if(ReadAdjustments.trimPolyGSequences(read))
            ++mReadStats.PolyGTrimmed;

        if(ReadAdjustments.convertEdgeIndelsToSoftClip(read))
            ++mReadStats.IndelSoftClipConverted;

        mCurrentJunctionGroup.addCandidateRead(read);

        ReadGroup readGroup = mReadGroupMap.get(read.getName());

        if(readGroup != null)
        {
            readGroup.addRead(read);
            return;
        }

        // link first and second in pair if within the same group
        boolean hasLocalMate = read.isMateMapped() && read.mateChromosome().equals(read.chromosome())
                && positionsOverlap(
                        read.mateAlignmentStart(), read.mateAlignmentEnd(),
                        mCurrentJunctionGroup.readRangeStart(), mCurrentJunctionGroup.readRangeEnd());

        // link first and second in pair if within the same group
        boolean hasLocalSupplementary = read.hasSupplementary() && read.supplementaryData().Chromosome.equals(read.chromosome())
                && positionWithin(
                        read.supplementaryData().Position, mCurrentJunctionGroup.readRangeStart(), mCurrentJunctionGroup.readRangeEnd());

        if(!hasLocalMate && !hasLocalSupplementary)
            return;

        int expectedCount = 1 + (hasLocalMate ? 1 : 0) + (hasLocalSupplementary ? 1 : 0); // approximate only for array size
        readGroup = new ReadGroup(read, expectedCount);
        mReadGroupMap.put(read.getName(), readGroup);
    }

    private class ReadGroup
    {
        private final List<Read> mReads;
        private final int mExpectedCount;

        public ReadGroup(final Read read, int expectedCount)
        {
            mExpectedCount = expectedCount;
            mReads = Lists.newArrayListWithCapacity(expectedCount);
            mReads.add(read);
        }

        public void addRead(final Read read)
        {
            mReads.add(read);
        }

        public void formReadLinks()
        {
            if(mReads.size() == 1)
                return;

            for(int i = 0; i < mReads.size() - 1; ++i)
            {
                Read first = mReads.get(i);

                for(int j = i + 1; j < mReads.size(); ++j)
                {
                    Read second = mReads.get(j);

                    first.makeReadLinks(second);
                }
            }
        }

        public String toString() { return format("reads(%d exp=%d) id(%s)", mReads.size(), mExpectedCount, mReads.get(0).getName()); }
    }
}
