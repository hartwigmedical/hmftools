package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.BAM_READ_JUNCTION_BUFFER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyDeduper.dedupProximateAssemblies;

import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.alignment.DecoyChecker;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionGroup;
import com.hartwig.hmftools.esvee.assembly.types.ThreadTask;
import com.hartwig.hmftools.esvee.assembly.output.ResultsWriter;
import com.hartwig.hmftools.esvee.assembly.read.BamReader;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments;
import com.hartwig.hmftools.esvee.assembly.read.ReadStats;

import htsjdk.samtools.SAMRecord;

public class JunctionGroupAssembler extends ThreadTask
{
    private final AssemblyConfig mConfig;

    private final Queue<JunctionGroup> mJunctionGroups;
    private final int mJunctionCount;

    private JunctionGroup mCurrentJunctionGroup;
    private final BamReader mBamReader;
    private final DecoyChecker mDecoyChecker;

    private final Map<String,ReadGroup> mReadGroupMap;
    private final Map<String,SAMRecord> mSupplementaryRepeats; // temporary to track an issue in SvPrep
    private final ReadStats mReadStats;

    public JunctionGroupAssembler(
            final AssemblyConfig config, final BamReader bamReader, final Queue<JunctionGroup> junctionGroups, final ResultsWriter resultsWriter)
    {
        super("PrimaryAssembly");
        mConfig = config;
        mBamReader = bamReader;
        mJunctionGroups = junctionGroups;
        mJunctionCount = junctionGroups.size();

        mDecoyChecker = new DecoyChecker(mConfig.DecoyGenome, resultsWriter.decoyMatchWriter());

        mReadGroupMap = Maps.newHashMap();
        mSupplementaryRepeats = Maps.newHashMap();
        mCurrentJunctionGroup = null;
        mReadStats = new ReadStats();
    }

    public static List<JunctionGroupAssembler> createThreadTasks(
            final List<JunctionGroup> junctionGroups, final List<BamReader> bamReaders, final AssemblyConfig config,
            final ResultsWriter resultsWriter, final int taskCount, final List<Thread> threadTasks)
    {
        List<JunctionGroupAssembler> primaryAssemblyTasks = Lists.newArrayList();

        Queue<JunctionGroup> junctionGroupQueue = new ConcurrentLinkedQueue<>();
        junctionGroupQueue.addAll(junctionGroups);

        int junctionGroupCount = junctionGroups.size();

        for(int i = 0; i < taskCount; ++i)
        {
            BamReader bamReader = bamReaders.get(i);

            JunctionGroupAssembler junctionGroupAssembler = new JunctionGroupAssembler(config, bamReader, junctionGroupQueue, resultsWriter);
            primaryAssemblyTasks.add(junctionGroupAssembler);
            threadTasks.add(junctionGroupAssembler);
        }

        if(taskCount > 1)
        {
            SV_LOGGER.debug("primary assembly splits {} junction groups across {} threads", junctionGroupCount, taskCount);
        }

        return primaryAssemblyTasks;
    }

    private static final int TASK_LOG_COUNT = 10000;

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
        mSupplementaryRepeats.clear();

        List<JunctionAssembly> junctionGroupAssemblies = Lists.newArrayList();

        RefBaseExtender refBaseExtender = new RefBaseExtender();
        DiscordantReads discordantReads = new DiscordantReads();

        // now pass applicable reads to each junction assembler - any read overlapping the junction
        // due to SvPrep filtering, most reads crossing the junction will have met soft-clip criteria
        for(int i = 0; i < junctionGroup.junctions().size(); ++i)
        {
            Junction junction = junctionGroup.junctions().get(i);

            JunctionAssembler junctionAssembler = new JunctionAssembler(junction);

            // doesn't seem to be making a big difference, but this is inefficient for long-range junction groups
            // since both the junctions and reads are ordered. Could consider re-ordering by unclipped start and comparing to junction position

            int junctionBoundaryStart = junction.isForward() ? junction.Position - BAM_READ_JUNCTION_BUFFER : junction.Position;
            int junctionBoundaryEnd = junction.isForward() ? junction.Position : junction.Position + BAM_READ_JUNCTION_BUFFER;

            List<Read> candidateReads = junctionGroup.candidateReads().stream()
                    .filter(x -> positionsOverlap(junctionBoundaryStart, junctionBoundaryEnd, x.unclippedStart(), x.unclippedEnd()))
                    .collect(Collectors.toList());

            if(candidateReads.isEmpty())
                continue;

            if(junction.DiscordantOnly)
            {
                if(mConfig.ProcessDiscordant)
                    discordantReads.processReads(junction, candidateReads);

                continue;
            }

            List<JunctionAssembly> candidateAssemblies = null;

            try
            {
                candidateAssemblies = junctionAssembler.processJunction(candidateReads);
            }
            catch(Exception e)
            {
                SV_LOGGER.error("failed to process junction({}) with {} reads", junction, candidateReads.size());
                e.printStackTrace();
                System.exit(1);
            }

            // dedup assemblies with close junction positions, same orientation
            dedupProximateAssemblies(junctionGroupAssemblies, candidateAssemblies);

            // extend assemblies with non-junction and discordant reads
            for(JunctionAssembly assembly : candidateAssemblies)
            {
                if(mDecoyChecker.enabled())
                {
                    if(mDecoyChecker.matchesDecoy(assembly))
                    {
                        SV_LOGGER.trace("assembly({}) matches decoy, excluding", assembly);
                        ++mReadStats.DecoySequences;
                        continue;
                    }
                }

                refBaseExtender.findAssemblyCandidateExtensions(assembly, junctionAssembler.nonJunctionReads());
                junctionGroupAssemblies.add(assembly);
            }
        }

        if(!discordantReads.groups().isEmpty())
        {
            discordantReads.mergeGroups();
            junctionGroup.addDiscordantGroups(discordantReads.groups());
        }

        junctionGroup.addJunctionAssemblies(junctionGroupAssemblies);

        // no longer needs to keep candidate reads since all have been assigned to assemblies
        junctionGroup.clearCandidateReads();

        // clear assembly read info for supporting fragments
        junctionGroupAssemblies.forEach(x -> x.clearSupportCachedReads());

        mCurrentJunctionGroup = null;
    }

    private void processRecord(final SAMRecord record)
    {
        mConfig.logReadId(record, "JunctionGroupAssembler:processRecord");

        // temporary checking of repeated (ie identical) supplementaries from SvPrep
        if(ignoreIdenticalSupplementary(record))
            return;

        Read read = new Read(record);

        ++mReadStats.TotalReads;

        if(mBamReader.currentIsReferenceSample())
            read.markReference();

        if(ReadAdjustments.trimPolyGSequences(read))
            ++mReadStats.PolyGTrimmed;

        if(ReadAdjustments.trimLowQualBases(read))
            ++mReadStats.LowBaseQualTrimmed;

        if(ReadAdjustments.convertEdgeIndelsToSoftClip(read))
            ++mReadStats.IndelSoftClipConverted;

        mCurrentJunctionGroup.addCandidateRead(read);

        ReadGroup readGroup = mReadGroupMap.get(read.id());

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
        mReadGroupMap.put(read.id(), readGroup);
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

        public String toString() { return format("reads(%d exp=%d) id(%s)", mReads.size(), mExpectedCount, mReads.get(0).id()); }
    }

    private boolean ignoreIdenticalSupplementary(final SAMRecord read)
    {
        if(!read.getSupplementaryAlignmentFlag())
            return false;

        SAMRecord previousRead = mSupplementaryRepeats.get(read.getReadName());

        if(previousRead == null)
        {
            mSupplementaryRepeats.put(read.getReadName(), read);
            return false;
        }

        if(previousRead.getAlignmentStart() != read.getAlignmentStart() || previousRead.getFlags() != read.getFlags())
            return false;

        if(!previousRead.getCigarString().equals(read.getCigarString()))
            return false;

        SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read);

        SV_LOGGER.trace("repeated supp({}) coords({}:{}) primary({}:{}) mate({}:{})",
                read.getReadName(), read.getReferenceName(), read.getAlignmentStart(),
                suppData.Chromosome, suppData.Position, read.getMateReferenceName(), read.getMateAlignmentStart());

        ++mReadStats.IdenticalSupplementaries;
        mSupplementaryRepeats.remove(read.getReadName());
        return true;
    }
}
