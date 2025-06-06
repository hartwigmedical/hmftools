package com.hartwig.hmftools.esvee.assembly;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.BAM_READ_JUNCTION_BUFFER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyDeduper.dedupProximateAssemblies;
import static com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments.markLineSoftClips;

import static htsjdk.samtools.CigarOperator.M;

import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.perf.TaskQueue;
import com.hartwig.hmftools.esvee.assembly.alignment.AlignmentChecker;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionGroup;
import com.hartwig.hmftools.esvee.assembly.types.ThreadTask;
import com.hartwig.hmftools.esvee.assembly.output.ResultsWriter;
import com.hartwig.hmftools.esvee.assembly.read.BamReader;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments;
import com.hartwig.hmftools.esvee.assembly.read.ReadStats;
import com.hartwig.hmftools.esvee.prep.ReadFilters;

import htsjdk.samtools.SAMRecord;

public class JunctionGroupAssembler extends ThreadTask
{
    private final AssemblyConfig mConfig;

    private final TaskQueue mJunctionGroups;

    private JunctionGroup mCurrentJunctionGroup;
    private final BamReader mBamReader;
    private final AlignmentChecker mAlignmentChecker;

    private final Map<String,ReadGroup> mReadGroupMap;
    private final Map<String,SAMRecord> mSupplementaryRepeats; // temporary to track an issue in SvPrep
    private final ReadStats mReadStats;
    private final List<JunctionAssembly> mDecoyAssemblies;

    public JunctionGroupAssembler(
            final AssemblyConfig config, final BamReader bamReader, final TaskQueue junctionGroups, final ResultsWriter resultsWriter)
    {
        super("PrimaryAssembly");
        mConfig = config;
        mBamReader = bamReader;
        mJunctionGroups = junctionGroups;

        mDecoyAssemblies = Lists.newArrayList();

        mAlignmentChecker = new AlignmentChecker(mConfig, resultsWriter.decoyMatchWriter());

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

        TaskQueue taskQueue = new TaskQueue(junctionGroupQueue, "junction groups", 10000);

        int junctionGroupCount = junctionGroups.size();

        for(int i = 0; i < taskCount; ++i)
        {
            BamReader bamReader = bamReaders.get(i);

            JunctionGroupAssembler junctionGroupAssembler = new JunctionGroupAssembler(config, bamReader, taskQueue, resultsWriter);
            primaryAssemblyTasks.add(junctionGroupAssembler);
            threadTasks.add(junctionGroupAssembler);
        }

        if(taskCount > 1)
        {
            SV_LOGGER.debug("primary assembly splits {} junction groups across {} threads", junctionGroupCount, taskCount);
        }

        return primaryAssemblyTasks;
    }

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                JunctionGroup junctionGroup = (JunctionGroup)mJunctionGroups.removeItem();

                mPerfCounter.start();
                processJunctionGroup(junctionGroup);

                stopCheckLog(format("juncGroup(%s)", junctionGroup), mConfig.PerfLogTime);
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
    public List<JunctionAssembly> decoyAssemblies() { return mDecoyAssemblies; }

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
        List<JunctionAssembly> dedupedIndels = Lists.newArrayList();

        RefBaseExtender refBaseExtender = new RefBaseExtender();

        // now pass applicable reads to each junction assembler - any read overlapping the junction
        // due to SvPrep filtering, most reads crossing the junction will have met soft-clip criteria
        for(int i = 0; i < junctionGroup.junctions().size(); ++i)
        {
            Junction junction = junctionGroup.junctions().get(i);

            JunctionAssembler junctionAssembler = new JunctionAssembler(junction, mConfig.RefGenome);

            // doesn't seem to be making a big difference, but this is inefficient for long-range junction groups
            // since both the junctions and reads are ordered. Could consider re-ordering by unclipped start and comparing to junction position

            int junctionBoundaryStart = junction.isForward() ? junction.Position - BAM_READ_JUNCTION_BUFFER : junction.Position;
            int junctionBoundaryEnd = junction.isForward() ? junction.Position : junction.Position + BAM_READ_JUNCTION_BUFFER;

            List<Read> candidateReads = junctionGroup.candidateReads().stream()
                    .filter(x -> positionsOverlap(junctionBoundaryStart, junctionBoundaryEnd, x.unclippedStart(), x.unclippedEnd()))
                    .collect(Collectors.toList());

            if(candidateReads.isEmpty())
                continue;

            if(junction.DiscordantOnly && mConfig.DiscordantOnlyDisabled)
                continue;

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
            dedupProximateAssemblies(junctionGroupAssemblies, candidateAssemblies, dedupedIndels);

            // extend assemblies with non-junction and discordant reads
            for(JunctionAssembly assembly : candidateAssemblies)
            {
                if(mAlignmentChecker.matchesDecoy(assembly))
                {
                    SV_LOGGER.trace("assembly({}) matches decoy, excluding", assembly);
                    ++mReadStats.DecoySequences;

                    mDecoyAssemblies.add(assembly);
                    continue;
                }

                if(mAlignmentChecker.failsMappability(assembly))
                {
                    SV_LOGGER.trace("assembly({}) fails ref-base alignment, excluding", assembly);
                    ++mReadStats.RefBaseAlignmentFails;
                    continue;
                }

                refBaseExtender.findAssemblyCandidateExtensions(assembly, junctionAssembler.nonJunctionReads());
                junctionGroupAssemblies.add(assembly);
            }
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
        // mConfig.logReadId(record, "JunctionGroupAssembler:processRecord");

        // temporary checking of repeated (ie identical) supplementaries from SvPrep
        if(ignoreIdenticalSupplementary(record))
            return;

        // old samples can be have invalid CIGARs
        if(!record.getReadUnmappedFlag() && record.getCigar().getCigarElements().stream().noneMatch(x -> x.getOperator() == M))
            return;

        if(ReadFilters.filterLowQualRead(record))
        {
            ++mReadStats.LowBaseQualFiltered;
            return;
        }

        Read read = new Read(record);

        ++mReadStats.TotalReads;

        if(mBamReader.currentIsReferenceSample())
            read.markReference();

        if(ReadAdjustments.trimPolyGSequences(read))
            ++mReadStats.PolyGTrimmed;

        markLineSoftClips(read);

        if(ReadAdjustments.trimLowQualSoftClipBases(read))
            ++mReadStats.LowBaseQualTrimmed;

        if(IndelBuilder.calcIndelInferredUnclippedPositions(read))
            ++mReadStats.IndelSoftClipConverted;

        mCurrentJunctionGroup.addCandidateRead(read);

        ReadGroup readGroup = mReadGroupMap.get(read.id());

        if(readGroup != null)
        {
            readGroup.addRead(read);
            return;
        }

        // link first and second in pair if within the same group
        boolean hasLocalMate = read.isMateUnmapped()
                || (read.isMateMapped() && read.mateChromosome().equals(read.chromosome())
                && positionsOverlap(
                        read.mateAlignmentStart(), read.mateAlignmentEnd(),
                        mCurrentJunctionGroup.readRangeStart(), mCurrentJunctionGroup.readRangeEnd()));

        // link first and second in pair if within the same group
        boolean hasLocalSupplementary = read.hasSupplementary() && read.supplementaryData().Chromosome.equals(read.chromosome())
                && positionWithin(
                        read.supplementaryData().Position, mCurrentJunctionGroup.readRangeStart(), mCurrentJunctionGroup.readRangeEnd());

        if(!hasLocalMate && !hasLocalSupplementary)
            return;

        int expectedCount = 1 + (hasLocalMate ? 1 : 0) + (hasLocalSupplementary ? 1 : 0); // approximate only for array size
        readGroup = new ReadGroup(read, expectedCount);
        mReadGroupMap.put(read.id(), readGroup);

        if(read.isUnmapped())
        {
            // check for an inconsistent mate which hasn't been cached in this read group
            for(Read cachedRead : mCurrentJunctionGroup.candidateReads())
            {
                if(cachedRead != read && cachedRead.isMateMapped() && cachedRead.id().equals(read.id()) && !readGroup.hasRead(cachedRead))
                {
                    readGroup.addRead(cachedRead);
                    break;
                }
            }
        }
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

        public boolean hasRead(final Read read) { return mReads.contains(read); }

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
        if(!mConfig.DevDebug)
            return false;

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
