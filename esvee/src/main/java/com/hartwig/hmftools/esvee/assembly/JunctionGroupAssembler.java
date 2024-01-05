package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.BAM_READ_JUNCTION_BUFFER;
import static com.hartwig.hmftools.esvee.SvConstants.TASK_LOG_COUNT;

import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.common.JunctionGroup;
import com.hartwig.hmftools.esvee.common.ThreadTask;
import com.hartwig.hmftools.esvee.read.BamReader;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.sequence.PrimaryAssembly;

import htsjdk.samtools.SAMRecord;

public class JunctionGroupAssembler extends ThreadTask
{
    private final SvConfig mConfig;

    private final Queue<JunctionGroup> mJunctionGroups;
    private final int mJunctionCount;

    private JunctionGroup mCurrentJunctionGroup;
    private final BamReader mBamReader;

    private final Map<String,Read> mReadGroupMap;

    private final List<PrimaryAssembly> mPrimaryAssemblies;

    public JunctionGroupAssembler(final SvConfig config, final BamReader bamReader, final Queue<JunctionGroup> junctionGroups)
    {
        super("PrimaryAssembly");
        mConfig = config;
        mBamReader = bamReader;
        mJunctionGroups = junctionGroups;
        mJunctionCount = junctionGroups.size();

        mReadGroupMap = Maps.newHashMap();
        mCurrentJunctionGroup = null;
        mPrimaryAssemblies = Lists.newArrayList();
    }

    public List<PrimaryAssembly> primaryAssemblies() { return mPrimaryAssemblies; }

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
                    SV_LOGGER.info("processed {} junction groups, remaining({})", processedCount, remainingCount);
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

    private void processJunctionGroup(final JunctionGroup junctionGroup)
    {
        mCurrentJunctionGroup = junctionGroup;

        int sliceStart = mCurrentJunctionGroup.minPosition() - BAM_READ_JUNCTION_BUFFER;
        int sliceEnd = mCurrentJunctionGroup.maxPosition() + BAM_READ_JUNCTION_BUFFER;

        SV_LOGGER.trace("junctionGroup({}:{}-{} count={}) slice",
                mCurrentJunctionGroup.chromosome(), sliceStart, sliceEnd, mCurrentJunctionGroup.count());

        mBamReader.sliceBam(mCurrentJunctionGroup.chromosome(), sliceStart, sliceEnd, this::processRecord);

        SV_LOGGER.debug("junctionGroup({}:{}-{} count={}) slice complete, readCount({})",
                mCurrentJunctionGroup.chromosome(), sliceStart, sliceEnd, mCurrentJunctionGroup.count(), mCurrentJunctionGroup.candidateReadCount());

        // now pass applicable reads to each primary assembler
        for(Junction junction : mCurrentJunctionGroup.junctions())
        {
            PrimaryAssembler primaryAssembler = new PrimaryAssembler(mConfig, junction);

            List<Read> junctionCandidateReads = mCurrentJunctionGroup.candidateReads().stream()
                    .filter(x -> AlignmentFilters.alignmentCrossesJunction(x, junction))
                    .collect(Collectors.toList());

            if(junctionCandidateReads.isEmpty())
                continue;

            List<PrimaryAssembly> candidateAssemblies = primaryAssembler.processJunction(junctionCandidateReads);
            mPrimaryAssemblies.addAll(candidateAssemblies);

            // mPrimaryAssemblyResults.add(new PrimaryAssemblyResult(junction, primaryAssembler.getCounters(), candidateAssemblies));
        }

        mCurrentJunctionGroup = null;
        mReadGroupMap.clear();
    }

    private static final int MATE_READ_BUFFER = 200;

    private void processRecord(final SAMRecord record)
    {
        // CHECK: do in SvPrep if worthwhile
        if(!AlignmentFilters.isRecordAverageQualityAbove(record.getBaseQualities(), SvConstants.AVG_BASE_QUAL_THRESHOLD))
            return;

        // mReadRescue::rescueRead) // CHECK: not required

        // mNormaliser::normalise() on each record

        Read read = new Read(record);

        mCurrentJunctionGroup.addCandidateRead(read);

        // link first and second in pair if within the same group
        if(read.isMateMapped() && read.mateChromosome().equals(read.getChromosome())
        && positionWithin(read.mateAlignmentStart(), mCurrentJunctionGroup.minPosition() - MATE_READ_BUFFER, mCurrentJunctionGroup.maxPosition()))
        {
            Read mateRead = mReadGroupMap.get(read.getName());

            if(mateRead != null)
            {
                mReadGroupMap.remove(read.getName());
                mateRead.setMateRead(read);
                read.setMateRead(mateRead);
            }
            else
            {
                mReadGroupMap.put(read.getName(), read);
            }
        }
    }

    public static List<PrimaryAssembly> mergePrimaryAssemblies(final List<JunctionGroupAssembler> assemblers)
    {
        List<PrimaryAssembly> primaryAssemblies = Lists.newArrayList();
        assemblers.forEach(x -> primaryAssemblies.addAll(x.primaryAssemblies()));
        return primaryAssemblies;
    }
}
