package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.SvConstants.BAM_READ_JUNCTION_BUFFER;
import static com.hartwig.hmftools.esvee.SvConstants.TASK_LOG_COUNT;
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
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.common.JunctionGroup;
import com.hartwig.hmftools.esvee.common.ThreadTask;
import com.hartwig.hmftools.esvee.output.ResultsWriter;
import com.hartwig.hmftools.esvee.read.BamReader;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.read.ReadAdjustments;
import com.hartwig.hmftools.esvee.read.ReadFilters;

import htsjdk.samtools.SAMRecord;

public class JunctionGroupAssembler extends ThreadTask
{
    private final SvConfig mConfig;
    private final ResultsWriter mResultsWriter;

    private final Queue<JunctionGroup> mJunctionGroups;
    private final int mJunctionCount;

    private JunctionGroup mCurrentJunctionGroup;
    private final BamReader mBamReader;

    private final Map<String,Read> mReadGroupMap;

    public JunctionGroupAssembler(
            final SvConfig config, final BamReader bamReader, final Queue<JunctionGroup> junctionGroups, final ResultsWriter resultsWriter)
    {
        super("PrimaryAssembly");
        mConfig = config;
        mResultsWriter = resultsWriter;
        mBamReader = bamReader;
        mJunctionGroups = junctionGroups;
        mJunctionCount = junctionGroups.size();

        mReadGroupMap = Maps.newHashMap();
        mCurrentJunctionGroup = null;
    }

    public static List<JunctionGroupAssembler> createThreadTasks(
            final List<JunctionGroup> junctionGroups, final List<BamReader> bamReaders, final SvConfig config,
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

        List<JunctionAssembly> junctionGroupAssemblies = Lists.newArrayList();

        // now pass applicable reads to each primary assembler
        for(int i = 0; i < mCurrentJunctionGroup.junctions().size(); ++i)
        {
            Junction junction = mCurrentJunctionGroup.junctions().get(i);

            JunctionAssembler junctionAssembler = new JunctionAssembler(mConfig, junction);

            // FIXME: doesn't seem to be making a big difference, but this is in efficient for long-range junction groups
            // since both the junctions and reads are ordered. Could consider re-ordering by unclipped start and comparing to junction position

            List<Read> junctionCandidateReads = mCurrentJunctionGroup.candidateReads().stream()
                    .filter(x -> ReadFilters.isCandidateJunctionRead(x, junction))
                    .collect(Collectors.toList());

            if(junctionCandidateReads.isEmpty())
                continue;

            List<JunctionAssembly> candidateAssemblies = junctionAssembler.processJunction(junctionCandidateReads);

            // dedup assemblies with close junction positions, same orientation
            dedupProximateAssemblies(junctionGroupAssemblies, candidateAssemblies);

            junctionGroupAssemblies.addAll(candidateAssemblies);

            // extend assemblies with non-junction and discordant reads
            // CHECK: skip likely germline reads here?
            for(JunctionAssembly assembly : candidateAssemblies)
            {
                AssemblyExtender assemblyExtender = new AssemblyExtender(assembly);
                assemblyExtender.extendAssembly(junctionAssembler.nonJunctionReads());
            }
        }

        mCurrentJunctionGroup.addJunctionAssemblies(junctionGroupAssemblies);

        junctionGroupAssemblies.forEach(x -> mResultsWriter.writeAssembly(x));

        mCurrentJunctionGroup = null;
        mReadGroupMap.clear();
    }

    private static final int MATE_READ_BUFFER = 200;

    private void processRecord(final SAMRecord record)
    {
        // CHECK: do in SvPrep if worthwhile
        if(!ReadFilters.isRecordAverageQualityAbove(record.getBaseQualities(), SvConstants.AVG_BASE_QUAL_THRESHOLD))
            return;

        Read read = new Read(record);

        if(mBamReader.currentIsReferenceSample())
            read.markReference();

        // CHECK: could track for stats
        ReadAdjustments.trimPolyGSequences(read);
        ReadAdjustments.convertEdgeIndelsToSoftClip(read);

        // mReadRescue::rescueRead) // CHECK: logic and purpose

        mCurrentJunctionGroup.addCandidateRead(read);

        // link first and second in pair if within the same group
        if(!read.isMateMapped() || !read.mateChromosome().equals(read.chromosome()))
            return;

        if(positionWithin(
                read.mateAlignmentStart(),
                mCurrentJunctionGroup.minPosition() - BAM_READ_JUNCTION_BUFFER,
                mCurrentJunctionGroup.maxPosition() + BAM_READ_JUNCTION_BUFFER))
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
}
