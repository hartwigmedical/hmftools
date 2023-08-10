package com.hartwig.hmftools.amber;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.AmberConstants.BAM_MIN_GAP_START;
import static com.hartwig.hmftools.amber.AmberConstants.CRAM_MIN_GAP_START;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamEvidenceReader
{
    final AmberConfig mConfig;
    private final BaseDepthFactory mBaseDepthFactory;

    public BamEvidenceReader(final AmberConfig config)
    {
        mConfig = config;
        mBaseDepthFactory = new BaseDepthFactory(mConfig.MinBaseQuality);
    }

    private class RegionTask
    {
        public final ChrBaseRegion Region;

        private final List<BaseDepth> mPositions;
        private int mCurrentIndex;
        private boolean mComplete;

        public RegionTask(final String chromosome, final BaseDepth baseDepth)
        {
            Region = new ChrBaseRegion(chromosome, baseDepth.Position, baseDepth.Position);
            mPositions = Lists.newArrayList(baseDepth);
            mCurrentIndex = 0;
            mComplete = false;
        }

        public void addPosition(final BaseDepth baseDepth)
        {
            mPositions.add(baseDepth);
            Region.setEnd(max(Region.end(), baseDepth.Position));
        }

        public void processRecord(final SAMRecord record)
        {
            int alignmentStart = record.getAlignmentStart();
            int alignmentEnd = record.getAlignmentEnd();

            int index = mCurrentIndex;
            for(; index < mPositions.size(); ++index)
            {
                BaseDepth baseDepth = mPositions.get(index);

                if(alignmentStart > baseDepth.Position)
                {
                    ++mCurrentIndex;
                    continue;
                }

                if(alignmentEnd < baseDepth.Position)
                    break;

                mBaseDepthFactory.addEvidence(baseDepth, record);
            }

            if(mCurrentIndex >= mPositions.size())
                mComplete = true;
        }

        public boolean isComplete() { return mComplete; }

        public String toString() { return format("region(%s) positions(%d) index(%d)", Region, mPositions.size(), mCurrentIndex); }
    }

    private class BamReaderThread extends Thread
    {
        final Queue<RegionTask> mTaskQueue;
        final SamReader mSamReader;
        final BamSlicer mBamSlicer;
        RegionTask mCurrentTask;

        BamReaderThread(
                final String bamFile, final SamReaderFactory samReaderFactory, final Queue<RegionTask> inTaskQueue,
                int minMappingQuality)
        {
            mTaskQueue = inTaskQueue;
            mSamReader = samReaderFactory.open(new File(bamFile));
            mBamSlicer = new BamSlicer(minMappingQuality, false, false, false);
            mCurrentTask = null;
        }

        @Override
        public void run()
        {
            // AMB_LOGGER.debug("bam reader thread start");

            while(true)
            {
                RegionTask task;
                try
                {
                    task = mTaskQueue.remove();
                    mCurrentTask = task;
                }
                catch (NoSuchElementException e)
                {
                    // finished processing
                    break;
                }

                mBamSlicer.slice(mSamReader, task.Region, this::processRecord);
            }

            try
            {
                mSamReader.close();
            }
            catch(IOException e)
            {
                AMB_LOGGER.error("IO exception in SamReader::close: {}", e.getMessage());
            }

            // AMB_LOGGER.debug("bam reader thread finish");
        }

        private void processRecord(final SAMRecord record)
        {
            if(mCurrentTask == null)
            {
                mBamSlicer.haltProcessing();
                return;
            }

            mCurrentTask.processRecord(record);

            if(mCurrentTask.isComplete())
                mBamSlicer.haltProcessing();
        }
    }

    public void processBam(
            final String bamFile, final SamReaderFactory samReaderFactory, final Map<Chromosome,List<BaseDepth>> chrBaseDepth)
            throws InterruptedException
    {
        AMB_LOGGER.debug("processing bam({})", bamFile);

        final Queue<RegionTask> taskQueue = new ConcurrentLinkedQueue<>();

        // create genome regions from the loci
        boolean limitRegions = bamFile.endsWith(".cram");
        populateTaskQueue(chrBaseDepth, taskQueue, limitRegions);

        // we create the consumer and producer
        List<BamReaderThread> bamReaders = new ArrayList<BamReaderThread>();

        for(int i = 0; i < max(mConfig.Threads, 1); ++i)
        {
            BamReaderThread thread = new BamReaderThread(bamFile, samReaderFactory, taskQueue, mConfig.MinMappingQuality);
            thread.setName(format("worker-%d", i));
            thread.start();
            bamReaders.add(thread);
        }

        AMB_LOGGER.trace("{} bam reader threads started", bamReaders.size());

        AmberTaskCompletion taskCompletion = new AmberTaskCompletion(taskQueue.size());
        for(BamReaderThread thread : bamReaders)
        {
            while(thread.isAlive())
            {
                // check status every 30 seconds
                thread.join(30_000);

                // check status
                taskCompletion.progress(taskQueue.size());
            }
        }

        AMB_LOGGER.trace("{} bam reader threads finished", bamReaders.size());
    }

    private void populateTaskQueue(
            final Map<Chromosome,List<BaseDepth>> chrBaseDepth, final Queue<RegionTask> taskQueue, boolean limitRegions)
    {
        int positionCount = chrBaseDepth.values().stream().mapToInt(x -> x.size()).sum();

        int minGap = limitRegions ? CRAM_MIN_GAP_START : BAM_MIN_GAP_START;

        // int maxPositionsPerRegion = limitRegions ? CRAM_REGION_GROUP_MAX : BAM_REGION_GROUP_MAX;
        // int gapIncrement = limitRegions ? CRAM_MIN_GAP_INCREMENT : BAM_MIN_GAP_INCREMENT;

        List<RegionTask> tasks = Lists.newArrayList();

        for(Map.Entry<Chromosome,List<BaseDepth>> entry : chrBaseDepth.entrySet())
        {
            String chromosome = entry.getKey().toString();
            List<BaseDepth> positions = entry.getValue();

            if(positions.isEmpty())
                continue;

            RegionTask currentTask = new RegionTask(chromosome, positions.get(0));
            tasks.add(currentTask);

            for(int i = 1; i < positions.size(); ++i)
            {
                BaseDepth position = positions.get(i);

                if(currentTask.Region.end() + minGap < position.Position) // or  || tasks.size() >= maxPositionsPerRegion
                {
                    // start a new region
                    currentTask = new RegionTask(chromosome, position);
                    tasks.add(currentTask);
                }
                else
                {
                    currentTask.addPosition(position);
                }
            }
        }

        taskQueue.addAll(tasks);

        AMB_LOGGER.debug("split {} sites across {} regions", positionCount, taskQueue.size());
    }
}
