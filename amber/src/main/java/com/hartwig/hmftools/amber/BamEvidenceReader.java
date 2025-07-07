package com.hartwig.hmftools.amber;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.AmberConstants.BAM_MIN_GAP_START;
import static com.hartwig.hmftools.amber.AmberConstants.CRAM_MIN_GAP_START;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.perf.PerformanceCounter;

import htsjdk.samtools.SamReaderFactory;

public class BamEvidenceReader
{
    private final AmberConfig mConfig;
    private final PositionEvidenceChecker mEvidenceChecker;

    public BamEvidenceReader(final AmberConfig config)
    {
        mConfig = config;
        mEvidenceChecker = new PositionEvidenceChecker(mConfig.MinBaseQuality);
    }

    public void processBam(
            final String bamFile, final SamReaderFactory samReaderFactory, final Map<Chromosome,List<PositionEvidence>> chrPositionEvidence)
            throws InterruptedException
    {
        AMB_LOGGER.trace("processing bam({})", bamFile);

        final Queue<RegionTask> taskQueue = new ConcurrentLinkedQueue<>();

        // create genome regions from the loci
        boolean limitRegions = bamFile.endsWith(".cram");
        populateTaskQueue(chrPositionEvidence, taskQueue, limitRegions);

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

        ProgressTracker taskCompletion = new ProgressTracker(taskQueue.size());
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

        if(AMB_LOGGER.isDebugEnabled())
        {
            PerformanceCounter combinedPc = new PerformanceCounter("Read");
            bamReaders.forEach(x -> combinedPc.merge(x.perfCounter()));
            combinedPc.logStats();
        }
    }

    private void populateTaskQueue(
            final Map<Chromosome,List<PositionEvidence>> chrBaseDepth, final Queue<RegionTask> taskQueue, boolean limitRegions)
    {
        int positionCount = chrBaseDepth.values().stream().mapToInt(x -> x.size()).sum();

        int minGap = mConfig.PositionGap > 0 ? mConfig.PositionGap : (limitRegions ? CRAM_MIN_GAP_START : BAM_MIN_GAP_START);

        // int maxPositionsPerRegion = limitRegions ? CRAM_REGION_GROUP_MAX : BAM_REGION_GROUP_MAX;
        // int gapIncrement = limitRegions ? CRAM_MIN_GAP_INCREMENT : BAM_MIN_GAP_INCREMENT;

        List<RegionTask> tasks = Lists.newArrayList();

        for(Map.Entry<Chromosome,List<PositionEvidence>> entry : chrBaseDepth.entrySet())
        {
            String chromosome = mConfig.RefGenVersion.versionedChromosome(entry.getKey().toString());

            if(!mConfig.SpecificChromosomes.isEmpty() && !mConfig.SpecificChromosomes.contains(chromosome))
                continue;

            List<PositionEvidence> positions = entry.getValue();

            if(positions.isEmpty())
                continue;

            RegionTask currentTask = new RegionTask(mEvidenceChecker, chromosome, positions.get(0));
            tasks.add(currentTask);

            for(int i = 1; i < positions.size(); ++i)
            {
                PositionEvidence posEvidence = positions.get(i);

                if(currentTask.Region.end() + minGap < posEvidence.Position) // or  || tasks.size() >= maxPositionsPerRegion
                {
                    // start a new region
                    currentTask = new RegionTask(mEvidenceChecker, chromosome, posEvidence);
                    tasks.add(currentTask);
                }
                else
                {
                    currentTask.addPosition(posEvidence);
                }
            }
        }

        if(AMB_LOGGER.isDebugEnabled())
        {
            RegionTask maxRegion = null;

            for(RegionTask task : tasks)
            {
                if(maxRegion == null || task.Region.length() > maxRegion.Region.length() || task.positionCount() > maxRegion.positionCount())
                    maxRegion = task;
            }

            AMB_LOGGER.debug("split {} sites across {} regions, max region({} size={} length={}) minGap({})",
                    positionCount, tasks.size(), maxRegion.Region, maxRegion.positionCount(), maxRegion.Region.length(), minGap);
        }

        taskQueue.addAll(tasks);
    }
}
