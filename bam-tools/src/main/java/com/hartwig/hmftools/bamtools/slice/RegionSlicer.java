package com.hartwig.hmftools.bamtools.slice;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.slice.SliceConfig.UNMAPPED_READS_DISABLED;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class RegionSlicer
{
    private final SliceConfig mConfig;

    public RegionSlicer(final ConfigBuilder configBuilder)
    {
        mConfig = new SliceConfig(configBuilder);
    }

    public void run()
    {
        if(!mConfig.isValid())
            System.exit(1);

        BT_LOGGER.info("starting BamSlicer");

        long startTimeMs = System.currentTimeMillis();

        SliceWriter sliceWriter = new SliceWriter(mConfig);
        ReadCache readCache = new ReadCache(mConfig);

        Queue<ChrBaseRegion> regionsQueue = new ConcurrentLinkedQueue<>();
        regionsQueue.addAll(mConfig.SliceRegions.Regions);

        List<Thread> threadTasks = Lists.newArrayList();
        List<RegionBamSlicer> regionBamSlicers = Lists.newArrayList();

        for(int i = 0; i < min(mConfig.SliceRegions.Regions.size(), mConfig.Threads); ++i)
        {
            RegionBamSlicer regionSlicer = new RegionBamSlicer(regionsQueue, mConfig, readCache, sliceWriter);
            regionBamSlicers.add(regionSlicer);
            threadTasks.add(regionSlicer);
        }

        BT_LOGGER.info("splitting {} regions across {} threads", regionsQueue.size(), mConfig.Threads);

        if(!runThreadTasks(threadTasks))
            System.exit(1);

        BT_LOGGER.info("initial slice complete");

        List<RemotePositions> remoteChrPositions = Lists.newArrayList();

        RefGenomeSource refGenome = loadRefGenome(mConfig.RefGenomeFile);

        for(Map.Entry<String,List<RemotePosition>> entry : readCache.chrRemotePositions().entrySet())
        {
            // no reason not to include MT and any alt contigs unless they are not defined in the ref genome dictionary
            String contig = entry.getKey();

            if(refGenome.refGenomeFile().getSequenceDictionary().getSequence(contig) == null)
            {
                BT_LOGGER.debug("skipping remote slice for contig({}) not in ref genome", contig);
                continue;
            }

            List<RemotePosition> remotePositions = entry.getValue();
            Collections.sort(remotePositions);

            remoteChrPositions.add(new RemotePositions(contig, remotePositions));
        }

        if(!remoteChrPositions.isEmpty())
        {
            Collections.sort(remoteChrPositions);

            List<RemoteReadSlicer> remoteReadSlicers = remoteChrPositions.stream()
                    .map(x -> new RemoteReadSlicer(x.Chromosome, x.Positions, mConfig, sliceWriter)).collect(Collectors.toList());

            List<Callable> callableTasks = remoteReadSlicers.stream().collect(Collectors.toList());

            if(!TaskExecutor.executeTasks(callableTasks, mConfig.Threads))
                System.exit(1);

            BT_LOGGER.info("secondary slice complete");
        }

        if(mConfig.MaxUnmappedReads != UNMAPPED_READS_DISABLED)
        {
            BT_LOGGER.info("slicing unmapped reads");

            sliceUnmappedReads(sliceWriter);

            BT_LOGGER.info("unmapped read slice complete");
        }

        sliceWriter.close();

        BT_LOGGER.info("Regions slice complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void sliceUnmappedReads(final SliceWriter sliceWriter)
    {
        SamReader samReader = !mConfig.RefGenomeFile.isEmpty() ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        long unmappedCount = 0;

        SAMRecordIterator iterator = samReader.queryUnmapped();

        while(iterator.hasNext())
        {
            sliceWriter.writeRead(iterator.next());
            ++unmappedCount;

            if(mConfig.MaxUnmappedReads > 0 && unmappedCount >= mConfig.MaxUnmappedReads)
                break;
        }
    }

    private class RemotePositions implements Comparable<RemotePositions>
    {
        public final String Chromosome;
        public final List<RemotePosition> Positions;

        public RemotePositions(final String chromosome, final List<RemotePosition> positions)
        {
            Chromosome = chromosome;
            Positions = positions;
        }

        public String toString() { return format("chr(%s) positions(%d)", Chromosome, Positions.size()); }

        @Override
        public int compareTo(final RemotePositions other)
        {
            if(!Chromosome.equals(other.Chromosome))
                return Chromosome.compareTo(other.Chromosome);

            if(Positions.size() != other.Positions.size())
                return Positions.size() > other.Positions.size() ? -1 : 1;

            return 0;
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        SliceConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        RegionSlicer regionSlicer = new RegionSlicer(configBuilder);
        regionSlicer.run();
    }
}
