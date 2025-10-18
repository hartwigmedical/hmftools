package com.hartwig.hmftools.redux.unmap;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.redux.PartitionReader.fullyUnmapped;
import static com.hartwig.hmftools.redux.PartitionReader.shouldFilterRead;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ReduxConstants.UNMAP_MAX_NON_OVERLAPPING_BASES;
import static com.hartwig.hmftools.redux.common.ReadInfo.readToString;
import static com.hartwig.hmftools.redux.unmap.ReadUnmapper.overlapsUnmapRegion;
import static com.hartwig.hmftools.redux.unmap.UnmapRegion.UNMAPPED_READS;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.perf.TaskQueue;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.mappability.UnmappingRegion;
import com.hartwig.hmftools.redux.BamReader;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.consensus.UltimaRoutines;
import com.hartwig.hmftools.redux.write.BamWriterSync;
import com.hartwig.hmftools.redux.write.FileWriterCache;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class RegionUnmapper extends Thread
{
    private final ReduxConfig mConfig;
    private final ReadUnmapper mReadUnmapper;
    private final TaskQueue mRegions;
    private final BamReader mBamReader;
    private final BamWriterSync mUnmappingBamWriter;
    private final BamWriterSync mFullyUnmappedBamWriter;

    private UnmapRegion mCurrentRegion;
    private UnmapRegionState mUnmapRegionState;

    private long mProcessedReads;
    private long mNextLogReadCount;
    private final boolean mLogReadIds;

    public RegionUnmapper(
            final ReduxConfig config, final TaskQueue regions, final BamWriterSync unmappingBamWriter, final BamWriterSync fullUnmappedBamWriter)
    {
        mConfig = config;
        mReadUnmapper = mConfig.UnmapRegions;
        mUnmappingBamWriter = unmappingBamWriter;
        mFullyUnmappedBamWriter = fullUnmappedBamWriter;
        mRegions = regions;

        mBamReader = new BamReader(config.BamFiles, config.RefGenomeFile);
        mUnmapRegionState = null;
        mCurrentRegion = null;

        mProcessedReads = 0;
        mNextLogReadCount = LOG_READ_COUNT;
        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    private static final int REGION_PROXIMITY_BUFFER = 50;

    public static List<RegionUnmapper> createThreadTasks(
            final ReduxConfig config, final FileWriterCache fileWriterCache, final List<Thread> threadTasks)
    {
        Map<String,List<UnmappingRegion>> allUnmappingRegions = config.UnmapRegions.getAllRegions();

        List<UnmapRegion> unmappingRegions = Lists.newArrayList();

        for(Map.Entry<String,List<UnmappingRegion>> entry : allUnmappingRegions.entrySet())
        {
            String chromosome = entry.getKey();

            if(config.SpecificChrRegions.excludeChromosome(chromosome))
                continue;

            List<UnmappingRegion> regions = entry.getValue();

            UnmapRegion currentRegion = null;

            // merge regions which are close together after factoring in the addition buffer
            for(UnmappingRegion region : regions)
            {
                if(config.SpecificChrRegions.excludeRegion(region.start(), region.end()))
                    continue;

                // build a buffer into the slice regions to account for unmapped mates starting earlier
                int regionStart = max(region.start() - (UNMAP_MAX_NON_OVERLAPPING_BASES * 2), 1);

                if(currentRegion == null || currentRegion.end() < regionStart - REGION_PROXIMITY_BUFFER)
                {
                    currentRegion = new UnmapRegion(chromosome, regionStart, region.end(), region);
                    unmappingRegions.add(currentRegion);
                }
                else
                {
                    currentRegion.setEnd(region.end());
                    currentRegion.Regions.add(region);
                }
            }
        }

        Collections.sort(unmappingRegions);

        if(config.SpecificChrRegions.hasFilters() && unmappingRegions.isEmpty())
            return Collections.emptyList();

        // add in an entry to extract fully unmapped reads
        if(!config.SkipFullyUnmappedReads)
            unmappingRegions.add(UNMAPPED_READS);

        List<RegionUnmapper> regionUnmapperTasks = Lists.newArrayList();

        Queue<ChrBaseRegion> regionsQueue = new ConcurrentLinkedQueue<>();
        regionsQueue.addAll(unmappingRegions);

        TaskQueue taskQueue = new TaskQueue(regionsQueue, "unmap regions", 10000);

        BamWriterSync unmappingBamWriter = fileWriterCache.getUnmappingBamWriter();
        BamWriterSync fullUnmappedBamWriter = fileWriterCache.getFullUnmappedBamWriter();

        for(int i = 0; i < config.Threads; ++i)
        {
            RegionUnmapper regionUnmapper = new RegionUnmapper(config, taskQueue, unmappingBamWriter, fullUnmappedBamWriter);
            regionUnmapperTasks.add(regionUnmapper);
            threadTasks.add(regionUnmapper);
        }

        if(config.Threads > 1)
        {
            RD_LOGGER.debug("splitting {} unmapping regions across {} threads", unmappingRegions.size(), config.Threads);
        }

        return regionUnmapperTasks;
    }

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                UnmapRegion unmapRegion = (UnmapRegion)mRegions.removeItem();

                processRegion(unmapRegion);
            }
            catch(NoSuchElementException e)
            {
                RD_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    public long processedReads() { return mProcessedReads; }

    private void processRegion(final UnmapRegion region)
    {
        if(region == UNMAPPED_READS)
        {
            RD_LOGGER.trace("extracting fully-unmapped region");
            processFullyUnmappedReads(mConfig, mFullyUnmappedBamWriter);
            return;
        }

        mCurrentRegion = region;
        RD_LOGGER.trace("unmapping region({})", region);

        mUnmapRegionState = new UnmapRegionState(mCurrentRegion, mCurrentRegion.Regions);
        mUnmapRegionState.LastMatchedIndex = 0;

        if(mBamReader != null)
        {
            mBamReader.sliceRegion(mCurrentRegion, this::processSamRecord);
        }
    }

    private static final int LOG_READ_COUNT = 1000000;

    private void processSamRecord(final SAMRecord read)
    {
        if(shouldFilterRead(read))
            return;

        // must overlap by minimum to be a candidate for unmapping
        int readStart = read.getAlignmentStart();

        if(mLogReadIds && mConfig.LogReadIds.contains(read.getReadName())) // debugging only
        {
            RD_LOGGER.debug("specific read: {}", readToString(read));
        }

        if(!read.getReadUnmappedFlag()) // early exit for a read which can't be unmapped
        {
            if(mCurrentRegion.Regions.stream().noneMatch(x -> overlapsUnmapRegion(x, read.getAlignmentStart(), read.getAlignmentEnd())))
                return;
        }

        ++mProcessedReads;

        if(mProcessedReads >= mNextLogReadCount)
        {
            double processedReads = mProcessedReads / 1000000.0;
            RD_LOGGER.debug("region({}) position({}) processed {}M reads",
                    mCurrentRegion, readStart, format("%.0f", processedReads));

            mNextLogReadCount += LOG_READ_COUNT;
        }

        read.setDuplicateReadFlag(false);

        boolean alreadyUnmapped = read.getReadUnmappedFlag();
        mReadUnmapper.checkTransformRead(read, mUnmapRegionState);

        boolean internallyUnmapped = !alreadyUnmapped && read.getReadUnmappedFlag();
        boolean fullyUnmapped = fullyUnmapped(read);

        // relocate any read which satisfies unmapping criteria - either to its mate's location or the fully-unmapped BAM
        if(internallyUnmapped || fullyUnmapped)
        {
            if(read.getSupplementaryAlignmentFlag() || read.isSecondaryAlignment())
                return;

            if(mConfig.RunChecks)
                mConfig.readChecker().addUnmappedRead(read, mCurrentRegion.Chromosome, readStart);

            if(fullyUnmapped)
                mFullyUnmappedBamWriter.writeRecordSync(read);
            else
                mUnmappingBamWriter.writeRecordSync(read);
        }
    }

    public static void processFullyUnmappedReads(final ReduxConfig config, final BamWriterSync fullyUnmappedBamWriter)
    {
        int totalUnmappedReads = 0;
        boolean isUltima = ReduxConfig.isUltima();

        for(String bamFilename : config.BamFiles)
        {
            SamReader samReader = SamReaderFactory.makeDefault()
                    .referenceSequence(new File(config.RefGenomeFile)).open(new File(bamFilename));

            SAMRecordIterator iterator = samReader.queryUnmapped();

            while(iterator.hasNext())
            {
                SAMRecord record = iterator.next();

                if(isUltima)
                    UltimaRoutines.stripAttributes(record);

                fullyUnmappedBamWriter.writeRecordSync(record);
                ++totalUnmappedReads;
            }
        }

        if(totalUnmappedReads > 0)
        {
            RD_LOGGER.debug("extracted {} fully-unmapped reads", totalUnmappedReads);
            config.UnmapRegions.stats().ExistingUnmapped.addAndGet(totalUnmappedReads);
        }
    }
}
