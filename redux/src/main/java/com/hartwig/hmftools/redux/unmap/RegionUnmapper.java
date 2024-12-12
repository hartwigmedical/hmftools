package com.hartwig.hmftools.redux.unmap;

import static java.lang.String.format;

import static com.hartwig.hmftools.redux.PartitionReader.fullyUnmapped;
import static com.hartwig.hmftools.redux.PartitionReader.shouldFilterRead;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.Constants.UNMAP_MAX_NON_OVERLAPPING_BASES;
import static com.hartwig.hmftools.redux.common.ReadInfo.readToString;
import static com.hartwig.hmftools.redux.unmap.ReadUnmapper.overlapsRegion;
import static com.hartwig.hmftools.redux.unmap.UnmapRegion.UNMAPPED_READS;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.HighDepthRegion;
import com.hartwig.hmftools.redux.BamReader;
import com.hartwig.hmftools.redux.ReduxConfig;
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

    public static List<RegionUnmapper> createThreadTasks(
            final ReduxConfig config, final FileWriterCache fileWriterCache, final List<Thread> threadTasks)
    {
        Map<String,List<HighDepthRegion>> allUnmappingRegions = config.UnmapRegions.getAllRegions();

        List<UnmapRegion> unmappingRegions = Lists.newArrayList();

        for(Map.Entry<String,List<HighDepthRegion>> entry : allUnmappingRegions.entrySet())
        {
            String chromosome = entry.getKey();

            if(config.SpecificChrRegions.excludeChromosome(chromosome))
                continue;

            List<HighDepthRegion> regions = entry.getValue();

            UnmapRegion currentRegion = null;

            // merge regions which are close together after factoring in the addition buffer
            for(HighDepthRegion region : regions)
            {
                if(config.SpecificChrRegions.excludeRegion(region.start(), region.end()))
                    continue;

                // build a buffer into the slice regions to account for unmapped mates starting earlier
                int regionStart = region.start() - (UNMAP_MAX_NON_OVERLAPPING_BASES * 2);
                int regionEnd = region.end();

                if(currentRegion == null || currentRegion.end() < regionStart - 50)
                {
                    currentRegion = new UnmapRegion(chromosome, regionStart, regionEnd, region);
                    unmappingRegions.add(currentRegion);
                }
                else
                {
                    currentRegion.setEnd(regionEnd);
                }
            }
        }

        Collections.sort(unmappingRegions);

        // add in an entry to extract fully unmapped reads
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
            processFullyUnmappedReads();
            return;
        }

        mCurrentRegion = region;

        mUnmapRegionState = new UnmapRegionState(mCurrentRegion, List.of(mCurrentRegion.Region));
        mUnmapRegionState.LastMatchedRegionIndex = 0;

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

        if(mCurrentRegion.IsStandardChromosome)
        {
            if(!read.getReadUnmappedFlag() && !overlapsRegion(mCurrentRegion.Region, read.getAlignmentStart(), read.getAlignmentEnd()))
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
        int readFlags = read.getFlags();
        mReadUnmapper.checkTransformRead(read, mUnmapRegionState);

        boolean internallyUnmapped = !alreadyUnmapped && read.getReadUnmappedFlag();
        boolean fullyUnmapped = fullyUnmapped(read);

        // relocate any read which satisfies unmapping criteria - either to its mate's location or the fully-unmapped BAM
        if(internallyUnmapped || fullyUnmapped)
        {
            if(read.getSupplementaryAlignmentFlag() || read.isSecondaryAlignment())
                return;

            if(mConfig.RunChecks)
                mReadUnmapper.addUnmappedRead(read, mCurrentRegion.Chromosome, readStart, readFlags);

            if(fullyUnmapped)
                mFullyUnmappedBamWriter.writeRecordSync(read);
            else
                mUnmappingBamWriter.writeRecordSync(read);
        }
    }

    private boolean processFullyUnmappedReads()
    {
        int totalUnmappedReads = 0;

        for(String bamFilename : mConfig.BamFiles)
        {
            SamReader samReader = SamReaderFactory.makeDefault()
                    .referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(bamFilename));

            SAMRecordIterator iterator = samReader.queryUnmapped();

            while(iterator.hasNext())
            {
                SAMRecord record = iterator.next();
                mFullyUnmappedBamWriter.writeRecordSync(record);
                ++totalUnmappedReads;
            }
        }

        if(totalUnmappedReads > 0)
        {
            RD_LOGGER.debug("extracted {} fully-unmapped reads", totalUnmappedReads);
            mReadUnmapper.stats().ExistingUnmapped.addAndGet(totalUnmappedReads);
        }

        return true;
    }
}
