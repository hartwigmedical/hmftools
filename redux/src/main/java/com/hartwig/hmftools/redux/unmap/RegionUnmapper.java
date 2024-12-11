package com.hartwig.hmftools.redux.unmap;

import static java.lang.String.format;

import static com.hartwig.hmftools.redux.PartitionReader.fullyUnmapped;
import static com.hartwig.hmftools.redux.PartitionReader.shouldFilterRead;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
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
    private final TaskQueue mRegions;
    private final BamReader mBamReader;
    private final BamWriterSync mUnmappingBamWriter;
    private final BamWriterSync mFullyUnmappedBamWriter;

    private UnmapRegion mCurrentRegion;
    private HighDepthRegion mHighDepthRegion;
    private UnmapRegionState mUnmapRegionState;

    private long mProcessedReads;
    private long mNextLogReadCount;
    private final boolean mLogReadIds;

    public RegionUnmapper(
            final ReduxConfig config, final TaskQueue regions, final BamWriterSync unmappingBamWriter, final BamWriterSync fullUnmappedBamWriter)
    {
        mConfig = config;
        mUnmappingBamWriter = unmappingBamWriter;
        mFullyUnmappedBamWriter = fullUnmappedBamWriter;
        mRegions = regions;

        mBamReader = new BamReader(config.BamFiles, config.RefGenomeFile);
        mUnmapRegionState = null;
        mCurrentRegion = null;
        mHighDepthRegion = null;

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

            for(HighDepthRegion region : regions)
            {
                if(config.SpecificChrRegions.excludeRegion(region.start(), region.end()))
                    continue;

                unmappingRegions.add(new UnmapRegion(chromosome, region.start(), region.end(), region.maxDepth()));
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
        mHighDepthRegion = new HighDepthRegion(region.start(), region.end(), region.MaxDepth);

        mUnmapRegionState = new UnmapRegionState(mCurrentRegion, List.of(mHighDepthRegion));
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
        if(mCurrentRegion.IsStandardChromosome)
        {
            int readPosEnd = read.getReadUnmappedFlag() ? read.getAlignmentStart() : read.getAlignmentEnd();

            if(!overlapsRegion(mHighDepthRegion, read.getAlignmentStart(), readPosEnd))
                return;
        }

        ++mProcessedReads;

        if(mProcessedReads >= mNextLogReadCount)
        {
            double processedReads = mProcessedReads / 1000000.0;
            RD_LOGGER.debug("region({}) position({}) processed {}M reads",
                    mCurrentRegion, read.getAlignmentStart(), format("%.0f", processedReads));

            mNextLogReadCount += LOG_READ_COUNT;
        }

        read.setDuplicateReadFlag(false);

        if(mLogReadIds && mConfig.LogReadIds.contains(read.getReadName())) // debugging only
        {
            RD_LOGGER.debug("specific read: {}", readToString(read));
        }

        boolean alreadyUnmapped = read.getReadUnmappedFlag();
        mConfig.UnmapRegions.checkTransformRead(read, mUnmapRegionState);

        boolean internallyUnmapped = !alreadyUnmapped && read.getReadUnmappedFlag();
        boolean fullyUnmapped = fullyUnmapped(read);

        // relocate any read which satisfies unmapping criteria - either to its mate's location or the fully-unmapped BAM
        if(internallyUnmapped || fullyUnmapped)
        {
            if(read.getSupplementaryAlignmentFlag() || read.isSecondaryAlignment())
                return;

            if(fullyUnmapped)
                mFullyUnmappedBamWriter.writeRecordSync(read);
            else
                mUnmappingBamWriter.writeRecordSync(read);
        }
    }

    private boolean processFullyUnmappedReads()
    {
        long totalUnmappedReads = 0;

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
        }

        return true;
    }
}
