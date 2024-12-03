package com.hartwig.hmftools.redux;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.secondsSinceNow;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.Constants.SUPP_ALIGNMENT_SCORE_MIN;
import static com.hartwig.hmftools.redux.common.FilterReadsType.NONE;
import static com.hartwig.hmftools.redux.common.FilterReadsType.readOutsideSpecifiedRegions;
import static com.hartwig.hmftools.redux.common.ReadInfo.readToString;

import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.HighDepthRegion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.redux.common.DuplicateGroup;
import com.hartwig.hmftools.redux.common.DuplicateGroupBuilder;
import com.hartwig.hmftools.redux.common.FragmentCoordReads;
import com.hartwig.hmftools.redux.common.ReadInfo;
import com.hartwig.hmftools.redux.common.FragmentStatus;
import com.hartwig.hmftools.redux.common.Statistics;
import com.hartwig.hmftools.redux.common.UnmapRegionState;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;
import com.hartwig.hmftools.redux.write.BamWriter;

import htsjdk.samtools.SAMRecord;

public class PartitionReader implements Callable
{
    private final ReduxConfig mConfig;

    private final BamReader mBamReader;
    private final BamWriter mBamWriter;
    private final BamWriter mUnsortedBamWriter;
    private final ReadCache mReadCache;
    private final DuplicateGroupBuilder mDuplicateGroupBuilder;
    private final ConsensusReads mConsensusReads;

    private final List<ChrBaseRegion> mSliceRegions;
    private ChrBaseRegion mCurrentRegion;
    private int mLastWriteLowerPosition;
    private UnmapRegionState mUnmapRegionState;
    private final Statistics mStats;

    // debug
    private final boolean mLogReadIds;
    private final PerformanceCounter mPcTotal;
    private final PerformanceCounter mPcProcessDuplicates;
    private int mNextLogReadCount;
    private int mProcessedReads;

    public PartitionReader(
            final ReduxConfig config, final List<ChrBaseRegion> regions, final BamWriter bamWriter, final BamWriter unsortedBamWriter)
    {
        mConfig = config;
        mBamWriter = bamWriter;
        mUnsortedBamWriter = unsortedBamWriter;
        mBamReader = new BamReader(config);
        mSliceRegions = regions;

        mReadCache = new ReadCache(ReadCache.DEFAULT_GROUP_SIZE, ReadCache.DEFAULT_MAX_SOFT_CLIP, mConfig.UMIs.Enabled);

        mDuplicateGroupBuilder = new DuplicateGroupBuilder(config);
        mStats = mDuplicateGroupBuilder.statistics();
        mConsensusReads = new ConsensusReads(config.RefGenome, mStats.ConsensusStats);
        mConsensusReads.setDebugOptions(config.RunChecks);

        mCurrentRegion = null;
        mUnmapRegionState = null;
        mLastWriteLowerPosition = 0;

        mProcessedReads = 0;
        mNextLogReadCount = LOG_READ_COUNT;
        mLogReadIds = !mConfig.LogReadIds.isEmpty();
        mPcTotal = new PerformanceCounter("Total");
        mPcProcessDuplicates = new PerformanceCounter("ProcessDuplicates");
    }

    public List<PerformanceCounter> perfCounters()
    {
        return List.of(mPcTotal, mPcProcessDuplicates);
    }

    public Statistics statistics() { return mStats; }

    @Override
    public Long call()
    {
        // RD_LOGGER.debug("processing {} regions", mSliceRegions.size());

        for(ChrBaseRegion region : mSliceRegions)
        {
            setupRegion(region);
            processRegion();
        }

        return (long)0;
    }

    @VisibleForTesting
    public void setupRegion(final ChrBaseRegion region)
    {
        mCurrentRegion = region;
        mConsensusReads.setChromosomeLength(mConfig.RefGenome.getChromosomeLength(region.Chromosome));
        mLastWriteLowerPosition = 0;

        perfCountersStart();

        setUnmappedRegions();

        mBamWriter.initialiseRegion(region.Chromosome, region.start());
    }

    public void processRegion()
    {
        if(mBamReader != null)
        {
            mBamReader.sliceRegion(mCurrentRegion, this::processSamRecord);
        }

        postProcessRegion();
    }

    @VisibleForTesting
    public void postProcessRegion()
    {
        // post-slice clean-up
        processReadGroups(mReadCache.evictAll());

        mBamWriter.onRegionComplete();

        RD_LOGGER.debug("region({}) complete, processed {} reads", mCurrentRegion, mProcessedReads);

        perfCountersStop();
    }

    private static final int LOG_READ_COUNT = 1000000;

    private void processSamRecord(final SAMRecord read)
    {
        if(!mCurrentRegion.containsPosition(read.getAlignmentStart())) // to avoid processing reads from the prior region again
            return;

        ++mProcessedReads;

        if(mProcessedReads >= mNextLogReadCount)
        {
            RD_LOGGER.debug("region({}) position({}) processed {} reads, cache(coords={} reads={})",
                    mCurrentRegion, read.getAlignmentStart(), mProcessedReads,
                    mReadCache.cachedFragCoordGroups(), mReadCache.cachedReadCount());

            mNextLogReadCount += LOG_READ_COUNT;
        }

        if(mConfig.JitterMsiOnly)
        {
            mBamWriter.processJitterRead(read);
            return;
        }

        if(read.hasAttribute(CONSENSUS_READ_ATTRIBUTE)) // drop any consensus reads from previous MarkDup-generated BAMs runs
            return;

        if(read.getSupplementaryAlignmentFlag())
        {
            Integer alignmentScore = read.getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE);
            if(alignmentScore != null && alignmentScore < SUPP_ALIGNMENT_SCORE_MIN)
                return;
        }

        read.setDuplicateReadFlag(false);

        if(mConfig.SpecificRegionsFilterType != NONE && readOutsideSpecifiedRegions(
                read, mConfig.SpecificChrRegions.Regions, mConfig.SpecificChrRegions.Chromosomes, mConfig.SpecificRegionsFilterType))
        {
            return;
        }

        ++mStats.TotalReads;

        if(mLogReadIds && mConfig.LogReadIds.contains(read.getReadName())) // debugging only
        {
            RD_LOGGER.debug("specific read: {}", readToString(read));
        }

        if(mConfig.UnmapRegions.enabled())
        {
            boolean readUnmapped = mConfig.UnmapRegions.checkTransformRead(read, mUnmapRegionState);

            if(read.getReadUnmappedFlag())
            {
                if(read.getSupplementaryAlignmentFlag() || read.isSecondaryAlignment())
                    return; // drop unmapped supplementaries and secondaries

                if(SamRecordUtils.mateUnmapped(read))
                {
                    mBamWriter.writeRead(new ReadInfo(read, null));
                    ++mStats.Unmapped;
                    return;
                }
            }

            if(readUnmapped)
            {
                mUnsortedBamWriter.writeRead(read, FragmentStatus.UNSET);
                return;
            }
        }

        if(read.isSecondaryAlignment())
        {
            mBamWriter.setBoundaryPosition(read.getAlignmentStart(), false);
            mBamWriter.writeRead(read, FragmentStatus.UNSET);
            return;
        }

        try
        {
            mBamWriter.setBoundaryPosition(read.getAlignmentStart(), false);

            mReadCache.processRead(read);

            processReadGroups(mReadCache.popReads());
        }
        catch(Exception e)
        {
            RD_LOGGER.error("read({}) exception: {}", readToString(read), e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private static final double LOG_PERF_TIME_SEC = 1;
    private static final int LOG_PERF_FRAG_COUNT = 3000;

    private void processReadGroups(final FragmentCoordReads fragmentCoordReads)
    {
        if(fragmentCoordReads == null)
            return;

        mPcProcessDuplicates.resume();

        int readCount = fragmentCoordReads.totalReadCount();

        int minCachedReadPosition = mReadCache.minCachedReadStart();
        int currentPosition = mReadCache.currentReadMinPosition();
        int minPoppedReadsPosition = fragmentCoordReads.minReadPositionStart();
        int readFlushPosition = minCachedReadPosition > 0 ? minCachedReadPosition - 1 : minPoppedReadsPosition - 1;

        if(mLastWriteLowerPosition > 0 && minPoppedReadsPosition < mLastWriteLowerPosition) // mConfig.RunChecks
        {
            RD_LOGGER.warn("position({}:{})) processing earlier popped read min start({}) vs last({}) minCachedPos({})",
                    mCurrentRegion.Chromosome, currentPosition, minPoppedReadsPosition, mLastWriteLowerPosition,
                    minCachedReadPosition);

            // find read and report details about why it may not have been popped previously
            // mReadCache.logNonPoppedReads(fragmentCoordReads, mLastWriteLowerPosition);
        }

        mLastWriteLowerPosition = readFlushPosition;

        if(readCount > LOG_PERF_FRAG_COUNT)
        {
            RD_LOGGER.debug("position({}:{}) processing {} frag-coord read groups with {} reads",
                    mCurrentRegion.Chromosome, currentPosition, fragmentCoordReads.coordinateCount(), readCount);
        }

        boolean logDetails = mConfig.perfDebug() && readCount > LOG_PERF_FRAG_COUNT;
        long startTimeMs = logDetails ? System.currentTimeMillis() : 0;

        List<DuplicateGroup> duplicateGroups = mDuplicateGroupBuilder.processDuplicateGroups(
                fragmentCoordReads.DuplicateGroups, fragmentCoordReads.SingleReads, true);

        if(logDetails)
        {
            double timeTakenSec = secondsSinceNow(startTimeMs);

            if(timeTakenSec > mConfig.PerfDebugTime)
            {
                RD_LOGGER.debug("position({}:{}-{}) singles({}) groups({} reads={}) processing time({})",
                        mCurrentRegion.Chromosome, minPoppedReadsPosition, currentPosition, fragmentCoordReads.SingleReads.size(),
                        fragmentCoordReads.DuplicateGroups.size(), fragmentCoordReads.duplicateGroupReadCount(),
                        format("%.1fs", timeTakenSec));
            }
        }

        // write single fragments and duplicate groups
        for(DuplicateGroup duplicateGroup : duplicateGroups)
        {
            duplicateGroup.formConsensusRead(mConsensusReads);
            mBamWriter.writeDuplicateGroup(duplicateGroup);
        }

        mBamWriter.writeReads(fragmentCoordReads.SingleReads, true);

        mStats.addNonDuplicateCounts(fragmentCoordReads.SingleReads.size());

        // mark the lowest read position for reads able to be sorted and written
        mBamWriter.setBoundaryPosition(readFlushPosition, true);

        mPcProcessDuplicates.pause();
    }

    private void setUnmappedRegions()
    {
        List<HighDepthRegion> chrUnmapRegions = mConfig.UnmapRegions.getRegions(mCurrentRegion.Chromosome);

        List<HighDepthRegion> partitionRegions;
        if(chrUnmapRegions != null)
        {
            partitionRegions = chrUnmapRegions.stream().filter(x -> x.overlaps(mCurrentRegion)).collect(Collectors.toList());
        }
        else
        {
            partitionRegions = Collections.emptyList();
        }

        mUnmapRegionState = new UnmapRegionState(mCurrentRegion, partitionRegions);
    }

    private void perfCountersStart()
    {
        if(mConfig.perfDebug())
            mPcTotal.start(format("%s", mCurrentRegion));
        else
            mPcTotal.start();

        mPcProcessDuplicates.startPaused();
    }

    private void perfCountersStop()
    {
        mPcTotal.stop();
        mPcProcessDuplicates.stop();
    }

    @VisibleForTesting
    public void processRead(final SAMRecord read) { processSamRecord(read); }

    @VisibleForTesting
    public void flushReadPositions()
    {
        processReadGroups(mReadCache.evictAll());
    }
}
