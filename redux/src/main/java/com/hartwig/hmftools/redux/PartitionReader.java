package com.hartwig.hmftools.redux;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.sequencing.SBXBamUtils.stripDuplexIndels;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAP_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.secondsSinceNow;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.Constants.SUPP_ALIGNMENT_SCORE_MIN;
import static com.hartwig.hmftools.redux.common.FilterReadsType.NONE;
import static com.hartwig.hmftools.redux.common.FilterReadsType.readOutsideSpecifiedRegions;
import static com.hartwig.hmftools.redux.common.ReadInfo.readToString;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.HighDepthRegion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.redux.common.DuplicateGroup;
import com.hartwig.hmftools.redux.common.DuplicateGroupBuilder;
import com.hartwig.hmftools.redux.common.FragmentCoordReads;
import com.hartwig.hmftools.redux.common.FragmentStatus;
import com.hartwig.hmftools.redux.common.Statistics;
import com.hartwig.hmftools.redux.unmap.UnmapRegionState;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;
import com.hartwig.hmftools.redux.write.BamWriter;
import com.hartwig.hmftools.redux.write.PartitionInfo;

import htsjdk.samtools.SAMRecord;

public class PartitionReader
{
    private final ReduxConfig mConfig;

    private final BamReader mBamReader;
    private final ReadCache mReadCache;
    private final DuplicateGroupBuilder mDuplicateGroupBuilder;
    private final ConsensusReads mConsensusReads;

    // state for current partition
    private BamWriter mBamWriter;
    private List<ChrBaseRegion> mSliceRegions;
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

    // for SBX preprocessing
    private final static int SBX_REF_BASE_BUFFER_START = 2_000;
    private final static int SBX_REF_BASE_BUFFER_END = 10_000;
    private byte[] mRefBases;
    private int mRefBasesStart;

    public PartitionReader(final ReduxConfig config, final BamReader bamReader)
    {
        mConfig = config;
        mBamReader = bamReader;

        mReadCache = new ReadCache(ReadCache.DEFAULT_GROUP_SIZE, ReadCache.DEFAULT_MAX_SOFT_CLIP, mConfig.UMIs.Enabled, mConfig.Sequencing);

        mDuplicateGroupBuilder = new DuplicateGroupBuilder(config);
        mStats = mDuplicateGroupBuilder.statistics();
        mConsensusReads = new ConsensusReads(config.RefGenome, config.Sequencing, mStats.ConsensusStats);
        mConsensusReads.setDebugOptions(config.RunChecks);

        mCurrentRegion = null;
        mUnmapRegionState = null;
        mLastWriteLowerPosition = 0;

        mProcessedReads = 0;
        mNextLogReadCount = LOG_READ_COUNT;
        mLogReadIds = !mConfig.LogReadIds.isEmpty();
        mPcTotal = new PerformanceCounter("Total");
        mPcProcessDuplicates = new PerformanceCounter("ProcessDuplicates");

        mRefBases = null;
        mRefBasesStart = 0;
    }

    public List<PerformanceCounter> perfCounters()
    {
        return List.of(mPcTotal, mPcProcessDuplicates);
    }

    public Statistics statistics() { return mStats; }

    public void processPartition(final PartitionInfo partition)
    {
        mBamWriter = partition.bamWriter();
        mSliceRegions = partition.regions();

        int startReadCount = mProcessedReads;

        for(ChrBaseRegion region : mSliceRegions)
        {
            setupRegion(region);
            processRegion();
        }

        int processedReads = mProcessedReads - startReadCount;
        RD_LOGGER.debug("partition({}) reader complete, total reads({})", partition, processedReads);
    }

    @VisibleForTesting
    public void setupRegion(final ChrBaseRegion region)
    {
        mCurrentRegion = region;
        int chromosomeLength = mConfig.RefGenome.getChromosomeLength(region.Chromosome);
        mConsensusReads.setChromosomeLength(chromosomeLength);
        mLastWriteLowerPosition = 0;

        perfCountersStart();

        setUnmappedRegions();

        mBamWriter.initialiseRegion(region.Chromosome, region.start());

        if(mConfig.Sequencing == SBX)
        {
            mRefBases = null;
            checkRefBases(region.start());
        }
    }

    private void checkRefBases(int refPositionStart)
    {
        if(mConfig.Sequencing != SBX)
            return;

        if(mRefBases != null && refPositionStart < mRefBasesStart + SBX_REF_BASE_BUFFER_END)
            return;

        int chromosomeLength = mConfig.RefGenome.getChromosomeLength(mCurrentRegion.Chromosome);
        mRefBasesStart = max(refPositionStart - SBX_REF_BASE_BUFFER_START, 1);
        int refBaseEnd = min(mRefBasesStart + SBX_REF_BASE_BUFFER_END, chromosomeLength);
        mRefBases = mConfig.RefGenome.getBases(mCurrentRegion.Chromosome, mRefBasesStart, refBaseEnd);
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

    public static boolean fullyUnmapped(final SAMRecord read)
    {
        return read.getReadUnmappedFlag() && (!read.getReadPairedFlag() || read.getMateUnmappedFlag());
    }

    public static boolean shouldFilterRead(final SAMRecord read)
    {
        if(read.hasAttribute(CONSENSUS_READ_ATTRIBUTE)) // drop any consensus reads from previous MarkDup-generated BAMs runs
            return true;

        if(read.getSupplementaryAlignmentFlag())
        {
            Integer alignmentScore = read.getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE);
            if(alignmentScore != null && alignmentScore < SUPP_ALIGNMENT_SCORE_MIN)
                return true;
        }

        return false;
    }

    private static final int LOG_READ_COUNT = 1000000;

    private void preprocessSamRecord(final SAMRecord read)
    {
        if(mConfig.Sequencing == SBX)
        {
            stripDuplexIndels(read, mRefBases, mRefBasesStart);
        }
    }

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

        if(shouldFilterRead(read))
            return;

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
            // any read in an unmapping region has already been tested by the RegionUnmapper - scenarios:
            // 1. If the read was unmapped, it will have been relocated to its mate's coordinates and marked as internally unmapped
            // 2. That same read will also be encountered again at its original location - it will again satisfy the criteria to be unmapped
            // and then should be discarded from any further processing since it now a duplicate instance of the read in scenario 1
            // 3. A read with its mate already unmapped (ie prior to running Redux) - with the read also now unmapped - both of these should
            // be dropped without further processing

            // Further more, any read which was unmapped by the RegionUnmapper and has now appeared here can skip being checked
            boolean isUnmapped = read.getReadUnmappedFlag();
            boolean internallyUnmapped = isUnmapped && read.hasAttribute(UNMAP_ATTRIBUTE);

            if(!internallyUnmapped)
            {
                mConfig.UnmapRegions.checkTransformRead(read, mUnmapRegionState);

                if(!isUnmapped && read.getReadUnmappedFlag()) // scenario 2 as described above
                    return;

                if(fullyUnmapped(read)) // scenario 3 as described above
                    return;
            }
        }

        checkRefBases(read.getAlignmentStart());

        preprocessSamRecord(read);

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
            if(mConfig.FormConsensus)
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
    public void setBamWriter(final BamWriter writer) { mBamWriter = writer; }

    @VisibleForTesting
    public void processRead(final SAMRecord read) { processSamRecord(read); }

    @VisibleForTesting
    public void flushReadPositions()
    {
        processReadGroups(mReadCache.evictAll());
    }
}
