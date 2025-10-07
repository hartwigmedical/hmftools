package com.hartwig.hmftools.redux;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNMAP_ATTRIBUTE;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.secondsSinceNow;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ReduxConfig.SEQUENCING_TYPE;
import static com.hartwig.hmftools.redux.ReduxConfig.isIllumina;
import static com.hartwig.hmftools.redux.ReduxConfig.isSbx;
import static com.hartwig.hmftools.redux.ReduxConfig.isUltima;
import static com.hartwig.hmftools.redux.ReduxConstants.SUPP_ALIGNMENT_SCORE_MIN;
import static com.hartwig.hmftools.redux.common.FilterReadsType.NONE;
import static com.hartwig.hmftools.redux.common.FilterReadsType.readOutsideSpecifiedRegions;
import static com.hartwig.hmftools.redux.common.ReadInfo.readToString;
import static com.hartwig.hmftools.redux.consensus.SbxRoutines.stripDuplexIndels;

import static org.apache.logging.log4j.Level.DEBUG;
import static org.apache.logging.log4j.Level.TRACE;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.perf.PerformanceCounter;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.mappability.UnmappingRegion;
import com.hartwig.hmftools.redux.common.DuplicateGroup;
import com.hartwig.hmftools.redux.common.DuplicateGroupBuilder;
import com.hartwig.hmftools.redux.common.FragmentCoordReads;
import com.hartwig.hmftools.redux.common.ReadInfo;
import com.hartwig.hmftools.redux.common.Statistics;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;
import com.hartwig.hmftools.redux.consensus.SbxRoutines;
import com.hartwig.hmftools.redux.consensus.UltimaRoutines;
import com.hartwig.hmftools.redux.unmap.ReadUnmapper;
import com.hartwig.hmftools.redux.unmap.UnmapRegionState;
import com.hartwig.hmftools.redux.write.BamWriter;
import com.hartwig.hmftools.redux.write.PartitionInfo;

import org.apache.logging.log4j.Level;

import htsjdk.samtools.SAMRecord;

public class PartitionReader
{
    private final ReduxConfig mConfig;

    private final BamReader mBamReader;
    private final IReadCache mReadCache;
    private final ReadUnmapper mReadUnmapper;
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

    public PartitionReader(final ReduxConfig config, final BamReader bamReader)
    {
        mConfig = config;
        mBamReader = bamReader;
        mReadUnmapper = mConfig.UnmapRegions;

        if(isIllumina())
        {
            if(mConfig.UMIs.Enabled)
            {
                mReadCache = new JitterReadCache(new ReadCache(
                        ReadCache.DEFAULT_GROUP_SIZE, ReadCache.DEFAULT_MAX_SOFT_CLIP, mConfig.UMIs.Enabled, mConfig.DuplicateConfig));
            }
            else
            {
                mReadCache = new ReadCache(
                        ReadCache.DEFAULT_GROUP_SIZE, ReadCache.DEFAULT_MAX_SOFT_CLIP, mConfig.UMIs.Enabled, mConfig.DuplicateConfig);
            }
        }
        else
        {
            // the sampled max read length is doubled, because it has been observed in non-Illumina bams that the max read length is usually
            // larger than the sampled max read length extra room is required
            mReadCache = new ReadCache(
                    3 * mConfig.readLength(),
                    2 * mConfig.readLength() - 1, mConfig.UMIs.Enabled, mConfig.DuplicateConfig);
        }

        mDuplicateGroupBuilder = new DuplicateGroupBuilder(config);
        mStats = mDuplicateGroupBuilder.statistics();
        mConsensusReads = new ConsensusReads(config.RefGenome, SEQUENCING_TYPE, mStats.ConsensusStats);

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

        mBamWriter.initialiseRegion(region);
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

        mBamWriter.regionComplete();

        Level logLevel = isAltContigRegion() ? TRACE : DEBUG;
        RD_LOGGER.log(logLevel, "region({}) complete, processed {} reads", mCurrentRegion, mProcessedReads);

        perfCountersStop();
    }

    public static boolean fullyUnmapped(final SAMRecord read)
    {
        return read.getReadUnmappedFlag() && (!read.getReadPairedFlag() || read.getMateUnmappedFlag());
    }

    public static boolean shouldFilterRead(final SAMRecord read)
    {
        if(read.hasAttribute(CONSENSUS_READ_ATTRIBUTE)) // drop any consensus reads from previous Redux-generated BAMs runs
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
        if(isSbx())
        {
            try
            {
                stripDuplexIndels(read);
            }
            catch(Exception e)
            {
                RD_LOGGER.error("preprocess read({}) error: {}", readToString(read), e.toString());
                e.printStackTrace();
            }
        }
    }

    private void postProcessSingleReads(final List<ReadInfo> singleReads)
    {
        singleReads.forEach(x -> finaliseRead(x.read()));
    }

    private void postProcessPrimaryRead(final DuplicateGroup duplicateGroup)
    {
        if(duplicateGroup.consensusRead() != null)
            finaliseRead(duplicateGroup.consensusRead());
        else if(duplicateGroup.primaryRead() != null)
            finaliseRead(duplicateGroup.primaryRead());
    }

    private void finaliseRead(final SAMRecord read)
    {
        // TODO: make an interface for seq-tech routines
        if(isSbx())
        {
            try
            {
                SbxRoutines.finaliseRead(mConfig.RefGenome, read);
            }
            catch(Exception e)
            {
                RD_LOGGER.error("finaluse read error: {}", readToString(read));
            }
        }
        else if(isUltima())
        {
            UltimaRoutines.finaliseRead(mConfig.RefGenome, read);
        }
    }

    private void processSamRecord(final SAMRecord read)
    {
        int readStart = read.getAlignmentStart();

        if(!mCurrentRegion.containsPosition(readStart)) // to avoid processing reads from the prior region again
            return;

        ++mProcessedReads;

        if(mProcessedReads >= mNextLogReadCount)
        {
            double processedReads = mProcessedReads / 1000000.0;
            RD_LOGGER.debug("region({}) position({}) processed {}M reads, cache(coords={} reads={})",
                    mCurrentRegion, readStart, format("%.0f", processedReads),
                    mReadCache.cachedFragCoordGroups(), mReadCache.cachedReadCount());

            mNextLogReadCount += LOG_READ_COUNT;
        }

        if(mLogReadIds && mConfig.LogReadIds.contains(read.getReadName())) // debugging only
        {
            RD_LOGGER.debug("specific read: {}", readToString(read));
        }

        if(mConfig.BqrAndJitterMsiOnly)
        {
            mBamWriter.captureReadInfo(read);
            return;
        }

        if(shouldFilterRead(read))
            return;

        if(!read.isSecondaryAlignment() && read.getReadPairedFlag() && !read.getMateUnmappedFlag() && !read.hasAttribute(MATE_CIGAR_ATTRIBUTE))
        {
            if(!read.getSupplementaryAlignmentFlag() || mConfig.FailOnMissingSuppMateCigar)
            {
                RD_LOGGER.error("read({}) missing mate CIGAR", readToString(read));
                System.exit(1);
            }
        }

        ++mStats.TotalReads;

        read.setDuplicateReadFlag(false);

        if(mConfig.SpecificRegionsFilterType != NONE && readOutsideSpecifiedRegions(
                read, mConfig.SpecificChrRegions.Regions, mConfig.SpecificChrRegions.Chromosomes, mConfig.SpecificRegionsFilterType))
        {
            return;
        }

        if(mReadUnmapper.enabled())
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
                mReadUnmapper.checkTransformRead(read, mUnmapRegionState);

                boolean fullyUnmapped = fullyUnmapped(read);
                boolean unmapped = !isUnmapped && read.getReadUnmappedFlag();

                if(unmapped || fullyUnmapped)
                {
                    mConfig.readChecker().checkRead(read, mCurrentRegion.Chromosome, readStart, fullyUnmapped);
                    return;
                }
            }
        }

        preprocessSamRecord(read);

        if(read.isSecondaryAlignment() || mConfig.SkipDuplicateMarking)
        {
            if(mConfig.SkipDuplicateMarking)
                finaliseRead(read);

            mBamWriter.setBoundaryPosition(readStart, false);
            mBamWriter.writeNonDuplicateRead(read);
            return;
        }

        try
        {
            mBamWriter.setBoundaryPosition(readStart, false);

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
            // do not form consensus if duplicateGroup only contains one non poly-g umi read
            if(mConfig.FormConsensus && duplicateGroup.readCount() - duplicateGroup.polyGUmiReads().size() >= 2)
            {
                duplicateGroup.formConsensusRead(mConsensusReads);
                mBamWriter.setBoundaryPosition(duplicateGroup.consensusRead().getAlignmentStart(), false);
            }

            postProcessPrimaryRead(duplicateGroup);
            mBamWriter.writeDuplicateGroup(duplicateGroup);

            if(mConfig.LogDuplicateGroupSize > 0 && duplicateGroup.readCount() >= mConfig.LogDuplicateGroupSize)
            {
                RD_LOGGER.debug("dup group size({}) coords({})", duplicateGroup.readCount(), duplicateGroup.fragmentCoordinates());
            }
        }

        postProcessSingleReads(fragmentCoordReads.SingleReads);
        mBamWriter.writeNonDuplicateReads(fragmentCoordReads.SingleReads);

        mStats.addNonDuplicateCounts(fragmentCoordReads.SingleReads.size());

        // mark the lowest read position for reads able to be sorted and written
        mBamWriter.setBoundaryPosition(readFlushPosition, true);

        mPcProcessDuplicates.pause();
    }

    private void setUnmappedRegions()
    {
        List<UnmappingRegion> chrUnmapRegions = mReadUnmapper.getRegions(mCurrentRegion.Chromosome);

        List<UnmappingRegion> partitionRegions;
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
        if(isAltContigRegion())
            return;

        if(mConfig.perfDebug())
            mPcTotal.start(format("%s", mCurrentRegion));
        else
            mPcTotal.start();

        mPcProcessDuplicates.startPaused();
    }

    private void perfCountersStop()
    {
        if(isAltContigRegion())
            return;

        mPcTotal.stop();
        mPcProcessDuplicates.stop();
    }

    private boolean isAltContigRegion() { return Chromosome.isAltRegionContig(mCurrentRegion.Chromosome); }

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
