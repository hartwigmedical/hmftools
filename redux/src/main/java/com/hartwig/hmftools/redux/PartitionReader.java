package com.hartwig.hmftools.redux;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.secondsSinceNow;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.Constants.SUPP_ALIGNMENT_SCORE_MIN;
import static com.hartwig.hmftools.redux.common.FilterReadsType.NONE;
import static com.hartwig.hmftools.redux.common.FilterReadsType.readOutsideSpecifiedRegions;
import static com.hartwig.hmftools.redux.old.FragmentUtils.readToString;

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
    private final ReadCache mReadCache;
    private final DuplicateGroupBuilder mDuplicateGroupBuilder;
    private final ConsensusReads mConsensusReads;

    private final List<ChrBaseRegion> mSliceRegions;
    private ChrBaseRegion mCurrentRegion;

    private UnmapRegionState mUnmapRegionState;

    private final boolean mLogReadIds;
    private final Statistics mStats;
    private final PerformanceCounter mPcTotal;
    private final PerformanceCounter mPcProcessDuplicates;

    public PartitionReader(
            final ReduxConfig config, final List<ChrBaseRegion> regions, final BamWriter bamWriter)
    {
        mConfig = config;
        mBamWriter = bamWriter;
        mBamReader = new BamReader(config);;
        mSliceRegions = regions;

        mReadCache = new ReadCache(ReadCache.DEFAULT_GROUP_SIZE, mConfig.UMIs.Enabled);

        mDuplicateGroupBuilder = new DuplicateGroupBuilder(config);
        mStats = mDuplicateGroupBuilder.statistics();
        mConsensusReads = new ConsensusReads(config.RefGenome, mStats.ConsensusStats);
        mConsensusReads.setDebugOptions(config.RunChecks);

        mCurrentRegion = null;
        mUnmapRegionState = null;

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

        perfCountersStop();
    }

    private void processSamRecord(final SAMRecord read)
    {
        if(!mCurrentRegion.containsPosition(read.getAlignmentStart())) // to avoid processing reads from the prior region again
            return;

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
            mConfig.UnmapRegions.checkTransformRead(read, mUnmapRegionState);

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
        int minProcessedReadsPosition = fragmentCoordReads.minReadPositionStart();
        int currentPosition = mReadCache.currentReadMinPosition();

        if(readCount > LOG_PERF_FRAG_COUNT)
        {
            RD_LOGGER.debug("position({}:{}-{}) processing {} frag-coord read groups with {} reads",
                    mCurrentRegion.Chromosome, minProcessedReadsPosition, currentPosition, fragmentCoordReads.coordinateCount(), readCount);
        }

        boolean logDetails = mConfig.PerfDebug && readCount > LOG_PERF_FRAG_COUNT;
        long startTimeMs = logDetails ? System.currentTimeMillis() : 0;

        List<DuplicateGroup> duplicateGroups = mDuplicateGroupBuilder.processDuplicateGroups(
                fragmentCoordReads.DuplicateGroups, fragmentCoordReads.SingleReads, true);

        if(logDetails)
        {
            double timeTakenSec = secondsSinceNow(startTimeMs);

            RD_LOGGER.debug("position({}:{}-{}) singles({}) groups({} reads={}) processing time({})",
                    mCurrentRegion.Chromosome, minProcessedReadsPosition, currentPosition, fragmentCoordReads.SingleReads.size(),
                    fragmentCoordReads.DuplicateGroups.size(),  fragmentCoordReads.duplicateGroupReadCount(),
                    format("%.1fs", timeTakenSec));
        }

        // write single fragments and duplicate groups
        for(DuplicateGroup duplicateGroup : duplicateGroups)
        {
            duplicateGroup.formConsensusRead(mConsensusReads);
            mBamWriter.writeDuplicateGroup(duplicateGroup);
        }

        mBamWriter.writeReads(fragmentCoordReads.SingleReads, true);

        mStats.addNonDuplicateCounts(fragmentCoordReads.SingleReads.size());

        mBamWriter.setBoundaryPosition(minProcessedReadsPosition, true);

        mPcProcessDuplicates.pause();
    }

    /*
    public void accept(final List<Fragment> positionFragments)
    {
        if(positionFragments.isEmpty())
            return;

        mPcProcessDuplicates.resume();
        List<Fragment> resolvedFragments = Lists.newArrayList();
        List<CandidateDuplicates> candidateDuplicatesList = Lists.newArrayList();
        List<List<Fragment>> positionDuplicateGroups = Lists.newArrayList();

        int posFragmentCount = positionFragments.size();
        int position = positionFragments.get(0).initialPosition();
        boolean logDetails = mConfig.PerfDebug && posFragmentCount > LOG_PERF_FRAG_COUNT; // was 10000
        long startTimeMs = logDetails ? System.currentTimeMillis() : 0;

        findDuplicateFragments(positionFragments, resolvedFragments, positionDuplicateGroups, candidateDuplicatesList, mConfig.UMIs.Enabled);

        List<Fragment> singleFragments = mConfig.UMIs.Enabled ?
                resolvedFragments.stream().filter(x -> x.status() == FragmentStatus.NONE).collect(Collectors.toList()) : Collections.EMPTY_LIST;

        List<DuplicateGroup> duplicateGroups = mDuplicateGroupBuilder.processDuplicateGroups(
                positionDuplicateGroups, true, singleFragments);

        if(logDetails)
        {
            double timeTakenSec = secondsSinceNow(startTimeMs);

            RD_LOGGER.debug("position({}:{}) fragments({}) resolved({}) dupGroups({}) candidates({}) processing time({})",
                    mCurrentRegion.Chromosome, position, posFragmentCount, resolvedFragments.size(),
                    duplicateGroups != null ? duplicateGroups.size() : 0,
                    candidateDuplicatesList.stream().mapToInt(x -> x.fragmentCount()).sum(),
                    format("%.1fs", timeTakenSec));
        }

        startTimeMs = logDetails ? System.currentTimeMillis() : 0;

        mCurrentPartitionData.processPrimaryFragments(resolvedFragments, candidateDuplicatesList, duplicateGroups);

        if(logDetails)
        {
            double timeTakenSec = secondsSinceNow(startTimeMs);

            if(timeTakenSec >= LOG_PERF_TIME_SEC)
            {
                RD_LOGGER.debug("position({}:{}) fragments({}) partition processing time({})",
                        mCurrentRegion.Chromosome, position, posFragmentCount, format("%.1fs", timeTakenSec));
            }
        }

        if(duplicateGroups != null)
            duplicateGroups.forEach(x -> processDuplicateGroup(x));

        if(!resolvedFragments.isEmpty())
        {
            mBamWriter.writeFragments(resolvedFragments, true);

            mStats.LocalComplete += resolvedFragments.stream().mapToInt(x -> x.readCount()).sum();

            int nonDuplicateFragments = (int)resolvedFragments.stream().filter(x -> x.status() == FragmentStatus.NONE).count();
            mStats.addNonDuplicateCounts(nonDuplicateFragments);
        }

        if(position > 0) // note the reversed fragment coordinate positions are skipped for updating sorted BAM writing routines
            mBamWriter.setBoundaryPosition(position, true);

        mPcProcessDuplicates.pause();
    }
    */

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
        if(mConfig.PerfDebug)
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
