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
import com.hartwig.hmftools.redux.old.DuplicateGroupOld;
import com.hartwig.hmftools.redux.old.FragmentOld;
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
    private int mPartitionRecordCount;
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

        mReadCache = new ReadCache(ReadCache.DEFAULT_GROUP_SIZE);

        mDuplicateGroupBuilder = new DuplicateGroupBuilder(config);
        mStats = mDuplicateGroupBuilder.statistics();
        mConsensusReads = new ConsensusReads(config.RefGenome, mStats.ConsensusStats);
        mConsensusReads.setDebugOptions(config.RunChecks);

        mCurrentRegion = null;
        mUnmapRegionState = null;

        mPartitionRecordCount = 0;

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
        // mReadPositions.evictAll();

        mBamWriter.onRegionComplete();

        perfCountersStop();

        RD_LOGGER.debug("partition({}) complete, reads({})", mCurrentRegion, mPartitionRecordCount);

        mPartitionRecordCount = 0;
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
        ++mPartitionRecordCount;

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

            processReadGroups();
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

    private void processReadGroups()
    {
        mPcProcessDuplicates.resume();

        FragmentCoordReads fragmentCoordReads = mReadCache.popReadsByFragmentCoordinates();

        if(fragmentCoordReads == null)
            return;

        int readCount = fragmentCoordReads.totalReadCount();

        if(readCount > LOG_PERF_FRAG_COUNT)
        {
            RD_LOGGER.debug("processing {} frag-coord read groups with {} reads", fragmentCoordReads.coordinateCount(), readCount);
        }

        boolean logDetails = mConfig.PerfDebug && readCount > LOG_PERF_FRAG_COUNT;
        long startTimeMs = logDetails ? System.currentTimeMillis() : 0;

        // TODO: check use of 'boolean requireOrientationMatch' when forming duplicate groups

        List<DuplicateGroup> duplicateGroups = mDuplicateGroupBuilder.processDuplicateGroups(
                fragmentCoordReads.DuplicateGroups, fragmentCoordReads.SingleReads, true);

        if(logDetails)
        {
            double timeTakenSec = secondsSinceNow(startTimeMs);

            /* TODO
            RD_LOGGER.debug("position({}:{}) fragments({}) resolved({}) dupGroups({}) candidates({}) processing time({})",
                    mCurrentRegion.Chromosome, position, posFragmentCount, nonDuplicateFragments.size(),
                    duplicateGroups != null ? duplicateGroups.size() : 0,
                    candidateDuplicatesList.stream().mapToInt(x -> x.fragmentCount()).sum(),
                    format("%.1fs", timeTakenSec));
            */
        }


        /*
        int posFragmentCount = positionFragments.size();
        int position = positionFragments.get(0).initialPosition();

        findDuplicateFragments(positionFragments, nonDuplicateFragments, positionDuplicateGroups, candidateDuplicatesList, mConfig.UMIs.Enabled);

        List<Fragment> singleFragments = mConfig.UMIs.Enabled ?
                nonDuplicateFragments.stream().filter(x -> x.status() == FragmentStatus.NONE).collect(Collectors.toList()) : Collections.EMPTY_LIST;

        List<DuplicateGroup> duplicateGroups = mDuplicateGroupBuilder.processDuplicateGroups(
                positionDuplicateGroups, true, singleFragments);

        if(logDetails)
        {
            double timeTakenSec = secondsSinceNow(startTimeMs);

            RD_LOGGER.debug("position({}:{}) fragments({}) resolved({}) dupGroups({}) candidates({}) processing time({})",
                    mCurrentRegion.Chromosome, position, posFragmentCount, nonDuplicateFragments.size(),
                    duplicateGroups != null ? duplicateGroups.size() : 0,
                    candidateDuplicatesList.stream().mapToInt(x -> x.fragmentCount()).sum(),
                    format("%.1fs", timeTakenSec));
        }

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

        if(!nonDuplicateFragments.isEmpty())
        {
            mBamWriter.writeReads(nonDuplicateFragments, true);

            mStats.LocalComplete += nonDuplicateFragments.stream().mapToInt(x -> x.readCount()).sum();

            int nonDuplicateFragments = (int)nonDuplicateFragments.stream().filter(x -> x.status() == FragmentStatus.NONE).count();
            mStats.addNonDuplicateCounts(nonDuplicateFragments);
        }
        */

        // TODO: what position should be set here to move the boundary for the sorted BAM writer?
        //  is it even required since this is now called directly from process read?
        int position = 0;

        if(position > 0) // note the reversed fragment coordinate positions are skipped for updating sorted BAM writing routines
            mBamWriter.setBoundaryPosition(position, true);

        mPcProcessDuplicates.pause();
    }

    /*

    private void processDuplicateGroup(final DuplicateGroupOld duplicateGroup)
    {
        // form consensus reads for any complete read leg groups and write reads
        List<SAMRecord> completeReads = duplicateGroup.popCompletedReads(mConsensusReads, false);
        mBamWriter.writeDuplicateGroup(duplicateGroup, completeReads);
    }

    private void processIncompleteRead(final SAMRecord read, final String basePartition)
    {
        if(basePartition == null)
        {
            // mate or supp is on a non-human chromosome, meaning it won't be retrieved - so write this immediately
            mBamWriter.writeRead(read, FragmentStatus.UNSET);
            ++mStats.LocalComplete;
            return;
        }

        if(basePartition.equals(mCurrentStrPartition))
        {
            PartitionResults partitionResults = mCurrentPartitionData.processIncompleteFragment(read);

            if(partitionResults != null)
            {
                if(partitionResults.umiGroups() != null || partitionResults.resolvedFragments() != null)
                {
                    if(partitionResults.umiGroups() != null)
                        partitionResults.umiGroups().forEach(x -> processDuplicateGroup(x));

                    if(partitionResults.resolvedFragments() != null)
                        mBamWriter.writeFragments(partitionResults.resolvedFragments(), true);
                }
                else if(partitionResults.fragmentStatus() != null && partitionResults.fragmentStatus().isResolved())
                {
                    mBamWriter.writeRead(read, partitionResults.fragmentStatus());
                }

                if(partitionResults.supplementaries() != null)
                {
                    partitionResults.supplementaries().forEach(x -> mBamWriter.writeSupplementary(x));
                }

                ++mStats.LocalComplete;
            }
            else
            {
                ++mStats.Incomplete;
            }
        }
        else
        {
            ++mStats.InterPartition;

            // cache this read and send through as groups when the partition is complete
            List<SAMRecord> pendingFragments = mPendingIncompleteReads.get(basePartition);

            if(pendingFragments == null)
            {
                pendingFragments = Lists.newArrayList();
                mPendingIncompleteReads.put(basePartition, pendingFragments);
            }

            pendingFragments.add(read);
        }
    }

    private void processPendingIncompletes()
    {
        if(mPendingIncompleteReads.isEmpty())
            return;

        if(mPendingIncompleteReads.size() > 100)
        {
            RD_LOGGER.debug("partition({}) processing {} remote reads from {} remote partitions",
                    mCurrentRegion, mPendingIncompleteReads.values().stream().mapToInt(x -> x.size()).sum(), mPendingIncompleteReads.size());
        }

        mPcPendingIncompletes.resume();

        for(Map.Entry<String,List<SAMRecord>> entry : mPendingIncompleteReads.entrySet())
        {
            String basePartition = entry.getKey();
            List<SAMRecord> reads = entry.getValue();

            PartitionData partitionData = mPartitionDataStore.getOrCreatePartitionData(basePartition);

            PartitionResults partitionResults = partitionData.processIncompleteFragments(reads);

            if(partitionResults.umiGroups() != null)
                partitionResults.umiGroups().forEach(x -> processDuplicateGroup(x));

            if(partitionResults.resolvedFragments() != null)
                mBamWriter.writeFragments(partitionResults.resolvedFragments(), true);

            if(partitionResults.supplementaries() != null)
                partitionResults.supplementaries().forEach(x -> mBamWriter.writeSupplementary(x));

        }

        mPendingIncompleteReads.clear();

        mPcPendingIncompletes.pause();
    }

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
        mReadCache.processRemaining();
        // mReadPositions.evictAll();
    }

    @VisibleForTesting
    public ConsensusReads consensusReads() { return mConsensusReads; }
}
