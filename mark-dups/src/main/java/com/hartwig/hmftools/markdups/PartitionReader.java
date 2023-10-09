package com.hartwig.hmftools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.secondsSinceNow;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;
import static com.hartwig.hmftools.markdups.common.DuplicateGroupBuilder.findDuplicateFragments;
import static com.hartwig.hmftools.markdups.common.FilterReadsType.NONE;
import static com.hartwig.hmftools.markdups.common.FilterReadsType.readOutsideSpecifiedRegions;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.formChromosomePartition;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.readToString;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.markdups.common.CandidateDuplicates;
import com.hartwig.hmftools.markdups.common.DuplicateGroup;
import com.hartwig.hmftools.markdups.common.DuplicateGroupBuilder;
import com.hartwig.hmftools.markdups.common.Fragment;
import com.hartwig.hmftools.markdups.common.FragmentStatus;
import com.hartwig.hmftools.markdups.common.PartitionData;
import com.hartwig.hmftools.markdups.common.PartitionResults;
import com.hartwig.hmftools.markdups.common.Statistics;
import com.hartwig.hmftools.markdups.consensus.ConsensusReads;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class PartitionReader implements Consumer<List<Fragment>>
{
    private final MarkDupsConfig mConfig;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final PartitionDataStore mPartitionDataStore;
    private final BamWriter mBamWriter;
    private final ReadPositionsCache mReadPositions;
    private final DuplicateGroupBuilder mDuplicateGroupBuilder;
    private final ConsensusReads mConsensusReads;

    private ChrBaseRegion mCurrentRegion;

    private String mCurrentStrPartition;
    private PartitionData mCurrentPartitionData;
    private final List<BaseRegion> mUnmapRegions;

    private final Map<String,List<SAMRecord>> mPendingIncompleteReads;

    private final boolean mLogReadIds;
    private int mPartitionRecordCount;
    private final Statistics mStats;
    private final PerformanceCounter mPcTotal;
    private final PerformanceCounter mPcAcceptPositions;
    private final PerformanceCounter mPcPendingIncompletes;

    public PartitionReader(
            final MarkDupsConfig config, final SamReader samReader,
            final BamWriter bamWriter, final PartitionDataStore partitionDataStore)
    {
        mConfig = config;
        mCurrentRegion = null;
        mPartitionDataStore = partitionDataStore;
        mBamWriter = bamWriter;

        mSamReader = samReader;

        mBamSlicer = new BamSlicer(0, true, true, false);
        mBamSlicer.setKeepUnmapped();

        mReadPositions = new ReadPositionsCache(config.BufferSize, !config.NoMateCigar, this);
        mDuplicateGroupBuilder = new DuplicateGroupBuilder(config);
        mConsensusReads = new ConsensusReads(config.RefGenome);
        mConsensusReads.setDebugOptions(config.RunChecks);

        mUnmapRegions = Lists.newArrayList();

        mPendingIncompleteReads = Maps.newHashMap();

        mPartitionRecordCount = 0;

        mStats = mDuplicateGroupBuilder.statistics();

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
        mPcTotal = new PerformanceCounter("Total");
        mPcAcceptPositions = new PerformanceCounter("AcceptPositions");
        mPcPendingIncompletes = new PerformanceCounter("PendingIncompletes");
    }

    public List<PerformanceCounter> perfCounters()
    {
        return List.of(mPcTotal, mPcAcceptPositions, mPcPendingIncompletes);
    }

    public Statistics statistics() { return mStats; }

    public void setupRegion(final ChrBaseRegion region)
    {
        mCurrentRegion = region;

        perfCountersStart();

        setUnmappedRegions();

        mCurrentStrPartition = formChromosomePartition(region.Chromosome, mCurrentRegion.start(), mConfig.PartitionSize);
        mCurrentPartitionData = mPartitionDataStore.getOrCreatePartitionData(mCurrentStrPartition);

        mReadPositions.setCurrentChromosome(region.Chromosome);

        // setExcludedRegion(ExcludedRegions.getPolyGRegion(mConfig.RefGenVersion));
    }

    public void processRegion()
    {
        if(mSamReader != null)
        {
            mBamSlicer.slice(mSamReader, mCurrentRegion, this::processSamRecord);
        }

        postProcessRegion();
    }

    @VisibleForTesting
    public void postProcessRegion()
    {
        // post-slice clean-up
        mReadPositions.evictAll();

        processPendingIncompletes();

        perfCountersStop();

        MD_LOGGER.debug("partition({}) complete, reads({})", mCurrentRegion, mPartitionRecordCount);

        if(mConfig.PerfDebug)
            mCurrentPartitionData.logCacheCounts();

        mPartitionRecordCount = 0;
    }

    private void processSamRecord(final SAMRecord read)
    {
        if(!mCurrentRegion.containsPosition(read.getAlignmentStart())) // to avoid processing reads from the prior region again
            return;

        if(mConfig.SpecificRegionsFilterType != NONE && readOutsideSpecifiedRegions(
                read, mConfig.SpecificChrRegions.Regions, mConfig.SpecificChrRegions.Chromosomes, mConfig.SpecificRegionsFilterType))
        {
            return;
        }

        ++mStats.TotalReads;
        ++mPartitionRecordCount;

        if(mLogReadIds && mConfig.LogReadIds.contains(read.getReadName())) // debugging only
        {
            MD_LOGGER.debug("specific read: {}", readToString(read));
        }

        if(mConfig.UnmapRegions.enabled())
        {
            mConfig.UnmapRegions.checkTransformRead(read, mUnmapRegions);

            if(read.getSupplementaryAlignmentFlag() && read.getReadUnmappedFlag())
                return; // drop unmapped supplementaries

            if(read.getReadUnmappedFlag() && read.getMateUnmappedFlag())
            {
                mBamWriter.writeFragment(new Fragment(read));
                ++mStats.Unmapped;
                return;
            }
        }

        try
        {
            if(!mReadPositions.processRead(read))
            {
                String basePartition = Fragment.getBasePartition(read, mConfig.PartitionSize);

                if(basePartition == null)
                {
                    // mate or supp is on a non-human chromsome, meaning it won't be retrieved - so write this immediately
                    mBamWriter.writeRead(read, FragmentStatus.UNSET);
                    ++mStats.LocalComplete;
                    return;
                }

                processIncompleteRead(read, basePartition);
            }
        }
        catch(Exception e)
        {
            MD_LOGGER.error("read({}) exception: {}", readToString(read), e.toString());
            e.printStackTrace();
            System.exit(1);
        }

        if(!mConfig.NoMateCigar && read.getReadPairedFlag() && !read.getMateUnmappedFlag() && !read.getSupplementaryAlignmentFlag()
                && !read.hasAttribute(MATE_CIGAR_ATTRIBUTE))
        {
            ++mStats.MissingMateCigar;
        }
    }

    private void processIncompleteRead(final SAMRecord read, final String basePartition)
    {
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
            MD_LOGGER.debug("partition({}) processing {} pending incomplete fragments",
                    mCurrentRegion, mPendingIncompleteReads.values().stream().mapToInt(x -> x.size()).sum());
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
        }

        mPendingIncompleteReads.clear();

        mPcPendingIncompletes.pause();
    }

    private void processDuplicateGroup(final DuplicateGroup duplicateGroup)
    {
        // form consensus reads for any complete read leg groups and write reads
        List<SAMRecord> completeReads = duplicateGroup.popCompletedReads(mConsensusReads, false);
        mBamWriter.writeDuplicateGroup(duplicateGroup, completeReads);
    }

    private static final double LOG_PERF_TIME_SEC = 1;
    private static final int LOG_PERF_FRAG_COUNT = 3000;

    public void accept(final List<Fragment> positionFragments)
    {
        if(positionFragments.isEmpty())
            return;

        mPcAcceptPositions.resume();
        List<Fragment> resolvedFragments = Lists.newArrayList();
        List<CandidateDuplicates> candidateDuplicatesList = Lists.newArrayList();
        List<List<Fragment>> positionDuplicateGroups = Lists.newArrayList();

        int posFragmentCount = positionFragments.size();
        int position = positionFragments.get(0).initialPosition();
        boolean logDetails = mConfig.PerfDebug && posFragmentCount > LOG_PERF_FRAG_COUNT; // was 10000
        long startTimeMs = logDetails ? System.currentTimeMillis() : 0;

        boolean inExcludedRegion = false; // dropped logic since added region unmapping logic

        //boolean inExcludedRegion = mExcludedRegion != null && positionFragments.stream()
        //        .anyMatch(x -> x.reads().stream().anyMatch(y -> FragmentUtils.overlapsExcludedRegion(mExcludedRegion, y)));

        findDuplicateFragments(positionFragments, resolvedFragments, positionDuplicateGroups, candidateDuplicatesList, mConfig.UMIs.Enabled);

        List<Fragment> singleFragments = mConfig.UMIs.Enabled && !inExcludedRegion ?
                resolvedFragments.stream().filter(x -> x.status() == FragmentStatus.NONE).collect(Collectors.toList()) : Collections.EMPTY_LIST;

        List<DuplicateGroup> duplicateGroups = mDuplicateGroupBuilder.processDuplicateGroups(
                positionDuplicateGroups, true, singleFragments, inExcludedRegion);

        if(logDetails)
        {
            double timeTakenSec = secondsSinceNow(startTimeMs);

            MD_LOGGER.debug("position({}:{}) fragments({}) resolved({}) dupGroups({}) candidates({}) processing time({})",
                    mCurrentRegion.Chromosome, position, posFragmentCount, resolvedFragments.size(), duplicateGroups.size(),
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
                MD_LOGGER.debug("position({}:{}) fragments({}) partition processing time({})",
                        mCurrentRegion.Chromosome, position, posFragmentCount, format("%.1fs", timeTakenSec));
            }
        }

        if(duplicateGroups != null)
            duplicateGroups.forEach(x -> processDuplicateGroup(x));

        if(!resolvedFragments.isEmpty())
        {
            mBamWriter.writeFragments(resolvedFragments, true);

            mStats.LocalComplete += resolvedFragments.stream().mapToInt(x -> x.readCount()).sum();
        }

        mPcAcceptPositions.pause();
    }

    private void setUnmappedRegions()
    {
        mUnmapRegions.clear();

        List<BaseRegion> chrRegions = mConfig.UnmapRegions.getRegions(mCurrentRegion.Chromosome);

        if(chrRegions != null)
            chrRegions.stream().filter(x -> x.overlaps(mCurrentRegion)).forEach(x -> mUnmapRegions.add(x));
    }

    private void perfCountersStart()
    {
        if(mConfig.PerfDebug)
            mPcTotal.start(format("%s", mCurrentRegion));
        else
            mPcTotal.start();

        mPcAcceptPositions.startPaused();
        mPcPendingIncompletes.startPaused();
    }

    private void perfCountersStop()
    {
        mPcTotal.stop();
        mPcAcceptPositions.stop();
        mPcPendingIncompletes.stop();
    }

    @VisibleForTesting
    public void processRead(final SAMRecord read) { processSamRecord(read); }

    @VisibleForTesting
    public void flushReadPositions() { mReadPositions.evictAll(); }

    @VisibleForTesting
    public void flushPendingIncompletes() { processPendingIncompletes(); }

    @VisibleForTesting
    public PartitionDataStore partitionDataStore() { return mPartitionDataStore; }

    @VisibleForTesting
    public ConsensusReads consensusReads() { return mConsensusReads; }
}
