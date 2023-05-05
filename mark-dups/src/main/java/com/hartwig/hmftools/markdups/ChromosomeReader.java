package com.hartwig.hmftools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.markdups.common.DuplicateGroups.findDuplicateFragments;
import static com.hartwig.hmftools.markdups.common.FilterReadsType.readOutsideSpecifiedRegions;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.formChromosomePartition;
import static com.hartwig.hmftools.markdups.common.FragmentUtils.readToString;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.function.Consumer;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.markdups.common.CandidateDuplicates;
import com.hartwig.hmftools.markdups.common.DuplicateGroups;
import com.hartwig.hmftools.markdups.common.Fragment;
import com.hartwig.hmftools.markdups.common.FragmentStatus;
import com.hartwig.hmftools.markdups.common.PartitionData;
import com.hartwig.hmftools.markdups.common.PartitionResults;
import com.hartwig.hmftools.markdups.common.Statistics;
import com.hartwig.hmftools.markdups.umi.ConsensusReads;
import com.hartwig.hmftools.markdups.umi.UmiGroup;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class ChromosomeReader implements Consumer<List<Fragment>>, Callable
{
    private final MarkDupsConfig mConfig;
    private final ChrBaseRegion mRegion;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final PartitionDataStore mPartitionDataStore;
    private final RecordWriter mRecordWriter;
    private final ReadPositionsCache mReadPositions;
    private final DuplicateGroups mDuplicateGroups;
    private final ConsensusReads mConsensusReads;

    private BaseRegion mCurrentPartition;
    private String mCurrentStrPartition;
    private PartitionData mCurrentPartitionData;

    private final Map<String,List<SAMRecord>> mPendingIncompleteReads;

    private final boolean mLogReadIds;
    private int mPartitionRecordCount;
    private final Statistics mStats;
    private final PerformanceCounter mPerfCounter;

    public ChromosomeReader(
            final ChrBaseRegion region, final MarkDupsConfig config, final RecordWriter recordWriter,
            final PartitionDataStore partitionDataStore)
    {
        mConfig = config;
        mRegion = region;
        mPartitionDataStore = partitionDataStore;
        mRecordWriter = recordWriter;

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, true, true, true);
        mBamSlicer.setKeepUnmapped();

        mReadPositions = new ReadPositionsCache(region.Chromosome, config.BufferSize, !config.NoMateCigar, this);
        mDuplicateGroups = new DuplicateGroups(config.UMIs);
        mConsensusReads = new ConsensusReads(config.UMIs, config.RefGenome);

        if(!mConfig.SpecificRegions.isEmpty())
        {
            // NOTE: doesn't currently handle multiple regions on the same chromosome
            ChrBaseRegion firstRegion = mConfig.SpecificRegions.stream().filter(x -> x.Chromosome.equals(mRegion.Chromosome)).findFirst().orElse(mRegion);
            int partitionStart = (firstRegion.start() / mConfig.PartitionSize) * mConfig.PartitionSize;
            mCurrentPartition = new BaseRegion(partitionStart, partitionStart + mConfig.PartitionSize - 1);
        }
        else
        {
            mCurrentPartition = new BaseRegion(1, mConfig.PartitionSize - 1);
        }

        mCurrentStrPartition = formChromosomePartition(mRegion.Chromosome, mCurrentPartition.start(), mConfig.PartitionSize);
        mCurrentPartitionData = mPartitionDataStore.getOrCreatePartitionData(mCurrentStrPartition);

        mPendingIncompleteReads = Maps.newHashMap();

        mPartitionRecordCount = 0;

        mStats = mDuplicateGroups.statistics();

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
        mPerfCounter = new PerformanceCounter("Slice");
    }

    public PerformanceCounter perfCounter() { return mPerfCounter; }
    public Statistics statistics() { return mStats; }

    @Override
    public Long call()
    {
        run();
        return (long)1;
    }

    public void run()
    {
        perfCounterStart();

        if(!mConfig.SpecificRegions.isEmpty())
        {
            for(ChrBaseRegion region : mConfig.SpecificRegions)
            {
                if(!region.Chromosome.equals(mRegion.Chromosome))
                    continue;

                MD_LOGGER.debug("processing specific region({})", region);
                mBamSlicer.slice(mSamReader, Lists.newArrayList(region), this::processSamRecord);
            }
        }
        else
        {
            MD_LOGGER.info("processing chromosome({})", mRegion.Chromosome);
            mBamSlicer.slice(mSamReader, Lists.newArrayList(mRegion), this::processSamRecord);
        }

        onPartitionComplete(false);

        MD_LOGGER.info("chromosome({}) complete, reads({})", mRegion.Chromosome, mStats.TotalReads);

        mConsensusReads.logStats(mRegion.Chromosome);
    }

    private void onPartitionComplete(boolean setupNext)
    {
        mReadPositions.evictAll();

        processPendingIncompletes();

        mPerfCounter.stop();

        MD_LOGGER.debug("partition({}:{}) complete, reads({})", mRegion.Chromosome, mCurrentPartition, mPartitionRecordCount);

        if(mConfig.PerfDebug)
            mCurrentPartitionData.logCacheCounts();

        mPartitionRecordCount = 0;

        if(setupNext)
        {
            // move ahead to the next partition, until the end of the chromosome is reached
            int regionStart = mCurrentPartition.end() + 1;

            if(regionStart > mRegion.end())
            {
                mCurrentPartition = null;
                return;
            }

            mCurrentPartition.setStart(regionStart);
            mCurrentPartition.setEnd(regionStart + mConfig.PartitionSize - 1);

            mCurrentStrPartition = formChromosomePartition(mRegion.Chromosome, mCurrentPartition.start(), mConfig.PartitionSize);
            mCurrentPartitionData = mPartitionDataStore.getOrCreatePartitionData(mCurrentStrPartition);

            perfCounterStart();
        }

        System.gc();
    }

    private void processSamRecord(final SAMRecord read)
    {
        if(readOutsideSpecifiedRegions(read, mConfig.SpecificRegions, mConfig.SpecificChromosomes, mConfig.SpecificRegionsFilterType))
            return;

        ++mStats.TotalReads;
        ++mPartitionRecordCount;

        if(mConfig.RunChecks)
            mRecordWriter.registerRead(read);

        int readStart = read.getAlignmentStart();

        while(mCurrentPartition != null && readStart > mCurrentPartition.end())
        {
            onPartitionComplete(true);

            if(mCurrentPartition == null)
            {
                mRecordWriter.writeFragment(new Fragment(read));
                return;
            }
        }

        if(mLogReadIds && mConfig.LogReadIds.contains(read.getReadName())) // debugging only
        {
            MD_LOGGER.debug("specific read: {}", readToString(read));
        }

        if(!mConfig.NoMateCigar && read.getReadPairedFlag() && !read.getMateUnmappedFlag() && !read.getSupplementaryAlignmentFlag()
        && !read.hasAttribute(MATE_CIGAR_ATTRIBUTE))
        {
            ++mStats.MissingMateCigar;
        }

        try
        {
            if(!mReadPositions.processRead(read))
            {
                ++mStats.Incomplete;

                String basePartition = Fragment.getBasePartition(read, mConfig.PartitionSize);

                if(basePartition == null)
                {
                    // mate or supp is on a non-human chromsome, meaning it won't be retrieved - so write this immediately
                    mRecordWriter.writeRead(read, FragmentStatus.UNSET);
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
                        partitionResults.umiGroups().forEach(x -> processUmiGroup(x));

                    if(partitionResults.resolvedFragments() != null)
                        mRecordWriter.writeFragments(partitionResults.resolvedFragments(), true);
                }
                else if(partitionResults.fragmentStatus() != null && partitionResults.fragmentStatus().isResolved())
                {
                    mRecordWriter.writeRead(read, partitionResults.fragmentStatus());
                }
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
            MD_LOGGER.debug("partition({}:{}) processing {} pending incomplete fragments",
                    mRegion.Chromosome, mCurrentPartition, mPendingIncompleteReads.values().stream().mapToInt(x -> x.size()).sum());
        }

        for(Map.Entry<String,List<SAMRecord>> entry : mPendingIncompleteReads.entrySet())
        {
            String basePartition = entry.getKey();
            List<SAMRecord> reads = entry.getValue();

            PartitionData partitionData = mPartitionDataStore.getOrCreatePartitionData(basePartition);

            PartitionResults partitionResults = partitionData.processIncompleteFragments(reads);

            if(partitionResults.umiGroups() != null)
                partitionResults.umiGroups().forEach(x -> processUmiGroup(x));

            if(partitionResults.resolvedFragments() != null)
                mRecordWriter.writeFragments(partitionResults.resolvedFragments(), true);
        }

        mPendingIncompleteReads.clear();
    }

    private void processUmiGroup(final UmiGroup umiGroup)
    {
        // form consensus reads for any complete read leg groups and write reads
        List<SAMRecord> completeReads = umiGroup.popCompletedReads(mConsensusReads, false);
        mRecordWriter.writeUmiReads(umiGroup, completeReads);
    }

    public void accept(final List<Fragment> positionFragments)
    {
        if(positionFragments.isEmpty())
            return;

        List<Fragment> resolvedFragments = Lists.newArrayList();
        List<CandidateDuplicates> candidateDuplicatesList = Lists.newArrayList();
        List<List<Fragment>> duplicateGroups = Lists.newArrayList();

        int posFragmentCount = positionFragments.size();
        boolean logDetails = mConfig.PerfDebug && posFragmentCount > 10000;
        long startTimeMs = logDetails ? System.currentTimeMillis() : 0;

        int position = positionFragments.get(0).initialPosition();

        findDuplicateFragments(positionFragments, resolvedFragments, duplicateGroups, candidateDuplicatesList);

        List<UmiGroup> umiGroups = null;

        if(mConfig.UMIs.Enabled)
        {
            umiGroups = mDuplicateGroups.processDuplicateUmiGroups(duplicateGroups);
        }
        else
        {
            mDuplicateGroups.processDuplicateGroups(duplicateGroups);
        }

        if(logDetails)
        {
            double timeTakenSec = (System.currentTimeMillis() - startTimeMs) / 1000.0;

            if(timeTakenSec >= 1.0)
            {
                MD_LOGGER.debug("position({}:{}) fragments({}) resolved({}) candidates({}) processing time({})",
                        mRegion.Chromosome, position, posFragmentCount, resolvedFragments.size(),
                        candidateDuplicatesList.stream().mapToInt(x -> x.fragmentCount()).sum(),
                        format("%.1fs", timeTakenSec));
            }
        }

        startTimeMs = System.currentTimeMillis();

        mCurrentPartitionData.processPrimaryFragments(resolvedFragments, candidateDuplicatesList, umiGroups);

        double timeTakenSec = (System.currentTimeMillis() - startTimeMs) / 1000.0;

        if(timeTakenSec >= 1.0)
        {
            MD_LOGGER.debug("position({}:{}) fragments({}) partition processing time({})",
                    mRegion.Chromosome, position, posFragmentCount, format("%.1fs", timeTakenSec));
        }

        if(umiGroups != null)
            umiGroups.forEach(x -> processUmiGroup(x));

        if(!resolvedFragments.isEmpty())
        {
            mRecordWriter.writeFragments(resolvedFragments, true);

            mStats.LocalComplete += (int)resolvedFragments.stream().filter(x -> x.allReadsPresent()).count();
        }
    }

    private void perfCounterStart()
    {
        if(mConfig.PerfDebug)
            mPerfCounter.start(format("%s:%s", mRegion.Chromosome, mCurrentPartition));
        else
            mPerfCounter.start();
    }

    @VisibleForTesting
    public void processRead(final SAMRecord read) { processSamRecord(read); }

    @VisibleForTesting
    public void flushReadPositions() { mReadPositions.evictAll(); }

    @VisibleForTesting
    public void flushPendingIncompletes() { processPendingIncompletes(); }

    @VisibleForTesting
    public void onChromosomeComplete() { onPartitionComplete(false); }

    @VisibleForTesting
    public PartitionDataStore partitionDataStore() { return mPartitionDataStore; }
}
