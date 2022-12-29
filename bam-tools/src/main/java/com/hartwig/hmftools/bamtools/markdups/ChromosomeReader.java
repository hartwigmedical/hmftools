package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.checkFragmentClassification;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.classifyFragments;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.formChromosomePartition;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.readInSpecifiedRegions;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.readToString;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class ChromosomeReader implements Consumer<CandidateDuplicates>, Callable
{
    private final MarkDupsConfig mConfig;
    private final ChrBaseRegion mRegion;
    private final BaseRegion mCurrentPartition;
    private String mCurrentStrPartition;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final GroupCombiner mRemoteGroupCombiner;
    private final GroupCombiner mLocalGroupCombiner;
    private final RecordWriter mRecordWriter;
    private final ReadPositionsCache mReadPositions;

    private final List<Fragment> mPartitionSupplementaries;
    private final Map<String,Fragment> mPartitionResolvedFragments;

    private final boolean mLogReadIds;
    private int mTotalRecordCount;
    private int mPartitionRecordCount;
    private int mMaxPositionFragments;
    private final Set<SAMRecord> mReadsProcessed;
    private final DuplicateStats mStats;
    private final PerformanceCounter mPerfCounter;

    public ChromosomeReader(
            final ChrBaseRegion region, final MarkDupsConfig config, final RecordWriter recordWriter,
            final GroupCombiner groupCombiner)
    {
        mConfig = config;
        mRegion = region;
        mRemoteGroupCombiner = groupCombiner;
        mRecordWriter = recordWriter;
        mLocalGroupCombiner = new GroupCombiner(mRecordWriter, true, false);

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, true, true, true);
        mBamSlicer.setKeepUnmapped();

        mReadPositions = new ReadPositionsCache(region.Chromosome, config.BufferSize, this);

        mPartitionSupplementaries = Lists.newArrayList();
        mPartitionResolvedFragments = Maps.newHashMap();

        if(!mConfig.SpecificRegions.isEmpty())
        {
            ChrBaseRegion firstRegion = mConfig.SpecificRegions.stream().filter(x -> x.Chromosome.equals(mRegion.Chromosome)).findFirst().orElse(mRegion);
            int partitionStart = (firstRegion.start() / mConfig.PartitionSize) * mConfig.PartitionSize;
            mCurrentPartition = new BaseRegion(partitionStart, partitionStart + mConfig.PartitionSize - 1);
        }
        else
        {
            mCurrentPartition = new BaseRegion(1, mConfig.PartitionSize);
        }

        mCurrentStrPartition = formChromosomePartition(mRegion.Chromosome, mCurrentPartition.start(), mConfig.PartitionSize);
        mTotalRecordCount = 0;
        mPartitionRecordCount = 0;
        mMaxPositionFragments = 0;

        mStats = new DuplicateStats();

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
        mReadsProcessed = Sets.newHashSet();
        mPerfCounter = new PerformanceCounter("Slice");
    }

    public int totalRecordCount() { return mTotalRecordCount; }
    public Set<SAMRecord> readsProcessed() { return mReadsProcessed; }
    public PerformanceCounter perfCounter() { return mPerfCounter; }
    public DuplicateStats duplicateStats() { return mStats; }

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

                BM_LOGGER.debug("processing specific region({})", region);
                mBamSlicer.slice(mSamReader, Lists.newArrayList(region), this::processSamRecord);
            }
        }
        else
        {
            BM_LOGGER.info("processing chromosome({})", mRegion.Chromosome);
            mBamSlicer.slice(mSamReader, Lists.newArrayList(mRegion), this::processSamRecord);
        }

        onPartitionComplete(false);

        BM_LOGGER.info("chromosome({}) complete, reads({})", mRegion.Chromosome, mTotalRecordCount);
    }

    private void onPartitionComplete(boolean setupNext)
    {
        String readPositionsStats = BM_LOGGER.isDebugEnabled() ? mReadPositions.cacheStatsStr() : "";

        mReadPositions.evictAll();

        List<Fragment> resolvedFragments = mPartitionResolvedFragments.values().stream().collect(Collectors.toList());

        List<CandidateDuplicates> candidateDuplicates = Lists.newArrayList();

        mLocalGroupCombiner.gatherPartitionFragments(mCurrentStrPartition, mCurrentPartition, resolvedFragments, candidateDuplicates);

        mStats.InterPartitionUnclear += candidateDuplicates.stream().mapToInt(x -> x.Fragments.size()).sum();
        mStats.ReadCount += mPartitionRecordCount;

        mRemoteGroupCombiner.processPartitionFragments(mCurrentStrPartition, resolvedFragments, candidateDuplicates, mPartitionSupplementaries);

        mPerfCounter.stop();

        BM_LOGGER.debug("partition({}:{}) complete, reads({}) remotes cached(supps={} resolved={} candidates={}) maxPosFrags({}) {}",
                mRegion.Chromosome, mCurrentPartition, mPartitionRecordCount, mPartitionSupplementaries.size(), mPartitionResolvedFragments.size(),
                candidateDuplicates.size(), mMaxPositionFragments, readPositionsStats);

        mPartitionResolvedFragments.clear();
        mPartitionSupplementaries.clear();
        mPartitionRecordCount = 0;
        mMaxPositionFragments = 0;

        if(setupNext)
        {
            mCurrentPartition.setStart(mCurrentPartition.end() + 1);
            mCurrentPartition.setEnd(mCurrentPartition.start() + mConfig.PartitionSize);
            mCurrentStrPartition = formChromosomePartition(mRegion.Chromosome, mCurrentPartition.start(), mConfig.PartitionSize);

            perfCounterStart();
        }

        System.gc();
    }

    private void processSamRecord(final SAMRecord read)
    {
        ++mTotalRecordCount;
        ++mPartitionRecordCount;

        if(mConfig.runReadChecks())
            mReadsProcessed.add(read);

        int readStart = read.getAlignmentStart();

        if(!readInSpecifiedRegions(read, mConfig.SpecificRegions, mConfig.SpecificChromosomes))
            return;

        if(readStart > mCurrentPartition.end())
        {
            onPartitionComplete(true);
        }

        if(mLogReadIds && mConfig.LogReadIds.contains(read.getReadName())) // debugging only
        {
            BM_LOGGER.debug("specific read: {}", readToString(read));
        }

        try
        {
            if(!mReadPositions.processRead(read))
            {
                // currently only supplementaries will not be stored against their initial fragment position or in an existing fragment
                processSupplementary(read);
            }
        }
        catch(Exception e)
        {
            BM_LOGGER.error("read({}) exception: {}", readToString(read), e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void processSupplementary(final SAMRecord read)
    {
        Fragment fragment = new Fragment(read);
        fragment.setRemotePartitions(mCurrentPartition);

        if(fragment.hasRemotePartitions())
            mPartitionSupplementaries.add(fragment);
        else
            mLocalGroupCombiner.localSupplementary(fragment, mCurrentStrPartition);
    }

    public void accept(final CandidateDuplicates candidateDuplicates)
    {
        List<Fragment> resolvedFragments = Lists.newArrayList();
        List<CandidateDuplicates> candidateDuplicatesList = Lists.newArrayList();

        mMaxPositionFragments = max(mMaxPositionFragments, candidateDuplicates.Fragments.size());

        classifyFragments(candidateDuplicates.Fragments, resolvedFragments, candidateDuplicatesList);

        if(mConfig.RunChecks)
            checkFragmentClassification(resolvedFragments, candidateDuplicatesList);

        // collapse unclear fragments with the same position since no benefit in keeping them separate
        if(candidateDuplicatesList.size() >= 1)
        {
            CandidateDuplicates first = candidateDuplicatesList.get(0);

            while(candidateDuplicatesList.size() > 1)
            {
                CandidateDuplicates next = candidateDuplicatesList.get(1);
                first.Fragments.addAll(next.Fragments);
                candidateDuplicatesList.remove(1);
            }
        }

        if(!resolvedFragments.isEmpty())
        {
            mStats.addDuplicateFrequencies(resolvedFragments);
            mRecordWriter.writeFragments(resolvedFragments);
        }

        for(Fragment fragment : resolvedFragments)
        {
            // the resolved fragment will be passed to the local group combiner and its state stored & applied

            // only store it for the remote group combiner if it has remote partitions
            fragment.setRemotePartitions(mCurrentPartition);

            if(fragment.hasRemotePartitions())
            {
                // need to combine any remote partitions from the 2 reads
                Fragment existingFragment = mPartitionResolvedFragments.get(fragment.id());

                if(existingFragment != null)
                {
                    fragment.reads().forEach(x -> existingFragment.addRead(x));
                    existingFragment.setRemotePartitions(mCurrentPartition);
                }
                else
                    mPartitionResolvedFragments.put(fragment.id(), fragment);
            }
        }

        // pass to the local combiner
        mLocalGroupCombiner.processPartitionFragments(mCurrentStrPartition, resolvedFragments, candidateDuplicatesList, Collections.EMPTY_LIST);
    }

    private void perfCounterStart()
    {
        if(mConfig.PerfDebug)
            mPerfCounter.start(format("%s:%s", mRegion.Chromosome, mCurrentPartition));
        else
            mPerfCounter.start();
    }
}
