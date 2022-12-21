package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.classifyFragments;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.formChromosomePartition;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.reconcileFragments;

import java.io.File;
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

public class ChromosomeReader implements Consumer<PositionFragments>, Callable
{
    private final MarkDupsConfig mConfig;
    private final ChrBaseRegion mRegion;
    private final BaseRegion mCurrentPartition;
    private String mCurrentStrPartition;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final GroupCombiner mRemoteGroupCombiner;
    private final RecordWriter mRecordWriter;
    private final ReadPositionArray mReadPositions;

    // private final GroupCombiner mLocalGroupCombiner;
    private final Map<String,Fragment> mPartitionSupplementaries;
    private final Map<String,Fragment> mPartitionResolvedFragments;
    private final List<PositionFragments> mPartitionIncompletePositionFragments;

    private final boolean mLogReadIds;
    private int mTotalRecordCount;
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

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(0, true, true, true);
        mBamSlicer.setKeepUnmapped();

        mRecordWriter = recordWriter;
        mReadPositions = new ReadPositionArray(region.Chromosome, config.BufferSize, this);

        // mLocalGroupCombiner = new GroupCombiner(mRecordWriter);
        mPartitionSupplementaries = Maps.newHashMap();
        mPartitionResolvedFragments = Maps.newHashMap();
        mPartitionIncompletePositionFragments = Lists.newArrayList();

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

        /*
        for(List<SAMRecord> supplementaries : mSupplementaries.values())
        {
            supplementaries.forEach(x -> mRecordWriter.writeRecord(x, FragmentStatus.NONE));
        }
        */

        onPartitionComplete(false);

        BM_LOGGER.info("chromosome({}) complete, reads({})", mRegion.Chromosome, mTotalRecordCount);
        mReadsProcessed.clear();
    }

    private void onPartitionComplete(boolean setupNext)
    {
        mReadPositions.evictAll();

        // first try to assign any supplementaries in this partition to their fragments
        reconcileFragments(mPartitionSupplementaries, mPartitionResolvedFragments, mPartitionIncompletePositionFragments);

        List<Fragment> remoteSupplementaries = mPartitionSupplementaries.values().stream().collect(Collectors.toList());
        List<Fragment> resolvedFragments = mPartitionResolvedFragments.values().stream().collect(Collectors.toList());

        int i = 0;
        while(i < resolvedFragments.size())
        {
            Fragment fragment = resolvedFragments.get(i);

            if(fragment.allReadsPresent())
            {
                mRecordWriter.writeFragment(fragment);
                resolvedFragments.remove(i);
            }
            else
            {
                ++i;
            }
        }

        mRemoteGroupCombiner.processPartitionFragments(
                mCurrentStrPartition, resolvedFragments, mPartitionIncompletePositionFragments, remoteSupplementaries);

        mPerfCounter.stop();

        BM_LOGGER.debug("partition({}:{}) complete, cached(supps={} resolved={} incompletes={})",
                mRegion.Chromosome, mCurrentPartition, mPartitionSupplementaries.size(), mPartitionResolvedFragments.size(),
                mPartitionIncompletePositionFragments.size());

        if(setupNext)
        {
            mCurrentPartition.setStart(mCurrentPartition.end() + 1);
            mCurrentPartition.setEnd(mCurrentPartition.start() + mConfig.PartitionSize);
            mCurrentStrPartition = formChromosomePartition(mRegion.Chromosome, mCurrentPartition.start(), mConfig.PartitionSize);

            perfCounterStart();
        }

        mPartitionResolvedFragments.clear();
        mPartitionIncompletePositionFragments.clear();
        mPartitionSupplementaries.clear();

        System.gc();
    }

    private void processSamRecord(final SAMRecord read)
    {
        ++mTotalRecordCount;

        if(BM_LOGGER.isTraceEnabled())
            mReadsProcessed.add(read);

        int readStart = read.getAlignmentStart();

        if(readStart > mCurrentPartition.end())
        {
            onPartitionComplete(true);
        }

        // if(!mCurrentPartition.containsPosition(readStart))
        //    return;

        if(mLogReadIds) // debugging only
        {
            if(mConfig.LogReadIds.contains(read.getReadName()))
            {
                BM_LOGGER.debug("specific readId({})", read.getReadName());
            }
        }

        if(!mReadPositions.processRead(read))
        {
            // currently only supplementaries will not be stored against their initial fragment position or in an existing fragment
            processSupplementary(read);
        }
    }

    private void processSupplementary(final SAMRecord read)
    {
        Fragment fragment = new Fragment(read);
        mPartitionSupplementaries.put(fragment.id(), fragment);
    }

    public void accept(final PositionFragments positionFragments)
    {
        List<Fragment> resolvedFragments = Lists.newArrayList();
        List<PositionFragments> incompletePositionFragments = Lists.newArrayList();

        classifyFragments(positionFragments, resolvedFragments, incompletePositionFragments);

        int unclearFragments = incompletePositionFragments.stream().mapToInt(x -> x.Fragments.size()).sum();
        if(positionFragments.Fragments.size() != resolvedFragments.size() + unclearFragments)
        {
            BM_LOGGER.error("failed to classify all fragments");
        }

        for(Fragment fragment : resolvedFragments)
        {
            if(!fragment.status().isResolved())
            {
                BM_LOGGER.error("fragment({}) incorrectly in resolved list", fragment);
            }
        }

        resolvedFragments.forEach(x -> mPartitionResolvedFragments.put(x.id(), x));
        mPartitionIncompletePositionFragments.addAll(incompletePositionFragments);
    }

    private void perfCounterStart()
    {
        if(mConfig.PerfDebug)
            mPerfCounter.start(format("%s:%s", mRegion.Chromosome, mCurrentPartition));
        else
            mPerfCounter.start();
    }
}
