package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.classifyFragments;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.formChromosomePartition;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.readInSpecifiedRegions;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.readToString;

import java.io.File;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

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

    private final BaseRegion mCurrentPartition;
    private String mCurrentStrPartition;
    private PartitionData mCurrentPartitionData;
    private final List<Fragment> mPendingIncompleteReads;

    private final boolean mLogReadIds;
    private int mTotalRecordCount;
    private int mPartitionRecordCount;
    private final DuplicateStats mStats;
    private final PerformanceCounter mPerfCounter;

    private final Set<SAMRecord> mReadsProcessed; // for debug only when running checks

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

        mReadPositions = new ReadPositionsCache(region.Chromosome, config.BufferSize, this);

        if(!mConfig.SpecificRegions.isEmpty())
        {
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

        mPendingIncompleteReads = Lists.newArrayList();

        mTotalRecordCount = 0;
        mPartitionRecordCount = 0;

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
        mReadPositions.evictAll();

        processIncompletes();

        mStats.ReadCount += mPartitionRecordCount;

        mPerfCounter.stop();

        BM_LOGGER.debug("partition({}:{}) complete, reads({})", mRegion.Chromosome, mCurrentPartition, mPartitionRecordCount);

        mCurrentPartitionData.logCacheCounts();

        mPartitionRecordCount = 0;

        if(setupNext)
        {
            mCurrentPartition.setStart(mCurrentPartition.end() + 1);
            mCurrentPartition.setEnd(mCurrentPartition.start() + mConfig.PartitionSize - 1);
            mCurrentStrPartition = formChromosomePartition(mRegion.Chromosome, mCurrentPartition.start(), mConfig.PartitionSize);
            mCurrentPartitionData = mPartitionDataStore.getOrCreatePartitionData(mCurrentStrPartition);

            perfCounterStart();
        }

        System.gc();
    }

    private void processSamRecord(final SAMRecord read)
    {
        int readStart = read.getAlignmentStart();

        if(!readInSpecifiedRegions(read, mConfig.SpecificRegions, mConfig.SpecificChromosomes, false))
            return;

        ++mTotalRecordCount;
        ++mPartitionRecordCount;

        if(mConfig.runReadChecks())
            mReadsProcessed.add(read);

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
            if(read.getSupplementaryAlignmentFlag() || !mReadPositions.processRead(read))
            {
                Fragment fragment = new Fragment(read);
                processIncompleteRead(fragment);
            }
        }
        catch(IllegalMonitorStateException e)
        {
            // cache this read and try again later on
            mPendingIncompleteReads.add(new Fragment(read));
        }
        catch(Exception e)
        {
            BM_LOGGER.error("read({}) exception: {}", readToString(read), e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void processIncompleteRead(final Fragment fragment)
    {
        String basePartition = Fragment.getBasePartition(fragment.reads().get(0), mConfig.PartitionSize);

        if(basePartition == null)
        {
            // mate or supp is on a non-human chromsome, meaning it won't be retrieved - so write this immediately
            mRecordWriter.writeFragment(fragment);
            return;
        }

        PartitionData partitionData = basePartition.equals(mCurrentStrPartition) ?
                mCurrentPartitionData : mPartitionDataStore.getOrCreatePartitionData(basePartition);

        List<Fragment> resolvedFragments = partitionData.processIncompleteFragment(fragment);

        if(fragment.status().isResolved())
            mRecordWriter.writeFragment(fragment);

        if(resolvedFragments != null)
            mRecordWriter.writeFragments(resolvedFragments);
    }

    private void processIncompletes()
    {
        if(mPendingIncompleteReads.isEmpty())
            return;

        BM_LOGGER.debug("partition({}:{}) processing {} pending incomplete fragment",
                mRegion.Chromosome, mCurrentPartition, mPendingIncompleteReads.size());

        for(Fragment fragment : mPendingIncompleteReads)
        {
            try
            {
                processIncompleteRead(fragment);
            }
            catch(Exception e)
            {
                BM_LOGGER.error("failed to process pending incomplete read({})", readToString(fragment.reads().get(0)));
                mRecordWriter.writeFragment(fragment);
            }
        }

        mPendingIncompleteReads.clear();
    }

    public void accept(final List<Fragment> positionFragments)
    {
        List<Fragment> resolvedFragments = Lists.newArrayList();
        List<CandidateDuplicates> candidateDuplicatesList = Lists.newArrayList();

        int posFragmentCount = positionFragments.size();
        boolean logDetails = mConfig.PerfDebug && posFragmentCount > 10000;
        long startTimeMs = 0;
        int position = 0;

        if(logDetails)
        {
            position = positionFragments.get(0).initialPosition();
            startTimeMs = System.currentTimeMillis();
        }

        classifyFragments(positionFragments, resolvedFragments, candidateDuplicatesList);

        if(logDetails)
        {
            double timeTakenSec = (System.currentTimeMillis() - startTimeMs) / 1000.0;

            if(timeTakenSec >= 1.0)
            {
                BM_LOGGER.debug("position({}:{}) fragments({}) resolved({}) candidates({}) time({})",
                        mRegion.Chromosome, position, posFragmentCount, resolvedFragments.size(),
                        candidateDuplicatesList.stream().mapToInt(x -> x.fragmentCount()).sum(),
                        format("%.1fs", timeTakenSec));
            }
        }

        if(logDetails)
            startTimeMs = System.currentTimeMillis();

        mCurrentPartitionData.processPrimaryFragments(resolvedFragments, candidateDuplicatesList);

        if(logDetails)
        {
            double timeTakenSec = (System.currentTimeMillis() - startTimeMs) / 1000.0;

            if(timeTakenSec >= 1.0)
            {
                BM_LOGGER.debug("position({}:{}) fragments({}) partition processing time({})",
                        mRegion.Chromosome, position, posFragmentCount, format("%.1fs", timeTakenSec));
            }
        }

        if(!resolvedFragments.isEmpty())
        {
            // no longer ordered in duplicate groups - either re-evaluate or store by coordinate key or track some other way
            mStats.addDuplicateInfo(resolvedFragments);
            mRecordWriter.writeFragments(resolvedFragments);
        }
    }

    private void perfCounterStart()
    {
        if(mConfig.PerfDebug)
            mPerfCounter.start(format("%s:%s", mRegion.Chromosome, mCurrentPartition));
        else
            mPerfCounter.start();
    }
}
