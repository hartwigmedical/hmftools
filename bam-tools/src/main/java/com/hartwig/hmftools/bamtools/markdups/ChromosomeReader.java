package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.DuplicateGroupUtils.findDuplicateFragments;
import static com.hartwig.hmftools.bamtools.markdups.DuplicateGroupUtils.processDuplicateGroups;
import static com.hartwig.hmftools.bamtools.markdups.FilterReadsType.readOutsideSpecifiedRegions;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.formChromosomePartition;
import static com.hartwig.hmftools.bamtools.markdups.FragmentUtils.readToString;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.function.Consumer;

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

public class ChromosomeReader implements Consumer<List<Fragment>>, Callable
{
    private final MarkDupsConfig mConfig;
    private final ChrBaseRegion mRegion;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final PartitionDataStore mPartitionDataStore;
    private final RecordWriter mRecordWriter;
    private final ReadPositionsCache mReadPositions;

    private BaseRegion mCurrentPartition;
    private String mCurrentStrPartition;
    private PartitionData mCurrentPartitionData;

    private final Map<String,List<Fragment>> mPendingIncompleteReads;

    private final boolean mLogReadIds;
    private int mPartitionRecordCount;
    private final Statistics mStats;
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

        mStats = new Statistics();

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
        mReadsProcessed = Sets.newHashSet();
        mPerfCounter = new PerformanceCounter("Slice");
    }

    public Set<SAMRecord> readsProcessed() { return mReadsProcessed; }
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

        BM_LOGGER.info("chromosome({}) complete, reads({})", mRegion.Chromosome, mStats.TotalReads);
    }

    private void onPartitionComplete(boolean setupNext)
    {
        mReadPositions.evictAll();

        processPendingIncompletes();

        mPerfCounter.stop();

        BM_LOGGER.debug("partition({}:{}) complete, reads({})", mRegion.Chromosome, mCurrentPartition, mPartitionRecordCount);

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

        if(mConfig.runReadChecks())
            mReadsProcessed.add(read);

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
            BM_LOGGER.debug("specific read: {}", readToString(read));
        }

        if(!mConfig.NoMateCigar && read.getReadPairedFlag() && !read.getMateUnmappedFlag() && !read.getSupplementaryAlignmentFlag()
        && !read.hasAttribute(MATE_CIGAR_ATTRIBUTE))
        {
            ++mStats.MissingMateCigar;
        }

        try
        {
            if(read.getSupplementaryAlignmentFlag() || !mReadPositions.processRead(read))
            {
                ++mStats.Incomplete;

                String basePartition = Fragment.getBasePartition(read, mConfig.PartitionSize);
                Fragment fragment = new Fragment(read);

                if(basePartition == null)
                {
                    // mate or supp is on a non-human chromsome, meaning it won't be retrieved - so write this immediately
                    mRecordWriter.writeFragment(fragment);
                    return;
                }

                processIncompleteRead(fragment, basePartition);
            }
        }
        catch(Exception e)
        {
            BM_LOGGER.error("read({}) exception: {}", readToString(read), e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void processIncompleteRead(final Fragment fragment, final String basePartition)
    {
        if(basePartition.equals(mCurrentStrPartition))
        {
            // String callerInfo = format("%s_incomplete_%d", mCurrentStrPartition, abs(fragment.initialPosition()));
            List<Fragment> resolvedFragments = mCurrentPartitionData.processIncompleteFragment(fragment);

            if(fragment.status().isResolved())
                mRecordWriter.writeFragment(fragment);

            if(resolvedFragments != null)
                mRecordWriter.writeFragments(resolvedFragments);
        }
        else
        {
            ++mStats.InterPartition;

            // cache this read and send through as groups when the partition is complete
            List<Fragment> pendingFragments = mPendingIncompleteReads.get(basePartition);

            if(pendingFragments == null)
            {
                pendingFragments = Lists.newArrayList();
                mPendingIncompleteReads.put(basePartition, pendingFragments);
            }

            pendingFragments.add(fragment);
        }
    }

    private void processPendingIncompletes()
    {
        if(mPendingIncompleteReads.isEmpty())
            return;

        if(mPendingIncompleteReads.size() > 100)
        {
            BM_LOGGER.debug("partition({}:{}) processing {} pending incomplete fragments",
                    mRegion.Chromosome, mCurrentPartition, mPendingIncompleteReads.values().stream().mapToInt(x -> x.size()).sum());
        }

        for(Map.Entry<String,List<Fragment>> entry : mPendingIncompleteReads.entrySet())
        {
            String basePartition = entry.getKey();
            List<Fragment> fragments = entry.getValue();

            PartitionData partitionData = mPartitionDataStore.getOrCreatePartitionData(basePartition);

            List<Fragment> resolvedFragments = partitionData.processIncompleteFragments(fragments);

            if(resolvedFragments != null)
            {
                mRecordWriter.writeFragments(resolvedFragments);
            }
        }

        mPendingIncompleteReads.clear();
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

        processDuplicateGroups(duplicateGroups, mConfig.UMIs, mStats);

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

        startTimeMs = System.currentTimeMillis();

        mCurrentPartitionData.processPrimaryFragments(resolvedFragments, candidateDuplicatesList);

        double timeTakenSec = (System.currentTimeMillis() - startTimeMs) / 1000.0;

        if(timeTakenSec >= 1.0)
        {
            BM_LOGGER.debug("position({}:{}) fragments({}) partition processing time({})",
                    mRegion.Chromosome, position, posFragmentCount, format("%.1fs", timeTakenSec));
        }

        if(!resolvedFragments.isEmpty())
        {
            mRecordWriter.writeFragments(resolvedFragments);

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
}
