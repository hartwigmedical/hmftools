package com.hartwig.hmftools.bamtools.tester;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import java.io.File;
import java.util.List;
import java.util.NoSuchElementException;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.perf.PerformanceCounter;
import com.hartwig.hmftools.common.perf.TaskQueue;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionTester extends Thread
{
    private final TesterConfig mConfig;
    private final TaskQueue mRegions;

    private final SamReader mSamReader;

    private ChrBaseRegion mCurrentRegion;

    private long mReadCount;

    private final List<PerformanceCounter> mPerfCounters;

    private final PerformanceCounter mQueryInternalPerfCounter;
    private final PerformanceCounter mRecordParsePerfCounter;

    private final boolean mLogReadIds;

    private enum PerfCounterType
    {
        QueryInterval,
        RecordParse;
    }

    public PartitionTester(final TesterConfig config, final TaskQueue regions)
    {
        mConfig = config;
        mRegions = regions;

        mCurrentRegion = null;

        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));

        mReadCount = 0;

        mPerfCounters = Lists.newArrayList();

        for(PerfCounterType perfCounterType : PerfCounterType.values())
        {
            mPerfCounters.add(new PerformanceCounter(perfCounterType.toString()));
        }

        mQueryInternalPerfCounter = mPerfCounters.get(PerfCounterType.QueryInterval.ordinal());
        mRecordParsePerfCounter = mPerfCounters.get(PerfCounterType.RecordParse.ordinal());

        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    public List<PerformanceCounter> perfCounters() { return mPerfCounters; }
    public long totalReads() { return mReadCount;}

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                ChrBaseRegion region = (ChrBaseRegion)mRegions.removeItem();

                processRegion(region);
            }
            catch(NoSuchElementException e)
            {
                BT_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    private static final int RECORD_PERF_COUNT = 100;

    private void processRegion(final ChrBaseRegion region)
    {
        mCurrentRegion = region;

        QueryInterval[] queryIntervals = BamSlicer.createIntervals(List.of(region), mSamReader.getFileHeader());

        mQueryInternalPerfCounter.start();
        SAMRecordIterator iterator = mSamReader.queryOverlapping(queryIntervals);
        mQueryInternalPerfCounter.stop();

        int recordCount = 0;
        int nextRecordCount = RECORD_PERF_COUNT;

        mRecordParsePerfCounter.start();

        while(iterator.hasNext())
        {
            processRecord(iterator.next());

            ++recordCount;

            if(recordCount >= nextRecordCount)
            {
                nextRecordCount += RECORD_PERF_COUNT;

                mRecordParsePerfCounter.stop();
                mRecordParsePerfCounter.start();
            }
        }

        iterator.close();

        mRecordParsePerfCounter.stop();

        if(mConfig.SpecificChrRegions.hasFilters())
        {
            BT_LOGGER.debug("region({}) processed {} reads", mCurrentRegion, recordCount);
        }
    }

    private void processRecord(final SAMRecord record)
    {
        // extract fields but nothing more
        ++mReadCount;

        int readStart = record.getAlignmentStart();
        int readEnd = record.getAlignmentEnd();

        if(!mCurrentRegion.containsPosition(readStart))
            return;

        if(mLogReadIds && mConfig.LogReadIds.contains(record.getReadName()))
        {
            BT_LOGGER.debug("specific read({}) coords({}:{}-{})", record.getContig(), readStart, readEnd);
        }
    }
}
