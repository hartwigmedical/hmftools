package com.hartwig.hmftools.teal;

import static com.hartwig.hmftools.teal.TeloConfig.TE_LOGGER;

import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.Set;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class BamReader implements Runnable
{
    public static class Task
    {
        private ChrBaseRegion mBaseRegion;
        private boolean mQueryUnmapped;
        public ChrBaseRegion baseRegion() { return mBaseRegion; }
        public boolean queryUnmapped() { return mQueryUnmapped; }
        public static Task fromBaseRegion(ChrBaseRegion baseRegion)
        {
            Task t = new Task();
            t.mBaseRegion = baseRegion;
            return t;
        }
        public static Task fromQueryUnmapped()
        {
            Task t = new Task();
            t.mQueryUnmapped = true;
            return t;
        }
    }

    private final TeloConfig mConfig;

    private final Queue<Task> mTaskQ;
    private final Queue<TelBamRecord> mTelBamRecordQ;
    private final Set<String> mIncompleteReadNames;

    private int mReadCount = 0;

    // one sam reader per thread
    private SamReader mSamReader;

    public BamReader(final TeloConfig config, Queue<Task> taskQ, Queue<TelBamRecord> telBamRecordQ, Set<String> incompleteReadNames)
    {
        mConfig = config;
        mTaskQ = taskQ;
        mTelBamRecordQ = telBamRecordQ;
        mIncompleteReadNames = incompleteReadNames;
        mSamReader = TeloUtils.openSamReader(mConfig);
    }

    @Override
    public void run()
    {
        while(true)
        {
            Task task;

            try
            {
                task = mTaskQ.remove();
            }
            catch (NoSuchElementException e)
            {
                // finished processing
                break;
            }
            findTelomereContent(task);
        }
    }

    public void findTelomereContent(Task task)
    {
        if(task.queryUnmapped())
        {
            processUnmappedReads();
        }
        if(task.baseRegion() != null)
        {
            processBamByRegion(task.baseRegion());
        }
    }

    // Important note:
    // here we must use SamReader::query. I have tested that this function would return unmapped read that has
    // chromosome and start position set. We need that cause unfortunately the queryUnmapped function is not going to return
    // them. See https://github.com/samtools/htsjdk/issues/278
    private void processBamByRegion(ChrBaseRegion baseRegion)
    {
        TE_LOGGER.info("processing region({})", baseRegion.toString());
        mReadCount = 0;

        // do not change the follow line to use functions other than query without testing that unmapped reads are returned
        try (final SAMRecordIterator iterator = mSamReader.query(baseRegion.Chromosome, baseRegion.start(), baseRegion.end(), false))
        {
            while (iterator.hasNext())
            {
                final SAMRecord record = iterator.next();
                processReadRecord(record);
            }
        }

        TE_LOGGER.info("processed region({}) read count({})", baseRegion.toString(), mReadCount);
    }

    // Important note:
    // queryUnmapped will not return unmapped read that has chromosome / pos start set. These are
    // unmapped read where the mate is mapped.
    // see https://github.com/samtools/htsjdk/issues/278
    private void processUnmappedReads()
    {
        TE_LOGGER.info("processing unmapped reads, incomplete={}", mIncompleteReadNames.size());
        mReadCount = 0;

        try(final SAMRecordIterator iterator = mSamReader.queryUnmapped())
        {
            while(iterator.hasNext())
            {
                final SAMRecord record = iterator.next();
                processReadRecord(record);
            }
        }

        TE_LOGGER.info("processed {} unmapped reads", mReadCount);
    }

    private boolean hasTelomericContent(@NotNull final SAMRecord record)
    {
        return TeloUtils.hasTelomericContent(record.getReadString());
    }

    private void processReadRecord(@NotNull final SAMRecord record)
    {
        boolean hasTeloContent = hasTelomericContent(record);

        if (!hasTeloContent && !mIncompleteReadNames.contains(record.getReadName()))
        {
            return;
        }

        ++mReadCount;

        // push it on to the queue
        TelBamRecord telBamRecord = new TelBamRecord();
        telBamRecord.samRecord = record;
        telBamRecord.hasTeloContent = hasTeloContent;

        mTelBamRecordQ.add(telBamRecord);
    }
}
