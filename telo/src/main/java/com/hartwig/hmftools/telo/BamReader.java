package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;

import java.util.NoSuchElementException;
import java.util.Queue;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class BamReader implements Runnable
{
    private final TeloConfig mConfig;

    private final Queue<ChrBaseRegion> mBaseRegionQ;
    private final Queue<TelBamRecord> mTelBamRecordQ;
    private final Set<String> mIncompleteReadNames;

    private int mReadCount = 0;

    // one sam reader per thread
    private SamReader mSamReader;

    public BamReader(final TeloConfig config, Queue<ChrBaseRegion> baseRegionQ, Queue<TelBamRecord> telBamRecordQ, Set<String> incompleteReadNames)
    {
        mConfig = config;
        mBaseRegionQ = baseRegionQ;
        mTelBamRecordQ = telBamRecordQ;
        mIncompleteReadNames = incompleteReadNames;
        mSamReader = TeloUtils.openSamReader(mConfig);
    }

    @Override
    public void run()
    {
        while(true)
        {
            ChrBaseRegion baseRegion;
            try
            {
                baseRegion = mBaseRegionQ.remove();
            }
            catch (NoSuchElementException e)
            {
                // finished processing
                break;
            }
            findTelomereContent(baseRegion);
        }
    }

    public void findTelomereContent(ChrBaseRegion baseRegion)
    {
        if(!baseRegion.equals(TeloConstants.UNMAPPED_BASE_REGION))
        {
            processBamByRegion(baseRegion);
        }
        else
        {
            processUnmappedReads();
        }
    }

    private void processBamByRegion(ChrBaseRegion baseRegion)
    {
        TE_LOGGER.info("processing region({})", baseRegion.toString());
        mReadCount = 0;

        BamSlicer bamSlicer = new BamSlicer(1, false, false, false);
        bamSlicer.slice(mSamReader, Lists.newArrayList(baseRegion), this::processReadRecord);

        TE_LOGGER.info("processed region({}) read count({})", baseRegion.toString(), mReadCount);
    }

    private void processUnmappedReads()
    {
        TE_LOGGER.info("processing unmapped reads, incomplete={}", mIncompleteReadNames.size());
        mReadCount = 0;

        try (final SAMRecordIterator iterator = mSamReader.queryUnmapped())
        {
            while (iterator.hasNext())
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
        // we discard any supplementary / secondary alignments
        if(record.isSecondaryOrSupplementary())
            return;

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
