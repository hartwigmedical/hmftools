package com.hartwig.hmftools.cobalt.metrics;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.cobalt.metrics.MetricsConfig.MIN_MAPPING_QUALITY;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.sv.SvUtils;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionReader extends Thread
{
    private final Queue<Partition> mPartitions;
    private final int mPartitionCount;
    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private Partition mCurrentPartition;

    public PartitionReader(final MetricsConfig config, final Queue<Partition> partitions)
    {
        mPartitions = partitions;
        mPartitionCount = partitions.size();
        mCurrentPartition = null;

        mSamReader = SamReaderFactory.makeDefault().open(new File(config.BamFile));

        mBamSlicer = new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, false, false, false);
        mBamSlicer.setKeepUnmapped();
    }

    @Override
    public void run()
    {
        if(mPartitions.isEmpty())
        {
            return;
        }

        while(true)
        {
            try
            {
                mCurrentPartition = mPartitions.remove();

                int processedCount = mPartitionCount - mPartitions.size();

                if((processedCount % 100) == 0)
                {
                    CB_LOGGER.info("processed {} partitions", processedCount);
                }
                processPartition();
            }
            catch(NoSuchElementException e)
            {
                CB_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                CB_LOGGER.error("unexpected exception", e);
                System.exit(1);
            }
        }

        try
        {
            mSamReader.close();
        }
        catch(IOException e)
        {
            CB_LOGGER.error("failed to close bam file: {}", e.toString());
        }
    }

    private void processPartition()
    {
        CB_LOGGER.trace("processing region({})", mCurrentPartition);
        mBamSlicer.slice(mSamReader, mCurrentPartition, this::processSamRecord);
        CB_LOGGER.trace("completed region({})", mCurrentPartition);
    }

    private void processSamRecord(final SAMRecord read)
    {
        int readStart = read.getAlignmentStart();
        if(!mCurrentPartition.containsPosition(readStart))
        {
            return;
        }
        if(read.getMappingQuality() < MIN_MAPPING_QUALITY)
        {
            return;
        }
        if(SvUtils.isDiscordant(read))
        {
            return;
        }

        if(!SamRecordUtils.firstInPair(read))
        {
            return;
        }
        if(SamRecordUtils.mateUnmapped(read))
        {
            return;
        }
        List<CigarElement> cigarElements = read.getCigar().getCigarElements();
        if(cigarElements.size() != 1)
        {
            return;
        }
        int readPosStart = read.getAlignmentStart();
        mCurrentPartition.recordFragment(readPosStart, SamRecordUtils.inferredInsertSizeAbs(read));
    }
}