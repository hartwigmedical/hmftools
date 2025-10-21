package com.hartwig.hmftools.cobalt.metrics;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.cobalt.metrics.MetricsConfig.MIN_MAPPING_QUALITY;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.sv.SvUtils;

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
    private final Map<String, SAMRecord> mReadGroupMap;

    public PartitionReader(final MetricsConfig config, final Queue<Partition> partitions)
    {

        mPartitions = partitions;
        mPartitionCount = partitions.size();
        mCurrentPartition = null;

        mSamReader = SamReaderFactory.makeDefault().open(new File(config.BamFile));

        mBamSlicer = new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, false, false, false);
        mBamSlicer.setKeepUnmapped();

        mReadGroupMap = Maps.newHashMap();
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
                e.printStackTrace();
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
        postSliceProcess();
    }

    @VisibleForTesting
    protected void postSliceProcess()
    {
        // process overlapping groups
        for(SAMRecord read : mReadGroupMap.values())
        {
            processSingleRecord(read);
        }

        mReadGroupMap.clear();
    }

    private void processSamRecord(final SAMRecord read)
    {
        processSamRecord(read, true);
    }

    private void processSamRecord(final SAMRecord read, boolean removeFragments)
    {
        int readStart = read.getAlignmentStart();

        if(!mCurrentPartition.containsPosition(readStart))
        {
            return;
        }

        if(SvUtils.isDiscordant(read))
        {
            processSingleRecord(read);
            return;
        }

        SAMRecord mateRead = mReadGroupMap.get(read.getReadName());

        if(mateRead == null)
        {
            mReadGroupMap.put(read.getReadName(), read);
            return;
        }

        processFragment(read, mateRead);

        if(removeFragments)
        {
            mReadGroupMap.remove(read.getReadName());
        }
    }

    private void processSingleRecord(final SAMRecord read)
    {
        if(read.getMappingQuality() < MIN_MAPPING_QUALITY)
        {
            return;
        }

        int readPosStart = read.getAlignmentStart();
        int readPosEnd = read.getAlignmentEnd();

        if(read.getReadNegativeStrandFlag())
        {
            readPosEnd += CigarUtils.rightSoftClipLength(read);
        }
        else
        {
            readPosStart -= CigarUtils.leftSoftClipLength(read);
        }

        int length = readPosEnd - readPosStart;
        addFragmentData(readPosStart, length);
    }

    private void processFragment(final SAMRecord read, final SAMRecord mate)
    {
        if(read.getMappingQuality() < MIN_MAPPING_QUALITY || mate.getMappingQuality() < MIN_MAPPING_QUALITY)
        {
            return;
        }

        int readPosStart = read.getAlignmentStart();

        if(!read.getReadNegativeStrandFlag())
        {
            readPosStart -= CigarUtils.leftSoftClipLength(read);
        }

        int matePosStart = mate.getAlignmentStart();

        if(!mate.getReadNegativeStrandFlag())
        {
            matePosStart -= CigarUtils.leftSoftClipLength(mate);
        }

        int fragPosStart = min(readPosStart, matePosStart);

        int insertSize = abs(read.getInferredInsertSize());

        addFragmentData(fragPosStart, insertSize);
    }

    private void addFragmentData(int fragmentPosStart, int rawFragmentLength)
    {
        mCurrentPartition.recordFragment(fragmentPosStart, rawFragmentLength);
    }
}
