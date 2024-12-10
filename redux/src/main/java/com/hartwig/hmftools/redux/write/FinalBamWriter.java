package com.hartwig.hmftools.redux.write;

import static java.lang.String.format;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.List;
import java.util.Queue;

import com.google.common.collect.Lists;
import com.google.common.io.Files;
import com.hartwig.hmftools.common.bamops.BamOperations;
import com.hartwig.hmftools.redux.ReduxConfig;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class FinalBamWriter extends Thread
{
    private final ReduxConfig mConfig;
    private FileWriterCache mFileWriterCache;
    private Queue<PartitionInfo> mCompletedPartitionsQueue;
    private int mProcessedPartitions;

    private String mFinalBamFilename;
    private SAMFileWriter mBamWriter;

    // perf tracking
    private long mTotalWaitTimeMs;
    private long mTotalWriteTimeMs;

    public FinalBamWriter(final ReduxConfig config, final FileWriterCache fileWriterCache)
    {
        mConfig = config;
        mFileWriterCache = fileWriterCache;
        mCompletedPartitionsQueue = fileWriterCache.completedPartitionsQueue();
        mFinalBamFilename = "";
        mBamWriter = null;
        mProcessedPartitions = 0;

        mTotalWaitTimeMs = 0;
        mTotalWriteTimeMs = 0;
    }

    @Override
    public void run()
    {
        int totalPartitions = mFileWriterCache.partitionCount();

        List<PartitionInfo> pendingPartitions = Lists.newArrayList();

        long startWaitTimeMs = System.currentTimeMillis();
        int completedPartitionReaders = 0;

        while(true)
        {
            PartitionInfo partitionInfo = mCompletedPartitionsQueue.poll();

            if(partitionInfo == null)
                continue;

            ++completedPartitionReaders;

            if((completedPartitionReaders % 10) == 0)
            {
                RD_LOGGER.info("completed {} partition readers", completedPartitionReaders);
            }

            // add in order to make checking more efficient
            int index = 0;
            while(index < pendingPartitions.size())
            {
                if(partitionInfo.regionIndex() < pendingPartitions.get(index).regionIndex())
                    break;

                ++index;
            }

            pendingPartitions.add(index, partitionInfo);

            PartitionInfo nextPartitionInfo = findNextPartition(pendingPartitions, mProcessedPartitions);
            boolean processedNew = false;

            while(nextPartitionInfo != null)
            {
                if(!processedNew)
                {
                    processedNew = true;
                    mTotalWaitTimeMs += System.currentTimeMillis() - startWaitTimeMs;
                }

                processPartition(nextPartitionInfo);
                ++mProcessedPartitions;

                RD_LOGGER.info("final partition processing completed({}/{}) pending({})",
                        mProcessedPartitions, mFileWriterCache.partitionCount(), pendingPartitions.size());

                if(mProcessedPartitions == totalPartitions)
                    break;

                nextPartitionInfo = findNextPartition(pendingPartitions, mProcessedPartitions);
            }

            if(mProcessedPartitions == totalPartitions)
                break;

            if(processedNew)
                startWaitTimeMs = System.currentTimeMillis();
        }
    }

    public void logTimes()
    {
        RD_LOGGER.info(format("final BAM writeTime(%.3f) waitTime(%.3f)",
                mTotalWriteTimeMs / 60000.0, mTotalWaitTimeMs / 60000.0));
    }

    private PartitionInfo findNextPartition(final List<PartitionInfo> pendingPartitions, int nextRequiredIndex)
    {
        if(pendingPartitions.isEmpty())
            return null;

        // partitions are cached in order so just check the first
        PartitionInfo partitionInfo = pendingPartitions.get(0);

        if(partitionInfo.regionIndex() == nextRequiredIndex)
        {
            pendingPartitions.remove(0);
            return partitionInfo;
        }

        return null;
    }

    private void processPartition(final PartitionInfo partitionInfo)
    {
        // ensure all other writers are closed before attempting to read their records
        partitionInfo.bamWriter().close();

        if(partitionInfo.regionIndex() == 0)
        {
            mFinalBamFilename = partitionInfo.bamWriter().filename();
            return;
        }

        long startWriteTimeMs = System.currentTimeMillis();

        RD_LOGGER.info("writing partition({}) to final BAM", partitionInfo);

        // move temporarily before concatenating
        String tmpConcatBam = mConfig.OutputBam + "tmp_concat.bam";

        try
        {
            Files.move(new File(mFinalBamFilename), new File(tmpConcatBam));

            List<String> inputBams = List.of(tmpConcatBam, partitionInfo.bamWriter().filename());

            if(!BamOperations.concatenateBams(mFileWriterCache.bamToolName(), mConfig.BamToolPath, mFinalBamFilename, inputBams, 1))
                System.exit(1);

            java.nio.file.Files.deleteIfExists(Paths.get(tmpConcatBam));
        }
        catch(Exception e)
        {
            RD_LOGGER.error("failed to move final BAM for concatentation", e.toString());
            System.exit(1);
        }

        mTotalWriteTimeMs += System.currentTimeMillis() - startWriteTimeMs;
    }

    private void processPartitionManual(final PartitionInfo partitionInfo)
    {
        long startWriteTimeMs = System.currentTimeMillis();

        RD_LOGGER.info("writing partition({}) to final BAM", partitionInfo);

        String bamFilename = partitionInfo.bamWriter().filename();

        if(partitionInfo.regionIndex() == 0)
        {
            // writing will be to this first BAM, so no need to write it again
            mBamWriter = partitionInfo.bamWriter().samFileWriter();
            return;
        }
        else
        {
            // ensure all other writers are closed before attempting to read their records
            partitionInfo.bamWriter().close();
        }

        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(bamFilename));

        SAMRecordIterator iterator = samReader.iterator();

        while(iterator.hasNext())
        {
            SAMRecord record = iterator.next();
            mBamWriter.addAlignment(record);
        }

        mTotalWriteTimeMs += System.currentTimeMillis() - startWriteTimeMs;
    }
}
