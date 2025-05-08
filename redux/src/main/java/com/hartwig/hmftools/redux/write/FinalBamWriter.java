package com.hartwig.hmftools.redux.write;

import static java.lang.String.format;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.io.File;
import java.nio.file.Paths;
import java.util.List;
import java.util.Queue;

import com.google.common.collect.Lists;
import com.google.common.io.Files;
import com.hartwig.hmftools.common.bamops.BamOperations;
import com.hartwig.hmftools.redux.ReduxConfig;

public class FinalBamWriter extends Thread
{
    private final ReduxConfig mConfig;
    private final FileWriterCache mFileWriterCache;

    private Queue<PartitionInfo> mCompletedPartitionsQueue;
    private int mProcessedPartitions;

    private String mFinalBamFilename;

    // perf tracking
    private long mTotalWaitTimeMs;
    private long mTotalWriteTimeMs;

    public FinalBamWriter(final ReduxConfig config, final FileWriterCache fileWriterCache)
    {
        mConfig = config;
        mFileWriterCache = fileWriterCache;
        mCompletedPartitionsQueue = fileWriterCache.completedPartitionsQueue();
        mFinalBamFilename = "";
        mProcessedPartitions = 0;

        mTotalWaitTimeMs = 0;
        mTotalWriteTimeMs = 0;
    }

    @Override
    public void run()
    {
        int totalPartitions = mFileWriterCache.partitionCount();
        int completedPartitions = 0;

        List<PartitionInfo> pendingPartitions = Lists.newArrayList();

        long startWaitTimeMs = System.currentTimeMillis();

        while(true)
        {
            boolean hasNew = false;

            while(true)
            {
                PartitionInfo partitionInfo = mCompletedPartitionsQueue.poll();

                if(partitionInfo == null)
                    break;

                ++completedPartitions;
                addNewPartition(pendingPartitions, partitionInfo);
                hasNew = true;
            }

            if(!hasNew)
                continue;

            List<PartitionInfo> nextPartitions = findNextPartitions(pendingPartitions, mProcessedPartitions);

            if(nextPartitions == null)
                continue;

            RD_LOGGER.debug("BAM concat partitions completed({}/{}) ready({}) pending({}) processed({})",
                    completedPartitions, mFileWriterCache.partitionCount(),
                    nextPartitions.size(), pendingPartitions.size(), mProcessedPartitions);

            mTotalWaitTimeMs += System.currentTimeMillis() - startWaitTimeMs;

            processPartitions(nextPartitions);

            mProcessedPartitions += nextPartitions.size();

            if(mProcessedPartitions == totalPartitions)
                break;

            startWaitTimeMs = System.currentTimeMillis(); // start the wait timer again
        }
    }

    private static void addNewPartition(final List<PartitionInfo> pendingPartitions, final PartitionInfo newPartition)
    {
        // add in order to make checking more efficient
        int index = 0;
        while(index < pendingPartitions.size())
        {
            if(newPartition.regionIndex() < pendingPartitions.get(index).regionIndex())
                break;

            ++index;
        }

        pendingPartitions.add(index, newPartition);
    }

    public void logTimes()
    {
        RD_LOGGER.info(format("BAM concatenate writeTime(%.3f) waitTime(%.3f)",
                mTotalWriteTimeMs / 60000.0, mTotalWaitTimeMs / 60000.0));
    }

    private List<PartitionInfo> findNextPartitions(final List<PartitionInfo> pendingPartitions, int nextRequiredIndex)
    {
        if(pendingPartitions.isEmpty())
            return null;

        if(pendingPartitions.get(0).regionIndex() != nextRequiredIndex)
            return null;

        List<PartitionInfo> nextPartitions = Lists.newArrayList();

        while(!pendingPartitions.isEmpty())
        {
            PartitionInfo partitionInfo = pendingPartitions.get(0);

            if(partitionInfo.regionIndex() == nextRequiredIndex)
            {
                nextPartitions.add(partitionInfo);
                pendingPartitions.remove(0);
                ++nextRequiredIndex;
            }
            else
            {
                break;
            }
        }

        return nextPartitions;
    }

    private void processPartitions(final List<PartitionInfo> partitions)
    {
        // ensure all other writers are closed before attempting to read their records
        partitions.forEach(x -> x.bamWriter().close());

        if(partitions.get(0).regionIndex() == 0)
        {
            mFinalBamFilename = partitions.get(0).bamWriter().filename();
        }

        long startWriteTimeMs = System.currentTimeMillis();

        if(partitions.size() == 1)
        {
            RD_LOGGER.info("concatenating partition({}) to final BAM", partitions.get(0));
        }
        else
        {
            PartitionInfo first = partitions.get(0);
            PartitionInfo last = partitions.get(partitions.size() - 1);

            RD_LOGGER.info("concatenating {} partitions({}:{} - {}:{}) to final BAM",
                    partitions.size(), first.regions().get(0).Chromosome, first.regions().get(0).start(),
                    last.regions().get(last.regions().size() - 1).Chromosome, last.regions().get(last.regions().size() - 1).end());
        }

        // move temporarily before concatenating
        String tmpConcatBam = mConfig.OutputDir + "tmp_concat.bam";

        try
        {
            Files.move(new File(mFinalBamFilename), new File(tmpConcatBam));

            List<String> inputBams = Lists.newArrayList(tmpConcatBam);
            partitions.stream().filter(x -> x.regionIndex() > 0).forEach(x -> inputBams.add(x.bamWriter().filename()));

            // add on the unmapped reads BAM as the final partition BAM is being processed
            if(mProcessedPartitions + partitions.size() >= mFileWriterCache.partitionCount())
            {
                BamWriter fullyUnmappedBamWriter = mFileWriterCache.getFullUnmappedBamWriter();
                if(fullyUnmappedBamWriter != null)
                {
                    fullyUnmappedBamWriter.close();
                    inputBams.add(fullyUnmappedBamWriter.filename());
                }
            }

            if(!BamOperations.concatenateBams(mFileWriterCache.bamToolName(), mConfig.BamToolPath, mFinalBamFilename, inputBams, 1))
                System.exit(1);

            RD_LOGGER.debug("concatenate of {} BAMs complete", partitions.size());

            if(!mConfig.KeepInterimBams)
            {
                for(String inputBam : inputBams)
                {
                    java.nio.file.Files.deleteIfExists(Paths.get(inputBam));
                }
            }
        }
        catch(Exception e)
        {
            RD_LOGGER.error("failed to move final BAM for concatentation", e.toString());
            System.exit(1);
        }

        mTotalWriteTimeMs += System.currentTimeMillis() - startWriteTimeMs;
    }
}
