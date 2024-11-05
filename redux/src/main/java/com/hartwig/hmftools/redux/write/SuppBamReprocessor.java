package com.hartwig.hmftools.redux.write;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.redux.PartitionDataStore;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.common.DuplicateGroup;
import com.hartwig.hmftools.redux.common.Fragment;
import com.hartwig.hmftools.redux.common.PartitionData;
import com.hartwig.hmftools.redux.common.PartitionResults;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class SuppBamReprocessor implements Callable
{
    private final ReduxConfig mConfig;
    private final BamWriter mBamWriter;
    private final PartitionDataStore mPartitionDataStore;
    private final ConsensusReads mConsensusReads;

    private final SamReader mBamReader;
    private final String mBamFilename;
    private final Map<String,List<SAMRecord>> mPendingReads;
    private int mCacheCount;
    private int mProcessedCount;

    private static final int READ_BATCH_SIZE = 10000;
    private static final int READ_LOG_COUNT = 100000;

    public SuppBamReprocessor(
            final ReduxConfig config, final String bamFilename, final BamWriter bamWriter,
            final PartitionDataStore partitionDataStore, final ConsensusReads consensusReads)
    {
        mConfig = config;
        mBamWriter = bamWriter;
        mBamFilename = bamFilename;
        mPartitionDataStore = partitionDataStore;
        mConsensusReads = consensusReads;

        mPendingReads = Maps.newHashMap();
        mCacheCount = 0;
        mProcessedCount = 0;

        mBamReader = SamReaderFactory.makeDefault().referenceSequence(new File(config.RefGenomeFile)).open(new File(bamFilename));
    }

    @Override
    public Long call()
    {
        SAMRecordIterator readIterator = mBamReader.iterator();

        while(readIterator.hasNext())
        {
            processSupplementary(readIterator.next());
        }

        processCachedReads();
        mPendingReads.clear();

        if(mProcessedCount > 0)
        {
            RD_LOGGER.debug("finished reprocessing {} supplementaries from {}", mProcessedCount, mBamFilename);
        }

        return (long) 0;
    }

    private void processSupplementary(final SAMRecord read)
    {
        String chrPartition = Fragment.getBasePartition(read, mConfig.PartitionSize);

        List<SAMRecord> pendingFragments = mPendingReads.get(chrPartition);

        if(pendingFragments == null)
        {
            pendingFragments = Lists.newArrayList();
            mPendingReads.put(chrPartition, pendingFragments);
        }

        pendingFragments.add(read);

        ++mCacheCount;

        if(mCacheCount >= READ_BATCH_SIZE)
        {
            processCachedReads();
            mCacheCount = 0;
        }

        ++mProcessedCount;

        if((mProcessedCount % READ_LOG_COUNT) == 0)
        {
            RD_LOGGER.debug("reprocessed {} supplementaries from {}", mProcessedCount, mBamFilename);
        }
    }

    private void processCachedReads()
    {
        // RD_LOGGER.trace("reprocessing {} supplementary reads", mCacheCount);

        for(Map.Entry<String,List<SAMRecord>> entry : mPendingReads.entrySet())
        {
            String basePartition = entry.getKey();
            List<SAMRecord> reads = entry.getValue();

            PartitionData partitionData = mPartitionDataStore.getOrCreatePartitionData(basePartition);

            PartitionResults partitionResults = partitionData.processIncompleteFragments(reads, false);

            if(partitionResults == null)
                return;

            if(partitionResults.umiGroups() != null)
            {
                for(DuplicateGroup duplicateGroup : partitionResults.umiGroups())
                {
                    // form consensus reads for any complete read leg groups and write reads
                    List<SAMRecord> completeReads = duplicateGroup.popCompletedReads(mConsensusReads, false);
                    mBamWriter.writeDuplicateGroup(duplicateGroup, completeReads);
                }
            }

            if(partitionResults.resolvedFragments() != null)
                mBamWriter.writeFragments(partitionResults.resolvedFragments(), true);

            reads.clear();
        }
    }

    public static void reprocessSupplementaries(
            final ReduxConfig config, final BamWriter bamWriter, final PartitionDataStore partitionDataStore,
            final List<SuppBamWriter> suppBamWriters, final ConsensusReads consensusReads)
    {
        suppBamWriters.forEach(x -> x.close());

        long totalSupplementaries = suppBamWriters.stream().mapToInt(x -> x.readCount()).sum();

        List<SuppBamReprocessor> reprocessTasks = suppBamWriters.stream()
                .map(x -> new SuppBamReprocessor(config, x.filename(), bamWriter, partitionDataStore, consensusReads))
                .collect(Collectors.toList());

        List<Callable> callableTasks = reprocessTasks.stream().collect(Collectors.toList());

        RD_LOGGER.info("reprocessing {} supplementary reads from {} BAMs", totalSupplementaries, callableTasks.size());

        if(!TaskExecutor.executeTasks(callableTasks, config.Threads))
        {
            RD_LOGGER.error("failed to reprocess supplementary BAMs");
            System.exit(1);
        }

        if(!config.KeepInterimBams)
        {
            try
            {
                for(SuppBamWriter suppWriter : suppBamWriters)
                {
                    Files.deleteIfExists(Paths.get(suppWriter.filename()));
                }
            }
            catch(IOException e)
            {
                RD_LOGGER.error("error deleting supplementary bams {}", e.toString());
            }
        }

        RD_LOGGER.debug("supplementary reprocessing complete");
    }
}
