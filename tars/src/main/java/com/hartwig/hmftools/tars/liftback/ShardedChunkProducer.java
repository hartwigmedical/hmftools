package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConfig.TARS_LOGGER;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;

import com.hartwig.hmftools.tars.liftback.shard.BamShardSplitter;
import com.hartwig.hmftools.tars.liftback.shard.ShardRecordIterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

// Parallel replacement for ChunkProducer: the single-thread BGZF parse is the producer-side bottleneck (the
// workers sit idle waiting), so the input is split into byte ranges that each begin on a read-name-group
// boundary (BamShardSplitter) and read by one thread apiece. Each shard cuts whole-fragment chunks into the
// shared queue exactly as ChunkProducer does; the fragment invariant holds because the splits land between
// groups. After every shard drains, one END_OF_STREAM per worker is enqueued.
public class ShardedChunkProducer extends Thread
{
    private final String mInputBam;
    private final String mRefGenomeFile;
    private final BlockingQueue<List<SAMRecord>> mQueue;
    private final int mWorkerCount;
    private final int mChunkTargetReads;
    private final int mShardCount;

    public ShardedChunkProducer(
            final String inputBam, final String refGenomeFile, final BlockingQueue<List<SAMRecord>> queue,
            final int workerCount, final int chunkTargetReads, final int shardCount)
    {
        mInputBam = inputBam;
        mRefGenomeFile = refGenomeFile;
        mQueue = queue;
        mWorkerCount = workerCount;
        mChunkTargetReads = chunkTargetReads;
        mShardCount = Math.max(1, shardCount);
    }

    @Override
    public void run()
    {
        try
        {
            final File bam = new File(mInputBam);
            final SAMFileHeader header = readHeader(bam);
            final List<BamShardSplitter.ShardRange> ranges = BamShardSplitter.computeSplits(bam, header, mShardCount);
            TARS_LOGGER.info("liftback reading input across {} shard(s)", ranges.size());

            final AtomicBoolean failed = new AtomicBoolean(false);
            final List<Thread> shardThreads = new ArrayList<>();
            for(int i = 0; i < ranges.size(); ++i)
            {
                final BamShardSplitter.ShardRange range = ranges.get(i);
                final boolean logProgress = i == 0; // one shard logs progress; the rest stay quiet
                final Thread shard = new Thread(() -> readShard(bam, header, range, logProgress, failed), "tars-shard");
                shard.start();
                shardThreads.add(shard);
            }

            for(final Thread shard : shardThreads)
                shard.join();

            if(failed.get())
                throw new IllegalStateException("a shard reader failed");

            for(int i = 0; i < mWorkerCount; ++i)
                mQueue.put(ChunkProducer.END_OF_STREAM);
        }
        catch(InterruptedException e)
        {
            Thread.currentThread().interrupt();
            TARS_LOGGER.error("liftback sharded producer interrupted");
            System.exit(1);
        }
        catch(Exception e)
        {
            TARS_LOGGER.error("liftback sharded producer failed: {}", e.toString());
            System.exit(1);
        }
    }

    private void readShard(
            final File bam, final SAMFileHeader header, final BamShardSplitter.ShardRange range,
            final boolean logProgress, final AtomicBoolean failed)
    {
        try(ShardRecordIterator iter = new ShardRecordIterator(bam, header, range))
        {
            ChunkProducer.streamChunks(iter, mChunkTargetReads, mQueue::put, mQueue::size, logProgress);
        }
        catch(InterruptedException e)
        {
            Thread.currentThread().interrupt();
            failed.set(true);
        }
        catch(Exception e)
        {
            TARS_LOGGER.error("liftback shard reader failed: {}", e.toString());
            failed.set(true);
        }
    }

    private SAMFileHeader readHeader(final File bam) throws java.io.IOException
    {
        try(SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .referenceSequence(new File(mRefGenomeFile))
                .open(bam))
        {
            return reader.getFileHeader();
        }
    }
}
