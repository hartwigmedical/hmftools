package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConfig.TARS_LOGGER;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.LongAdder;

import com.hartwig.hmftools.tars.liftback.shard.BamShardSplitter;
import com.hartwig.hmftools.tars.liftback.shard.ShardRecordIterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

// Parallel input reader: a single-thread BGZF parse starves the workers, so the input is split into byte
// ranges that each begin on a read-name-group boundary (BamShardSplitter) and read by one thread apiece. Each
// shard cuts whole-fragment chunks into the shared queue; the fragment invariant holds because the splits land
// between groups. A monitor thread reports a single aggregate progress bar (compressed bytes consumed across
// all shards / file size). After every shard drains, one END_OF_STREAM per worker is enqueued.
public class ShardedChunkProducer extends Thread
{
    private final String mInputBam;
    private final String mRefGenomeFile;
    private final BlockingQueue<List<SAMRecord>> mQueue;
    private final int mWorkerCount;
    private final int mChunkTargetReads;
    private final int mShardCount;

    private static final long PROGRESS_INTERVAL_MS = 15_000;
    private static final int PROGRESS_BAR_WIDTH = 24;

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

            // iterators stay open until all shards finish so the monitor can sample their offsets safely.
            final List<ShardRecordIterator> iterators = new ArrayList<>();
            for(final BamShardSplitter.ShardRange range : ranges)
                iterators.add(new ShardRecordIterator(bam, header, range));

            final LongAdder readsCounter = new LongAdder();
            final AtomicBoolean done = new AtomicBoolean(false);
            final AtomicBoolean failed = new AtomicBoolean(false);

            final Thread monitor = new Thread(() -> runMonitor(iterators, bam.length(), readsCounter, done), "tars-progress");
            monitor.setDaemon(true);
            monitor.start();

            final List<Thread> shardThreads = new ArrayList<>();
            for(final ShardRecordIterator iter : iterators)
            {
                final Thread shard = new Thread(() -> readShard(iter, readsCounter, failed), "tars-shard");
                shard.start();
                shardThreads.add(shard);
            }

            for(final Thread shard : shardThreads)
                shard.join();

            done.set(true);
            monitor.interrupt();
            closeQuietly(iterators);

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

    private void readShard(final ShardRecordIterator iter, final LongAdder readsCounter, final AtomicBoolean failed)
    {
        try
        {
            ChunkProducer.streamChunks(iter, mChunkTargetReads, mQueue::put, readsCounter);
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

    // one aggregate progress line per tick: a bar of compressed bytes consumed across all shards, plus the
    // read count and throughput since the last tick. Logged to a file in production, so each tick is its own
    // line (no in-place redraw).
    private void runMonitor(
            final List<ShardRecordIterator> iterators, final long fileLength, final LongAdder readsCounter,
            final AtomicBoolean done)
    {
        long lastReads = 0;
        long lastNanos = System.nanoTime();
        while(!done.get())
        {
            try
            {
                Thread.sleep(PROGRESS_INTERVAL_MS);
            }
            catch(InterruptedException e)
            {
                return;
            }
            if(done.get())
                return;

            long consumed = 0;
            for(final ShardRecordIterator iter : iterators)
                consumed += iter.consumedBytes();

            final int percent = fileLength > 0 ? (int) Math.min(100, consumed * 100 / fileLength) : 0;
            final long reads = readsCounter.sum();
            final long nowNanos = System.nanoTime();
            final double intervalSecs = (nowNanos - lastNanos) / 1e9;
            final long rate = intervalSecs > 0 ? (long) ((reads - lastReads) / intervalSecs) : 0;

            TARS_LOGGER.info("liftback {} {}% | {} reads | {}/s",
                    bar(percent), percent, formatCount(reads), formatCount(rate));

            lastReads = reads;
            lastNanos = nowNanos;
        }
    }

    private static String bar(final int percent)
    {
        final int filled = Math.min(PROGRESS_BAR_WIDTH, percent * PROGRESS_BAR_WIDTH / 100);
        final StringBuilder sb = new StringBuilder(PROGRESS_BAR_WIDTH + 2);
        sb.append('[');
        for(int i = 0; i < PROGRESS_BAR_WIDTH; ++i)
            sb.append(i < filled ? '#' : '.');
        sb.append(']');
        return sb.toString();
    }

    // compact human count: 412M, 851K, 1.2B.
    private static String formatCount(final long value)
    {
        if(value >= 1_000_000_000)
            return String.format("%.1fB", value / 1e9);
        if(value >= 1_000_000)
            return String.format("%dM", value / 1_000_000);
        if(value >= 1_000)
            return String.format("%dK", value / 1_000);
        return String.valueOf(value);
    }

    private static void closeQuietly(final List<ShardRecordIterator> iterators)
    {
        for(final ShardRecordIterator iter : iterators)
        {
            try
            {
                iter.close();
            }
            catch(Exception ignored)
            {
            }
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
