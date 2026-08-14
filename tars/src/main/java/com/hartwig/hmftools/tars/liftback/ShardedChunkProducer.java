package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.LongAdder;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.BlockCompressedInputStream;

// Reads a name-sorted BAM in parallel. With no coordinate index to slice by, it shards by raw byte range: snap a
// target offset to the next BGZF block, advance to the first read-name change, decode each range on its own thread,
// then cut those records into name-aligned chunks on the shared queue for the workers.
// Every boundary is a read-name boundary, so a shard holds whole fragments: mates/supplementaries never split and
// there is no cross-shard reconciliation.
// BGZF block layout (section 4.1) and the virtual file offset shifted here (coffset<<16, section 4.1.1) are defined
// in the SAM/BAM spec: https://samtools.github.io/hts-specs/SAMv1.pdf
public class ShardedChunkProducer extends Thread
{
    private final String mInputBam;
    private final String mRefGenomeFile;
    private final BlockingQueue<List<SAMRecord>> mQueue;
    private final int mWorkerCount;
    private final int mChunkTargetReads;
    private final int mShardCount;

    private static final long PROGRESS_INTERVAL_MS = 15_000;

    // enqueued once per worker to signal end-of-stream; compared by reference.
    public static final List<SAMRecord> END_OF_STREAM = new ArrayList<>();

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
            File bam = new File(mInputBam);
            SAMFileHeader header = readHeader(bam);
            List<ShardRange> ranges = computeSplits(bam, header, mShardCount);
            TARS_LOGGER.info("liftback reading input across {} shard(s)", ranges.size());

            // kept open until all shards finish so the monitor can read their offsets.
            List<ShardRecordIterator> iterators = new ArrayList<>();
            for(ShardRange range : ranges)
            {
                iterators.add(new ShardRecordIterator(bam, header, range));
            }

            LongAdder readsCounter = new LongAdder();
            AtomicBoolean done = new AtomicBoolean(false);

            Thread monitor = new Thread(() -> runMonitor(iterators, bam.length(), readsCounter, done), "tars-progress");
            monitor.setDaemon(true);
            monitor.start();

            List<Thread> shardThreads = new ArrayList<>();
            for(ShardRecordIterator iter : iterators)
            {
                Thread shard = new Thread(() -> readShard(iter, readsCounter), "tars-shard");
                shard.start();
                shardThreads.add(shard);
            }

            for(Thread shard : shardThreads)
            {
                shard.join();
            }

            done.set(true);
            monitor.interrupt();
            closeQuietly(iterators);

            for(int i = 0; i < mWorkerCount; ++i)
            {
                mQueue.put(END_OF_STREAM);
            }
        }
        catch(Exception e)
        {
            TARS_LOGGER.error("liftback sharded producer failed: {}", e.toString());
            System.exit(1);
        }
    }

    private void readShard(final ShardRecordIterator iter, final LongAdder readsCounter)
    {
        try
        {
            streamChunks(iter, mChunkTargetReads, mQueue, readsCounter);
        }
        catch(Exception e)
        {
            TARS_LOGGER.error("liftback shard reader failed: {}", e.toString());
            System.exit(1);
        }
    }

    static final long EOF = Long.MAX_VALUE;

    // BGZF block-header magic: gzip magic (1f 8b), deflate method (08), FEXTRA flag (04)
    private static final byte[] BGZF_MAGIC = { 0x1f, (byte) 0x8b, 0x08, 0x04 };
    private static final int SCAN_BUFFER = 64 * 1024;

    // [startVptr, endVptr); endVptr == EOF means read to the end of the file.
    record ShardRange(long startVptr, long endVptr)
    {
    }

    static List<ShardRange> computeSplits(final File bam, final SAMFileHeader header, final int shardCount)
            throws IOException
    {
        long firstRecordVptr = headerEndVptr(bam, header);
        if(shardCount <= 1)
        {
            return List.of(new ShardRange(firstRecordVptr, EOF));
        }

        long compressedLength = bam.length();
        List<Long> splits = new ArrayList<>();
        splits.add(firstRecordVptr);

        long lastBlockOffset = firstRecordVptr >>> 16;
        for(int k = 1; k < shardCount; ++k)
        {
            long targetOffset = compressedLength * k / shardCount;
            if(targetOffset <= lastBlockOffset)
                continue; // shards collapsing onto the same block: fewer, larger shards is fine

            long blockStart = findBlockStart(bam, targetOffset);
            if(blockStart < 0)
                break; // no further block boundary: the remaining file is one shard

            long splitVptr = firstGroupBoundaryVptr(bam, header, blockStart << 16);
            if(splitVptr == EOF)
                break; // no group boundary before EOF: the rest is one shard, and EOF as a shard start would seek to Long.MAX_VALUE

            if(splitVptr > splits.get(splits.size() - 1))
            {
                splits.add(splitVptr);
                lastBlockOffset = splitVptr >>> 16;
            }
        }

        List<ShardRange> ranges = new ArrayList<>();
        for(int i = 0; i < splits.size(); ++i)
        {
            long end = i + 1 < splits.size() ? splits.get(i + 1) : EOF;
            ranges.add(new ShardRange(splits.get(i), end));
        }
        return ranges;
    }

    private static long headerEndVptr(final File bam, final SAMFileHeader header) throws IOException
    {
        try(BlockCompressedInputStream stream = new BlockCompressedInputStream(bam))
        {
            BinaryCodec codec = new BinaryCodec(stream);
            codec.readBytes(new byte[4]);             // "BAM\1"
            int textLength = codec.readInt();
            codec.readBytes(new byte[textLength]);    // header text
            int refCount = codec.readInt();
            for(int i = 0; i < refCount; ++i)
            {
                int nameLength = codec.readInt();
                codec.readBytes(new byte[nameLength]); // reference name
                codec.readInt();                      // reference length
            }
            return stream.getFilePointer();
        }
    }

    // next BGZF block header at/after targetOffset, as a compressed offset; -1 if none.
    private static long findBlockStart(final File bam, final long targetOffset) throws IOException
    {
        try(RandomAccessFile raf = new RandomAccessFile(bam, "r"))
        {
            long length = raf.length();
            byte[] buffer = new byte[SCAN_BUFFER];
            long offset = targetOffset;

            while(offset < length - 18)
            {
                raf.seek(offset);
                int read = raf.read(buffer);
                if(read <= 0)
                {
                    return -1;
                }

                for(int i = 0; i + 18 <= read; ++i)
                {
                    if(buffer[i] == BGZF_MAGIC[0] && buffer[i + 1] == BGZF_MAGIC[1]
                            && buffer[i + 2] == BGZF_MAGIC[2] && buffer[i + 3] == BGZF_MAGIC[3]
                            && buffer[i + 12] == 0x42 && buffer[i + 13] == 0x43) // 'BC' subfield id
                        return offset + i;
                }

                offset += read - 17; // overlap so a header straddling the buffer edge is not missed
            }
        }
        return -1;
    }

    // virtual pointer of the first read-name change at/after blockVptr; EOF if none.
    private static long firstGroupBoundaryVptr(final File bam, final SAMFileHeader header, final long blockVptr)
            throws IOException
    {
        try(BlockCompressedInputStream stream = new BlockCompressedInputStream(bam))
        {
            stream.seek(blockVptr);
            BAMRecordCodec codec = new BAMRecordCodec(header);
            codec.setInputStream(stream);

            String firstName = null;
            long recordVptr = stream.getFilePointer();
            SAMRecord record;
            while((record = codec.decode()) != null)
            {
                if(firstName == null)
                {
                    firstName = record.getReadName();
                }
                else if(!record.getReadName().equals(firstName))
                {
                    return recordVptr;
                }

                recordVptr = stream.getFilePointer();
            }
            return EOF;
        }
    }

    // Both range ends are read-name boundaries, so stopping once the next record reaches endVptr keeps fragments whole.
    static final class ShardRecordIterator implements Iterator<SAMRecord>, Closeable
    {
        private final BlockCompressedInputStream mStream;
        private final BAMRecordCodec mCodec;
        private final long mEndVptr;
        private final long mStartOffset;
        private SAMRecord mNext;

        ShardRecordIterator(final File bam, final SAMFileHeader header, final ShardRange range)
        {
            try
            {
                mStream = new BlockCompressedInputStream(bam);
                mStream.seek(range.startVptr());
            }
            catch(IOException e)
            {
                throw new UncheckedIOException(e);
            }

            mCodec = new BAMRecordCodec(header);
            mCodec.setInputStream(mStream);
            mEndVptr = range.endVptr();
            mStartOffset = range.startVptr() >>> 16;
            mNext = readNext();
        }

        // compressed bytes consumed so far; the monitor sums these across shards.
        long consumedBytes()
        {
            return Math.max(0, (mStream.getFilePointer() >>> 16) - mStartOffset);
        }

        private SAMRecord readNext()
        {
            if(mEndVptr != EOF && mStream.getFilePointer() >= mEndVptr)
            {
                return null;
            }

            return mCodec.decode(); // null at EOF
        }

        @Override
        public boolean hasNext()
        {
            return mNext != null;
        }

        @Override
        public SAMRecord next()
        {
            if(mNext == null)
            {
                throw new NoSuchElementException();
            }

            SAMRecord record = mNext;
            mNext = readNext();
            return record;
        }

        @Override
        public void close() throws IOException
        {
            mStream.close();
        }
    }

    // Flush only at a name boundary once the chunk reaches targetReads: splitting a fragment would orphan its
    // mates/supplementaries and defeat the per-group cache. readsCounter is nullable.
    static void streamChunks(
            final Iterator<SAMRecord> iter, final int targetReads,
            final BlockingQueue<List<SAMRecord>> queue, final LongAdder readsCounter)
            throws InterruptedException
    {
        List<SAMRecord> chunk = new ArrayList<>(targetReads);
        String currentName = null;

        while(iter.hasNext())
        {
            SAMRecord record = iter.next();
            String name = record.getReadName();

            if(currentName != null && !name.equals(currentName) && chunk.size() >= targetReads)
            {
                if(readsCounter != null)
                {
                    readsCounter.add(chunk.size());
                }
                queue.put(chunk);
                chunk = new ArrayList<>(targetReads);
            }

            chunk.add(record);
            currentName = name;
        }

        if(!chunk.isEmpty())
        {
            if(readsCounter != null)
            {
                readsCounter.add(chunk.size());
            }
            queue.put(chunk);
        }
    }

    // periodic progress: reads processed and rough % of the input consumed across all shards.
    private void runMonitor(
            final List<ShardRecordIterator> iterators, final long fileLength, final LongAdder readsCounter,
            final AtomicBoolean done)
    {
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
            {
                return;
            }

            long consumed = 0;
            for(ShardRecordIterator iter : iterators)
            {
                consumed += iter.consumedBytes();
            }

            int percent = fileLength > 0 ? (int) Math.min(100, consumed * 100 / fileLength) : 0;
            TARS_LOGGER.info("liftback processed {} reads ({}% of input)", readsCounter.sum(), percent);
        }
    }

    private static void closeQuietly(final List<ShardRecordIterator> iterators)
    {
        for(ShardRecordIterator iter : iterators)
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

    private SAMFileHeader readHeader(final File bam) throws IOException
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
