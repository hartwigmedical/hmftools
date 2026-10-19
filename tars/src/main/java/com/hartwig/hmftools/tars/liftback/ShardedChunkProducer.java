package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import static htsjdk.samtools.util.BlockCompressedStreamConstants.BGZF_ID1;
import static htsjdk.samtools.util.BlockCompressedStreamConstants.BGZF_ID2;
import static htsjdk.samtools.util.BlockCompressedStreamConstants.BGZF_LEN;
import static htsjdk.samtools.util.BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH;
import static htsjdk.samtools.util.BlockCompressedStreamConstants.GZIP_CM_DEFLATE;
import static htsjdk.samtools.util.BlockCompressedStreamConstants.GZIP_FLG;
import static htsjdk.samtools.util.BlockCompressedStreamConstants.GZIP_ID1;
import static htsjdk.samtools.util.BlockCompressedStreamConstants.GZIP_ID2;
import static htsjdk.samtools.util.BlockCompressedStreamConstants.GZIP_XLEN;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.io.UncheckedIOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.LongAdder;

import com.hartwig.hmftools.common.utils.MemoryCalcs;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.BlockCompressedInputStream;

// Shards a name-sorted BAM by raw byte range, there being no coordinate index to slice by. Every shard boundary is a
// read-name boundary, so mates and supplementaries never split across shards and no reconciliation is needed.
// BGZF block layout and the virtual offset shifted here (coffset<<16): https://samtools.github.io/hts-specs/SAMv1.pdf
public class ShardedChunkProducer extends Thread
{
    private final List<String> mInputBams;
    private final String mRefGenomeFile;
    private final BlockingQueue<List<SAMRecord>> mQueue;
    private final int mWorkerCount;
    private final int mChunkTargetReads;
    private final int mShardCount;

    private static final long PROGRESS_INTERVAL_MS = 15_000;

    private final AtomicInteger mPeakMemoryMb = new AtomicInteger();

    // compared by reference
    public static final List<SAMRecord> END_OF_STREAM = new ArrayList<>();

    public ShardedChunkProducer(
            final List<String> inputBams, final String refGenomeFile, final BlockingQueue<List<SAMRecord>> queue,
            final int workerCount, final int chunkTargetReads, final int shardCount)
    {
        mInputBams = inputBams;
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
            LongAdder readsCounter = new LongAdder();
            long totalBytes = totalInputBytes();
            long completedBytes = 0;

            for(String inputBam : mInputBams)
            {
                File bam = new File(inputBam);
                readInput(bam, readsCounter, completedBytes, totalBytes);
                completedBytes += bam.length();
            }

            mPeakMemoryMb.accumulateAndGet(MemoryCalcs.calcMemoryUsage(), Math::max);

            TARS_LOGGER.info("liftback read {} reads from {} input BAM(s), peak memory({}mb)",
                    readsCounter.sum(), mInputBams.size(), mPeakMemoryMb.get());

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

    private long totalInputBytes()
    {
        long total = 0;
        for(String inputBam : mInputBams)
        {
            total += new File(inputBam).length();
        }
        return total;
    }

    private void readInput(
            final File bam, final LongAdder readsCounter, final long completedBytes, final long totalBytes)
            throws IOException, InterruptedException
    {
        SAMFileHeader header = readHeader(bam);
        List<ShardRange> ranges = computeSplits(bam, header, mShardCount);
        TARS_LOGGER.info("liftback reading {} across {} shard(s)", bam.getName(), ranges.size());

        // kept open until all shards finish so the monitor can read their offsets.
        List<ShardRecordIterator> iterators = new ArrayList<>();
        for(ShardRange range : ranges)
        {
            iterators.add(new ShardRecordIterator(bam, header, range));
        }

        AtomicBoolean done = new AtomicBoolean(false);

        Thread monitor = new Thread(
                () -> runMonitor(iterators, completedBytes, totalBytes, readsCounter, done), "tars-progress");
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

    private static final int BLOCK_SCAN_LENGTH = 64 * 1024;
    private static final int RECORD_SCAN_LENGTH = 4 * 64 * 1024;
    private static final int VALIDATED_RECORD_COUNT = 4;

    // [startVptr, endVptr); endVptr == EOF means read to the end of the file.
    record ShardRange(long startVptr, long endVptr, String endReadName)
    {
    }

    private record GroupBoundary(long vptr, String readName)
    {
    }

    static List<ShardRange> computeSplits(final File bam, final SAMFileHeader header, final int shardCount)
            throws IOException
    {
        long firstRecordVptr = headerEndVptr(bam, header);
        if(shardCount <= 1)
        {
            return List.of(new ShardRange(firstRecordVptr, EOF, null));
        }

        long compressedLength = bam.length();
        List<GroupBoundary> boundaries = new ArrayList<>();
        long lastBlockOffset = firstRecordVptr >>> 16;
        for(int k = 1; k < shardCount; ++k)
        {
            long targetOffset = compressedLength * k / shardCount;
            if(targetOffset <= lastBlockOffset)
            {
                continue;
            }

            long blockStart = findBlockStart(bam, targetOffset);
            if(blockStart < 0)
            {
                break;
            }

            GroupBoundary boundary = firstGroupBoundary(bam, header, blockStart);
            if(boundary == null)
            {
                break;
            }

            boundaries.add(boundary);
            lastBlockOffset = boundary.vptr() >>> 16;
        }

        List<ShardRange> ranges = new ArrayList<>();
        long startVptr = firstRecordVptr;
        for(GroupBoundary boundary : boundaries)
        {
            ranges.add(new ShardRange(startVptr, boundary.vptr(), boundary.readName()));
            startVptr = boundary.vptr();
        }
        ranges.add(new ShardRange(startVptr, EOF, null));
        return ranges;
    }

    private static long findBlockStart(final File bam, final long targetOffset) throws IOException
    {
        try(RandomAccessFile file = new RandomAccessFile(bam, "r"))
        {
            byte[] buffer = new byte[BLOCK_SCAN_LENGTH];
            long offset = targetOffset;
            while(offset <= file.length() - BLOCK_HEADER_LENGTH)
            {
                file.seek(offset);
                int read = file.read(buffer);
                for(int i = 0; i + BLOCK_HEADER_LENGTH <= read; ++i)
                {
                    if(isBlockHeader(buffer, i))
                    {
                        return offset + i;
                    }
                }

                offset += read - BLOCK_HEADER_LENGTH + 1;
            }
        }
        return -1;
    }

    private static boolean isBlockHeader(final byte[] buffer, final int offset)
    {
        return buffer[offset] == GZIP_ID1 && (buffer[offset + 1] & 0xFF) == GZIP_ID2
                && buffer[offset + 2] == GZIP_CM_DEFLATE && buffer[offset + 3] == GZIP_FLG
                && buffer[offset + 10] == GZIP_XLEN && buffer[offset + 11] == 0
                && buffer[offset + 12] == BGZF_ID1 && buffer[offset + 13] == BGZF_ID2
                && buffer[offset + 14] == BGZF_LEN && buffer[offset + 15] == 0;
    }

    private static GroupBoundary firstGroupBoundary(final File bam, final SAMFileHeader header, final long blockStart)
            throws IOException
    {
        try(BlockCompressedInputStream stream = new BlockCompressedInputStream(bam))
        {
            stream.seek(blockStart << 16);
            byte[] window = stream.readNBytes(RECORD_SCAN_LENGTH);
            int recordOffset = firstRecordOffset(
                    window, window.length < RECORD_SCAN_LENGTH, header.getSequenceDictionary().size());
            if(recordOffset < 0)
            {
                return null;
            }

            stream.seek(blockStart << 16);
            stream.skipNBytes(recordOffset);
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
                    return new GroupBoundary(recordVptr, record.getReadName());
                }

                recordVptr = stream.getFilePointer();
            }
            return null;
        }
    }

    private static int firstRecordOffset(final byte[] window, final boolean endOfFile, final int referenceCount)
    {
        ByteBuffer buffer = ByteBuffer.wrap(window).order(ByteOrder.LITTLE_ENDIAN);
        for(int offset = 0; offset < window.length; ++offset)
        {
            if(startsRecordChain(buffer, offset, endOfFile, referenceCount))
            {
                return offset;
            }
        }
        return -1;
    }

    private static boolean startsRecordChain(
            final ByteBuffer buffer, final int offset, final boolean endOfFile, final int referenceCount)
    {
        int recordStart = offset;
        for(int i = 0; i < VALIDATED_RECORD_COUNT; ++i)
        {
            if(recordStart == buffer.limit())
            {
                return endOfFile && i > 0;
            }

            recordStart = recordEnd(buffer, recordStart, referenceCount);
            if(recordStart < 0)
            {
                return false;
            }
        }
        return true;
    }

    private static int recordEnd(final ByteBuffer buffer, final int offset, final int referenceCount)
    {
        int nameStart = offset + 36;
        if(nameStart > buffer.limit())
        {
            return -1;
        }

        int blockSize = buffer.getInt(offset);
        int referenceIndex = buffer.getInt(offset + 4);
        int readNameLength = buffer.get(offset + 12) & 0xFF;
        int mateReferenceIndex = buffer.getInt(offset + 24);
        long recordEnd = offset + 4L + blockSize;
        if(referenceIndex < -1 || referenceIndex >= referenceCount || mateReferenceIndex < -1 || mateReferenceIndex >= referenceCount
                || readNameLength < 2 || blockSize < 32 + readNameLength || recordEnd > buffer.limit())
        {
            return -1;
        }

        int nameEnd = nameStart + readNameLength - 1;
        for(int i = nameStart; i < nameEnd; ++i)
        {
            if(buffer.get(i) < '!' || buffer.get(i) > '~')
            {
                return -1;
            }
        }
        return buffer.get(nameEnd) == 0 ? (int) recordEnd : -1;
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

    // Both range ends are read-name boundaries, so stopping once the next record reaches endVptr keeps fragments whole.
    static final class ShardRecordIterator implements Iterator<SAMRecord>, Closeable
    {
        private final BlockCompressedInputStream mStream;
        private final BAMRecordCodec mCodec;
        private final long mEndVptr;
        private final String mEndReadName;
        private final long mStartOffset;
        private String mLastReadName;
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
            mEndReadName = range.endReadName();
            mStartOffset = range.startVptr() >>> 16;
            mNext = readNext();
        }

        long consumedBytes()
        {
            return Math.max(0, (mStream.getFilePointer() >>> 16) - mStartOffset);
        }

        private SAMRecord readNext()
        {
            long filePointer = mStream.getFilePointer();
            if(mEndVptr != EOF && filePointer >= mEndVptr)
            {
                if(filePointer != mEndVptr || mEndReadName.equals(mLastReadName))
                {
                    throw new IllegalStateException(String.format(
                            "shard end(%d) is not a read-name boundary: reached(%d) lastRead(%s) nextRead(%s)",
                            mEndVptr, filePointer, mLastReadName, mEndReadName));
                }
                return null;
            }

            SAMRecord record = mCodec.decode();
            if(record != null)
            {
                mLastReadName = record.getReadName();
            }
            return record;
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

    private void runMonitor(
            final List<ShardRecordIterator> iterators, final long completedBytes, final long totalBytes,
            final LongAdder readsCounter, final AtomicBoolean done)
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

            long consumed = completedBytes;
            for(ShardRecordIterator iter : iterators)
            {
                consumed += iter.consumedBytes();
            }

            int percent = totalBytes > 0 ? (int) Math.min(100, consumed * 100 / totalBytes) : 0;
            int memoryMb = MemoryCalcs.calcMemoryUsage();
            mPeakMemoryMb.accumulateAndGet(memoryMb, Math::max);

            TARS_LOGGER.info("liftback processed {} reads ({}% of input) memory({}mb)", readsCounter.sum(), percent, memoryMb);
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
