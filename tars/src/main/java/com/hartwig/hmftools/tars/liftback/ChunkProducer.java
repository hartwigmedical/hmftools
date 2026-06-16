package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bamops.BamToolName.fromPath;
import static com.hartwig.hmftools.tars.common.TarsConfig.TARS_LOGGER;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.BlockingQueue;

import com.hartwig.hmftools.common.bamops.BamToolName;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

// Single producer over bwa's name-grouped output (mates + supplementaries already contiguous, FASTQ order --
// no sort or index needed). Streams the BAM once, batching contiguous read-name groups into chunks of
// ~CHUNK_TARGET_READS reads, cut only at a name boundary so every chunk holds whole fragments. Each chunk goes
// to a worker that lifts and emits it -- no fragment cache, no re-fetch. One END_OF_STREAM per worker at the end.
//
// To stop a single-threaded BGZF decompress from starving the workers on a whole-sample BAM, the input is piped
// through the bam tool's MULTI-THREADED decompression (uncompressed BAM out) when a bamtool is configured, so
// htsjdk parses with no decompression on this thread. Falls back to a direct htsjdk read otherwise.
public class ChunkProducer extends Thread
{
    // reference-compared sentinel the producer enqueues (one per worker) to signal end-of-stream.
    public static final List<SAMRecord> END_OF_STREAM = new ArrayList<>();

    private final String mInputBam;
    private final String mRefGenomeFile;
    private final BlockingQueue<List<SAMRecord>> mQueue;
    private final int mWorkerCount;
    private final int mChunkTargetReads;
    private final String mBamToolPath; // nullable: when set, decompress via the tool's threads
    private final int mDecompressThreads;

    // buffer over the decompress pipe: `-u` emits high-volume uncompressed BGZF and htsjdk reads blocks straight
    // off the raw pipe, so a large buffer batches the otherwise tiny per-block read() syscalls.
    private static final int PIPE_BUFFER_BYTES = 1 << 20;

    public ChunkProducer(
            final String inputBam, final String refGenomeFile, final BlockingQueue<List<SAMRecord>> queue,
            final int workerCount, final int chunkTargetReads, final String bamToolPath, final int decompressThreads)
    {
        mInputBam = inputBam;
        mRefGenomeFile = refGenomeFile;
        mQueue = queue;
        mWorkerCount = workerCount;
        mChunkTargetReads = chunkTargetReads;
        mBamToolPath = bamToolPath;
        mDecompressThreads = decompressThreads;
    }

    @Override
    public void run()
    {
        try
        {
            if(mBamToolPath != null)
                streamViaToolDecompress();
            else
                streamDirect();

            for(int i = 0; i < mWorkerCount; ++i)
                mQueue.put(END_OF_STREAM);
        }
        catch(InterruptedException e)
        {
            Thread.currentThread().interrupt();
            TARS_LOGGER.error("liftback chunk producer interrupted");
            System.exit(1);
        }
        catch(Exception e)
        {
            TARS_LOGGER.error("liftback chunk producer failed: {}", e.toString());
            System.exit(1);
        }
    }

    private void streamDirect() throws InterruptedException, IOException
    {
        try(SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .referenceSequence(new File(mRefGenomeFile))
                .open(new File(mInputBam)))
        {
            streamChunks(reader.iterator());
        }
    }

    // pipe the name-grouped BAM through `samtools view -u -@ N` (or `sambamba view -f bam -l 0 -t N`): the tool
    // decompresses across N threads and emits uncompressed BAM, so htsjdk parses it with no decompression here.
    private void streamViaToolDecompress() throws InterruptedException, IOException
    {
        final List<String> command = decompressCommand();
        TARS_LOGGER.info("liftback streaming input via {} ({} decompress threads)", fromPath(mBamToolPath), mDecompressThreads);

        final Process process = new ProcessBuilder(command).redirectError(ProcessBuilder.Redirect.INHERIT).start();
        try(SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT)
                .open(SamInputResource.of(new BufferedInputStream(process.getInputStream(), PIPE_BUFFER_BYTES))))
        {
            streamChunks(reader.iterator());
        }
        catch(Exception e)
        {
            // only kill the tool on a real failure; do NOT destroy() on the happy path -- on a whole-sample BAM the
            // tool is still flushing when we finish reading, and a SIGTERM there surfaces as a spurious exit 143.
            process.destroyForcibly();
            throw e;
        }

        final int exitCode = process.waitFor();
        if(exitCode != 0)
            throw new IOException("decompress process exit code " + exitCode);
    }

    private List<String> decompressCommand()
    {
        final List<String> command = new ArrayList<>();
        command.add(mBamToolPath);
        command.add("view");
        if(fromPath(mBamToolPath) == BamToolName.SAMTOOLS)
        {
            command.add("-u"); // uncompressed BAM out
            command.add("-@");
            command.add(String.valueOf(mDecompressThreads));
        }
        else
        {
            command.add("-f");
            command.add("bam");
            command.add("-l");
            command.add("0"); // uncompressed
            command.add("-t");
            command.add(String.valueOf(mDecompressThreads));
        }
        command.add(mInputBam);
        return command;
    }

    private void streamChunks(final SAMRecordIterator iter) throws InterruptedException
    {
        streamChunks(iter, mChunkTargetReads, mQueue::put);
    }

    // sink that may block (the bounded queue) and so propagates InterruptedException.
    interface ChunkSink
    {
        void accept(List<SAMRecord> chunk) throws InterruptedException;
    }

    // progress is logged every this many reads streamed (whole-sample runs are otherwise silent for many minutes).
    private static final long PROGRESS_LOG_INTERVAL = 10_000_000;

    // cut the name-grouped record stream into chunks of >= targetReads, never splitting a read-name group:
    // a chunk is only flushed at a name boundary, so every chunk holds whole fragments.
    static void streamChunks(final Iterator<SAMRecord> iter, final int targetReads, final ChunkSink sink)
            throws InterruptedException
    {
        List<SAMRecord> chunk = new ArrayList<>(targetReads);
        String currentName = null;
        long readsStreamed = 0;
        int chunksQueued = 0;
        long nextProgressLog = PROGRESS_LOG_INTERVAL;

        while(iter.hasNext())
        {
            final SAMRecord record = iter.next();
            final String name = record.getReadName();

            if(currentName != null && !name.equals(currentName) && chunk.size() >= targetReads)
            {
                sink.accept(chunk);
                chunk = new ArrayList<>(targetReads);
                ++chunksQueued;
            }

            chunk.add(record);
            currentName = name;

            if(++readsStreamed == nextProgressLog)
            {
                TARS_LOGGER.info("liftback streamed {} reads ({} chunks queued)", readsStreamed, chunksQueued);
                nextProgressLog += PROGRESS_LOG_INTERVAL;
            }
        }

        if(!chunk.isEmpty())
        {
            sink.accept(chunk);
            ++chunksQueued;
        }

        TARS_LOGGER.info("liftback stream complete: {} reads in {} chunks", readsStreamed, chunksQueued);
    }
}
