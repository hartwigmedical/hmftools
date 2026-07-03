package com.hartwig.hmftools.tars.liftback;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.atomic.LongAdder;

import htsjdk.samtools.SAMRecord;

// Cuts a name-grouped record stream into chunks of whole read-name groups, the unit handed to each worker.
// A chunk is flushed only at a name boundary, so a fragment's records are never split across chunks (which
// would defeat the per-group cache and orphan mates/supps). Used by ShardedChunkProducer's per-shard readers.
public class ChunkProducer
{
    private ChunkProducer() { }

    // reference-compared sentinel enqueued (one per worker) to signal end-of-stream.
    public static final List<SAMRecord> END_OF_STREAM = new ArrayList<>();

    // sink that may block (the bounded queue) and so propagates InterruptedException.
    interface ChunkSink
    {
        void accept(List<SAMRecord> chunk) throws InterruptedException;
    }

    // cut the stream into chunks of >= targetReads, never splitting a name group. readsCounter (nullable) is
    // bumped by each chunk's size as it is flushed, so a shared monitor can report aggregate throughput.
    static void streamChunks(
            final Iterator<SAMRecord> iter, final int targetReads, final ChunkSink sink, final LongAdder readsCounter)
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
                sink.accept(chunk);
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
            sink.accept(chunk);
        }
    }
}
