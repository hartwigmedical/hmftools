package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;
import static com.hartwig.hmftools.tars.liftback.ShardedChunkProducer.END_OF_STREAM;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;

import com.hartwig.hmftools.tars.common.TarsCigarUtils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;

// Drains chunks, splits each into contiguous name-groups for the processor and writes lifted records to its own
// unsorted shard. Owns its processor and ref genome handle so nothing serialises across workers.
public class LiftBackWorker extends Thread
{
    private final BlockingQueue<List<SAMRecord>> mQueue;
    private final LiftBackGroupProcessor mProcessor;
    private final SAMFileWriter mShardWriter;

    public LiftBackWorker(
            final BlockingQueue<List<SAMRecord>> queue, final LiftBackResources resources,
            final SAMFileHeader header, final String shardBam)
    {
        mQueue = queue;
        mProcessor = resources.createProcessor();
        mShardWriter = new SAMFileWriterFactory().makeBAMWriter(header, false, new File(shardBam));
    }

    // counters are incremented on the pre-lift record so they count inputs not emitted records; read by TarsApplication
    // only after the worker threads join
    public LiftBackStats stats() { return mProcessor.stats(); }

    @Override
    public void run()
    {
        try
        {
            while(true)
            {
                List<SAMRecord> chunk = mQueue.take();
                if(chunk == END_OF_STREAM)
                    break;
                processChunk(chunk);
            }
        }
        catch(Exception e)
        {
            TARS_LOGGER.error("liftback worker failed: {}", e.toString(), e);
            System.exit(1);
        }
        finally
        {
            mShardWriter.close();
        }
    }

    private void processChunk(final List<SAMRecord> chunk)
    {
        List<SAMRecord> group = new ArrayList<>();
        String currentName = null;

        for(SAMRecord record : chunk)
        {
            String name = record.getReadName();
            if(currentName != null && !name.equals(currentName))
            {
                mProcessor.processNameGroup(group, this::write);
                group.clear();
            }
            group.add(record);
            currentName = name;
        }

        if(!group.isEmpty())
        {
            mProcessor.processNameGroup(group, this::write);
        }
    }

    private void write(final SAMRecord record)
    {
        sanitizeForOutput(record);
        mShardWriter.addAlignment(record);
    }

    // htsjdk and REDUX reject zero-length CIGAR elements or a SEQ length that disagrees with the CIGAR. Normalise safe
    // structural noise, but fail on a length mismatch rather than inventing an alignment that the input did not support.
    static void sanitizeForOutput(final SAMRecord record)
    {
        if(record.getReadUnmappedFlag())
        {
            return;
        }

        Cigar cigar = record.getCigar();
        if(cigar == null || cigar.isEmpty())
        {
            return;
        }

        Cigar cleaned = TarsCigarUtils.normalize(cigar);

        byte[] bases = record.getReadBases();
        int seqLength = bases != null ? bases.length : 0;

        if(seqLength > 0 && cleaned.getReadLength() != seqLength)
        {
            throw new IllegalStateException(String.format(
                    "read(%s) SEQ length(%d) disagrees with CIGAR(%s) read length(%d)",
                    record.getReadName(), seqLength, cleaned, cleaned.getReadLength()));
        }

        if(!cleaned.equals(cigar))
        {
            record.setCigar(cleaned);
        }
    }
}
