package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;
import static com.hartwig.hmftools.tars.liftback.ShardedChunkProducer.END_OF_STREAM;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;

// Drains read-name chunks, transforms each contiguous name-group with a per-group LiftedMatePair, and writes
// lifted records to its own unsorted shard. Owns its engines and ref genome handle so nothing serialises across
// workers. Exits on END_OF_STREAM.
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

    // counted by the processor, which sees each record's pre-lift state. Read by TarsApplication after the join.
    public int recordsSeen() { return mProcessor.recordsSeen(); }

    public int primariesUnmapped() { return mProcessor.primariesUnmapped(); }

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

    // Every kept record lands here: normalise it, then write it to this worker's shard.
    private void write(final SAMRecord record)
    {
        sanitizeForOutput(record);
        mShardWriter.addAlignment(record);
    }

    // No malformed record leaves tars: htsjdk and redux reject a zero-length CIGAR element or a SEQ length that
    // disagrees with the CIGAR. Normalise the CIGAR and, if SEQ still disagrees - eg a failed-lift supplementary
    // mirrored onto its primary's coords - write a matching all-M placeholder so the read stays valid and in the SA
    // chain rather than being dropped or crashing downstream.
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

        Cigar cleaned = new Cigar(
                ContigTranslator.mergeAdjacentSameOp(ContigTranslator.dropZeroLength(cigar.getCigarElements())));

        byte[] bases = record.getReadBases();
        int seqLength = bases != null ? bases.length : 0;

        if(seqLength > 0 && cleaned.getReadLength() != seqLength)
        {
            TARS_LOGGER.warn("read({}) seq length({}) disagrees with cigar({}); writing {}M placeholder",
                    record.getReadName(), seqLength, cleaned, seqLength);
            cleaned = new Cigar(List.of(new CigarElement(seqLength, CigarOperator.M)));
        }

        if(!cleaned.equals(cigar))
        {
            record.setCigar(cleaned);
        }
    }
}
