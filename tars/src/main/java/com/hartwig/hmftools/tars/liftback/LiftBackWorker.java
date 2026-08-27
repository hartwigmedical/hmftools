package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;
import static com.hartwig.hmftools.tars.liftback.ShardedChunkProducer.END_OF_STREAM;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;

import com.hartwig.hmftools.tars.common.TarsCigarUtils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
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
    private final RegionPerfTracker mRegionPerf;

    private String mGroupContig;
    private int mGroupPosition;

    public LiftBackWorker(
            final BlockingQueue<List<SAMRecord>> queue, final LiftBackResources resources,
            final SAMFileHeader header, final String shardBam, final RegionPerfTracker regionPerf)
    {
        mQueue = queue;
        mProcessor = resources.createProcessor();
        mShardWriter = new SAMFileWriterFactory().makeBAMWriter(header, false, new File(shardBam));
        mRegionPerf = regionPerf;
    }

    public RegionPerfTracker regionPerf() { return mRegionPerf; }

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
                processGroup(group);
                group.clear();
            }
            group.add(record);
            currentName = name;
        }

        if(!group.isEmpty())
        {
            processGroup(group);
        }
    }

    private void processGroup(final List<SAMRecord> group)
    {
        if(mRegionPerf == null)
        {
            mProcessor.processNameGroup(group, this::write);
            return;
        }

        mGroupContig = null;
        mGroupPosition = 0;
        int readCount = group.size();

        long startTimeNanos = System.nanoTime();
        mProcessor.processNameGroup(group, this::write);
        mRegionPerf.add(mGroupContig, mGroupPosition, System.nanoTime() - startTimeNanos, readCount);
    }

    private void write(final SAMRecord record)
    {
        sanitizeForOutput(record);

        if(mRegionPerf != null && mGroupContig == null && !record.getReadUnmappedFlag())
        {
            mGroupContig = record.getReferenceName();
            mGroupPosition = record.getAlignmentStart();
        }

        mShardWriter.addAlignment(record);
    }

    // htsjdk and redux reject a zero-length CIGAR element or a SEQ length that disagrees with the CIGAR. Normalise the
    // CIGAR and, if SEQ still disagrees (e.g. a failed-lift supplementary mirrored onto its primary's coords), write a
    // matching all-M placeholder so the read stays valid and in the SA chain.
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
