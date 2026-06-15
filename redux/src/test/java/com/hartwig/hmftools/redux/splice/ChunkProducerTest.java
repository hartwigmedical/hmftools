package com.hartwig.hmftools.redux.splice;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;

// the producer's one safety-critical invariant: a chunk is only cut at a read-name boundary, so no
// fragment's records are ever split across two chunks (a split fragment would defeat the per-group cache
// and orphan mates/supps for dedup).
public class ChunkProducerTest
{
    private static SAMRecord read(final String name)
    {
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName(name);
        return record;
    }

    private static List<List<SAMRecord>> chunk(final List<SAMRecord> records, final int targetReads)
            throws InterruptedException
    {
        final List<List<SAMRecord>> chunks = new ArrayList<>();
        ChunkProducer.streamChunks(records.iterator(), targetReads, chunks::add);
        return chunks;
    }

    private static void assertWholeFragments(final List<List<SAMRecord>> chunks)
    {
        for(int i = 0; i + 1 < chunks.size(); ++i)
        {
            final String lastName = chunks.get(i).get(chunks.get(i).size() - 1).getReadName();
            final String nextName = chunks.get(i + 1).get(0).getReadName();
            assertTrue("fragment split across chunk boundary: " + lastName, !lastName.equals(nextName));
        }
    }

    @Test
    public void cutsOnlyAtNameBoundary()
    {
        // three records per fragment, target 2 reads: the cut must wait for the name to change rather than
        // firing mid-fragment, so each chunk holds whole fragments.
        final List<SAMRecord> records = new ArrayList<>();
        for(int f = 0; f < 5; ++f)
        {
            records.add(read("frag" + f));
            records.add(read("frag" + f));
            records.add(read("frag" + f));
        }

        final List<List<SAMRecord>> chunks;
        try
        {
            chunks = chunk(records, 2);
        }
        catch(InterruptedException e)
        {
            throw new AssertionError(e);
        }

        assertWholeFragments(chunks);

        int total = 0;
        for(final List<SAMRecord> c : chunks)
            total += c.size();
        assertEquals(15, total);
    }

    @Test
    public void emitsSingleChunkUnderTarget() throws InterruptedException
    {
        final List<SAMRecord> records = List.of(read("a"), read("a"), read("b"));
        final List<List<SAMRecord>> chunks = chunk(records, 5000);
        assertEquals(1, chunks.size());
        assertEquals(3, chunks.get(0).size());
    }

    @Test
    public void oversizeFragmentStaysWhole() throws InterruptedException
    {
        // a single fragment larger than the target must not be split -- the boundary check gates the cut.
        final List<SAMRecord> records = new ArrayList<>();
        for(int i = 0; i < 10; ++i)
            records.add(read("big"));
        records.add(read("next"));

        final List<List<SAMRecord>> chunks = chunk(records, 3);
        assertWholeFragments(chunks);
        assertEquals(10, chunks.get(0).size());
    }
}
