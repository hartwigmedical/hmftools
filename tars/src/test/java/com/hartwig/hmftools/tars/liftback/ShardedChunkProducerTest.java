package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.SAM_DICTIONARY_V37;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import com.hartwig.hmftools.tars.liftback.ShardedChunkProducer.ShardRange;
import com.hartwig.hmftools.tars.liftback.ShardedChunkProducer.ShardRecordIterator;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;

public class ShardedChunkProducerTest
{
    @Rule
    public TemporaryFolder mFolder = new TemporaryFolder();

    private static SAMFileHeader nameGroupedHeader()
    {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        header.setSequenceDictionary(SAM_DICTIONARY_V37);
        return header;
    }

    private File writeNameGroupedBam(final SAMFileHeader header, final int fragments) throws IOException
    {
        File bam = mFolder.newFile("small.bam");
        try(SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, true, bam))
        {
            for(int i = 0; i < fragments; ++i)
            {
                String name = String.format("frag%04d", i);
                writer.addAlignment(createSamRecord(name, CHR_1, 100 + i, "ACGTACGTAC", "10M", CHR_1, 300 + i, false, false, null));
                writer.addAlignment(createSamRecord(name, CHR_1, 300 + i, "ACGTACGTAC", "10M", CHR_1, 100 + i, false, false, null));
            }
        }
        return bam;
    }

    // regression for the small-input EOF bug: a high shard count once produced EOF as a shard start.
    @Test
    public void testSmallInputHighShardCountNeverStartsAtEof() throws IOException
    {
        SAMFileHeader header = nameGroupedHeader();
        File bam = writeNameGroupedBam(header, 40);   // far fewer BGZF blocks than shardCount

        List<ShardRange> ranges = ShardedChunkProducer.computeSplits(bam, header, 64);

        assertFalse(ranges.isEmpty());
        for(ShardRange range : ranges)
            assertFalse("EOF used as a shard start", range.startVptr() == ShardedChunkProducer.EOF);

        List<String> names = new ArrayList<>();
        for(ShardRange range : ranges)
        {
            try(ShardRecordIterator iterator = new ShardRecordIterator(bam, header, range))
            {
                while(iterator.hasNext())
                    names.add(iterator.next().getReadName());
            }
        }

        assertEquals(80, names.size());   // 40 fragments x 2 reads
        for(int i = 0; i < names.size(); i += 2)
            assertEquals("fragment split across a shard boundary", names.get(i), names.get(i + 1));
    }

    @Test
    public void testSingleShardSpansWholeFile() throws IOException
    {
        SAMFileHeader header = nameGroupedHeader();
        File bam = writeNameGroupedBam(header, 10);

        List<ShardRange> ranges = ShardedChunkProducer.computeSplits(bam, header, 1);
        assertEquals(1, ranges.size());
        assertTrue(ranges.get(0).endVptr() == ShardedChunkProducer.EOF);
    }

    private static SAMRecord read(final String name)
    {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName(name);
        return record;
    }

    private static List<List<SAMRecord>> chunk(final List<SAMRecord> records, final int targetReads) throws InterruptedException
    {
        BlockingQueue<List<SAMRecord>> queue = new LinkedBlockingQueue<>();
        ShardedChunkProducer.streamChunks(records.iterator(), targetReads, queue, null);
        return new ArrayList<>(queue);
    }

    private static void assertWholeFragments(final List<List<SAMRecord>> chunks)
    {
        for(int i = 0; i + 1 < chunks.size(); ++i)
        {
            String lastName = chunks.get(i).get(chunks.get(i).size() - 1).getReadName();
            String nextName = chunks.get(i + 1).get(0).getReadName();
            assertTrue("fragment split across chunk boundary: " + lastName, !lastName.equals(nextName));
        }
    }

    @Test
    public void testCutsOnlyAtNameBoundary() throws InterruptedException
    {
        // three records per fragment against a target of 2 reads: every cut must wait for the name to change.
        List<SAMRecord> records = new ArrayList<>();
        for(int fragment = 0; fragment < 5; ++fragment)
        {
            records.add(read("frag" + fragment));
            records.add(read("frag" + fragment));
            records.add(read("frag" + fragment));
        }

        List<List<SAMRecord>> chunks = chunk(records, 2);

        assertWholeFragments(chunks);

        int total = 0;
        for(List<SAMRecord> chunk : chunks)
        {
            total += chunk.size();
        }
        assertEquals(15, total);
    }

    @Test
    public void testEmitsSingleChunkUnderTarget() throws InterruptedException
    {
        List<SAMRecord> records = List.of(read("a"), read("a"), read("b"));
        List<List<SAMRecord>> chunks = chunk(records, 5000);
        assertEquals(1, chunks.size());
        assertEquals(3, chunks.get(0).size());
    }

    @Test
    public void testOversizeFragmentStaysWhole() throws InterruptedException
    {
        List<SAMRecord> records = new ArrayList<>();
        for(int i = 0; i < 10; ++i)
        {
            records.add(read("big"));
        }
        records.add(read("next"));

        List<List<SAMRecord>> chunks = chunk(records, 3);
        assertWholeFragments(chunks);
        assertEquals(10, chunks.get(0).size());
    }

    @Test
    public void testProducerEmitsAllRecordsThenEndOfStream() throws Exception
    {
        SAMFileHeader header = nameGroupedHeader();
        File bam = writeNameGroupedBam(header, 50);
        File ref = mFolder.newFile("ref.fa");   // reading the BAM header does not touch the reference

        int workers = 2;
        BlockingQueue<List<SAMRecord>> queue = new LinkedBlockingQueue<>();
        ShardedChunkProducer producer = new ShardedChunkProducer(
                bam.getAbsolutePath(), ref.getAbsolutePath(), queue, workers, 8, 4);
        producer.start();

        int records = 0;
        int sentinels = 0;
        while(sentinels < workers)
        {
            List<SAMRecord> chunk = queue.take();
            if(chunk == ShardedChunkProducer.END_OF_STREAM)
                ++sentinels;
            else
                records += chunk.size();
        }
        producer.join();

        assertEquals(100, records);   // 50 fragments x 2 reads
        assertEquals(workers, sentinels);
    }
}
