package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.nameGroupedHeader;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.pairedRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.unpairedRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.writeBam;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.Collectors;

import com.hartwig.hmftools.tars.liftback.ShardedChunkProducer.ShardRange;
import com.hartwig.hmftools.tars.liftback.ShardedChunkProducer.ShardRecordIterator;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

public class ShardedChunkProducerTest
{
    @Rule
    public TemporaryFolder mFolder = new TemporaryFolder();

    private static SAMFileHeader header()
    {
        return nameGroupedHeader(new SAMSequenceRecord(CHR_1, 100_000));
    }

    private File writeNameGroupedBam(final SAMFileHeader header, final int fragments) throws IOException
    {
        return writeNameGroupedBam(header, fragments, "small.bam");
    }

    private File writeNameGroupedBam(final SAMFileHeader header, final int fragments, final String filename)
            throws IOException
    {
        File bam = mFolder.newFile(filename);
        writeBam(bam, header, records(fragments));
        return bam;
    }

    private static List<SAMRecord> records(final int fragments)
    {
        List<SAMRecord> records = new ArrayList<>();
        for(int i = 0; i < fragments; ++i)
        {
            String name = String.format("frag%04d", i);
            SAMRecord first = pairedRecord(name, true, CHR_1, 100 + i, "10M", false);
            SAMRecord second = pairedRecord(name, false, CHR_1, 300 + i, "10M", true);
            first.setMateReferenceName(CHR_1);
            first.setMateAlignmentStart(300 + i);
            first.setMateNegativeStrandFlag(true);
            second.setMateReferenceName(CHR_1);
            second.setMateAlignmentStart(100 + i);
            for(SAMRecord record : List.of(first, second))
            {
                record.setReadString("ACGTACGTAC");
                record.setBaseQualityString("IIIIIIIIII");
                records.add(record);
            }
        }
        return records;
    }

    // regression for the small-input EOF bug: a high shard count once produced EOF as a shard start.
    @Test
    public void testSmallInputHighShardCountNeverStartsAtEof() throws IOException
    {
        SAMFileHeader header = header();
        File bam = writeNameGroupedBam(header, 40);   // far fewer BGZF blocks than shardCount

        List<ShardRange> ranges = ShardedChunkProducer.computeSplits(bam, header, 64);

        assertFalse(ranges.isEmpty());
        for(ShardRange range : ranges)
            assertFalse("EOF used as a shard start", range.startVptr() == ShardedChunkProducer.EOF);

        List<SAMRecord> actual = new ArrayList<>();
        for(ShardRange range : ranges)
        {
            try(ShardRecordIterator iterator = new ShardRecordIterator(bam, header, range))
            {
                while(iterator.hasNext())
                    actual.add(iterator.next());
            }
        }

        assertEquals(samStrings(records(40)), samStrings(actual));
    }

    @Test
    public void testSingleShardSpansWholeFile() throws IOException
    {
        SAMFileHeader header = header();
        File bam = writeNameGroupedBam(header, 10);

        List<ShardRange> ranges = ShardedChunkProducer.computeSplits(bam, header, 1);
        assertEquals(1, ranges.size());
        assertTrue(ranges.get(0).endVptr() == ShardedChunkProducer.EOF);
    }

    private static SAMRecord read(final String name)
    {
        return unpairedRecord(name);
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
        assertEquals(records, chunks.stream().flatMap(List::stream).collect(Collectors.toList()));
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
        assertEquals(records, chunks.stream().flatMap(List::stream).collect(Collectors.toList()));
    }

    @Test(timeout = 15_000)
    public void testProducerEmitsAllRecordsThenEndOfStream() throws Exception
    {
        SAMFileHeader header = header();
        File bam = writeNameGroupedBam(header, 50);
        assertRecordCounts(records(50), produce(List.of(bam)));
    }

    @Test(timeout = 15_000)
    public void testProducerReadsEveryInputBamAndSignalsOnce() throws Exception
    {
        SAMFileHeader header = header();
        File firstBam = writeNameGroupedBam(header, 50, "lane_01.bam");
        File secondBam = writeNameGroupedBam(header, 30, "lane_02.bam");

        List<SAMRecord> expected = new ArrayList<>(records(50));
        expected.addAll(records(30));
        assertRecordCounts(expected, produce(List.of(firstBam, secondBam)));
    }

    private List<SAMRecord> produce(final List<File> bams) throws Exception
    {
        File ref = mFolder.newFile("ref.fa");

        int workers = 2;
        BlockingQueue<List<SAMRecord>> queue = new LinkedBlockingQueue<>(2);
        ShardedChunkProducer producer = new ShardedChunkProducer(
                bams.stream().map(File::getAbsolutePath).collect(Collectors.toList()), ref.getAbsolutePath(), queue, workers, 8, 4);
        producer.setDaemon(true);
        producer.start();

        int sentinels = 0;
        List<List<SAMRecord>> chunks = new ArrayList<>();
        while(sentinels < workers)
        {
            List<SAMRecord> chunk = queue.poll(5, TimeUnit.SECONDS);
            assertNotNull("producer did not deliver a chunk or end-of-stream", chunk);
            if(chunk == ShardedChunkProducer.END_OF_STREAM)
            {
                ++sentinels;
            }
            else
            {
                assertEquals("data after end-of-stream", 0, sentinels);
                assertFalse("empty data chunk", chunk.isEmpty());
                chunks.add(chunk);
            }
        }
        producer.join(5_000);

        assertFalse("producer did not terminate", producer.isAlive());
        assertEquals(workers, sentinels);
        assertTrue("unexpected chunks or end-of-stream markers", queue.isEmpty());

        for(List<SAMRecord> chunk : chunks)
        {
            for(int i = 0; i < chunk.size(); i += 2)
            {
                assertTrue("orphaned mate in chunk", i + 1 < chunk.size());
                assertEquals(chunk.get(i).getReadName(), chunk.get(i + 1).getReadName());
                assertTrue(chunk.get(i).getFirstOfPairFlag());
                assertTrue(chunk.get(i + 1).getSecondOfPairFlag());
            }
        }
        return chunks.stream().flatMap(List::stream).collect(Collectors.toList());
    }

    private static List<String> samStrings(final List<SAMRecord> records)
    {
        return records.stream().map(SAMRecord::getSAMString).collect(Collectors.toList());
    }

    private static void assertRecordCounts(final List<SAMRecord> expected, final List<SAMRecord> actual)
    {
        assertEquals(
                samStrings(expected).stream().collect(Collectors.groupingBy(Function.identity(), Collectors.counting())),
                samStrings(actual).stream().collect(Collectors.groupingBy(Function.identity(), Collectors.counting())));
    }

    @Test(timeout = 15_000)
    public void testMultipleShardsPreserveEveryRecordAndKeepNameGroupsWhole() throws Exception
    {
        SAMFileHeader header = header();
        List<SAMRecord> expected = records(5_000);
        File bam = mFolder.newFile("multi_block.bam");
        writeBam(bam, header, expected);
        List<ShardRange> ranges = ShardedChunkProducer.computeSplits(bam, header, 8);
        assertTrue("fixture must exercise multiple BGZF shard boundaries", ranges.size() > 1);

        List<SAMRecord> actual = new ArrayList<>();
        Map<String, Integer> shardByName = new HashMap<>();
        for(int shard = 0; shard < ranges.size(); ++shard)
        {
            try(ShardRecordIterator iterator = new ShardRecordIterator(bam, header, ranges.get(shard)))
            {
                while(iterator.hasNext())
                {
                    SAMRecord record = iterator.next();
                    Integer previousShard = shardByName.putIfAbsent(record.getReadName(), shard);
                    if(previousShard != null)
                    {
                        assertEquals("fragment split across shards", shard, previousShard.intValue());
                    }
                    actual.add(record);
                }
            }
        }
        assertEquals(samStrings(expected), samStrings(actual));
        assertRecordCounts(expected, produce(List.of(bam)));
    }

    @Test
    public void testNameGroupSpanningMultipleBgzfBlocksStaysInOneShard() throws Exception
    {
        SAMFileHeader header = header();
        List<SAMRecord> expected = records(5_000);
        for(int i = 0; i < 4_000; ++i)
        {
            expected.get(i).setReadName("largeGroup");
        }
        File bam = mFolder.newFile("large_group.bam");
        writeBam(bam, header, expected);
        List<ShardRange> ranges = ShardedChunkProducer.computeSplits(bam, header, 16);
        assertTrue("fixture must have multiple shards", ranges.size() > 1);
        List<List<SAMRecord>> shards = new ArrayList<>();
        for(ShardRange range : ranges)
        {
            List<SAMRecord> records = new ArrayList<>();
            try(ShardRecordIterator iterator = new ShardRecordIterator(bam, header, range))
            {
                iterator.forEachRemaining(records::add);
            }
            shards.add(records);
        }
        assertWholeFragments(shards);
        assertEquals(samStrings(expected), samStrings(shards.stream().flatMap(List::stream).collect(Collectors.toList())));
    }

    private static List<SAMRecord> readShard(final File bam, final SAMFileHeader header, final ShardRange range)
            throws IOException
    {
        List<SAMRecord> records = new ArrayList<>();
        try(ShardRecordIterator iterator = new ShardRecordIterator(bam, header, range))
        {
            iterator.forEachRemaining(records::add);
        }
        return records;
    }

    @Test(expected = IllegalStateException.class)
    public void testShardEndingInsideRecordFails() throws IOException
    {
        SAMFileHeader header = header();
        File bam = writeNameGroupedBam(header, 5_000);
        ShardRange first = ShardedChunkProducer.computeSplits(bam, header, 8).get(0);

        readShard(bam, header, new ShardRange(first.startVptr(), first.endVptr() + 1, first.endReadName()));
    }

    @Test(expected = IllegalStateException.class)
    public void testShardEndingInsideReadGroupFails() throws IOException
    {
        SAMFileHeader header = header();
        File bam = writeNameGroupedBam(header, 5_000);
        ShardRange first = ShardedChunkProducer.computeSplits(bam, header, 8).get(0);
        List<SAMRecord> firstShard = readShard(bam, header, first);
        String lastReadName = firstShard.get(firstShard.size() - 1).getReadName();

        readShard(bam, header, new ShardRange(first.startVptr(), first.endVptr(), lastReadName));
    }

    @Test
    public void testEmptyInputDoesNotEmitChunks() throws Exception
    {
        assertTrue(chunk(List.of(), 2).isEmpty());
    }

    @Test(timeout = 15_000)
    public void testEmptyBamStillSignalsEveryWorker() throws Exception
    {
        File bam = writeNameGroupedBam(header(), 0);
        assertTrue(produce(List.of(bam)).isEmpty());
    }

    @Test
    public void testReadCounterIncludesTheFinalPartialChunk() throws Exception
    {
        List<SAMRecord> records = List.of(read("a"), read("a"), read("b"));
        BlockingQueue<List<SAMRecord>> queue = new LinkedBlockingQueue<>();
        LongAdder counter = new LongAdder();
        ShardedChunkProducer.streamChunks(records.iterator(), 2, queue, counter);

        assertEquals(3, counter.sum());
        assertEquals(records, queue.stream().flatMap(List::stream).collect(Collectors.toList()));
    }
}
