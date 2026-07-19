package com.hartwig.hmftools.tars.liftback.shard;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

// Regression for the small-input shard bug: firstGroupBoundaryVptr returns the EOF sentinel when a probe lands
// on the BGZF end marker (common when the input is smaller than shardCount blocks). EOF must never become a
// shard START -- ShardRecordIterator would then seek to Long.MAX_VALUE ("Invalid file pointer"). The read-back
// also asserts the classic invariant: every fragment's records stay within one shard.
public class BamShardSplitterTest
{
    @Rule
    public TemporaryFolder mFolder = new TemporaryFolder();

    private static SAMFileHeader header()
    {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(new SAMSequenceRecord("chr1", 1_000_000))));
        return header;
    }

    private static SAMRecord read(final SAMFileHeader header, final String name, final int start)
    {
        SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setReferenceIndex(0);
        record.setAlignmentStart(start);
        record.setReadUnmappedFlag(false);
        record.setCigarString("10M");
        record.setReadBases("ACGTACGTAC".getBytes());
        byte[] quals = new byte[10];
        Arrays.fill(quals, (byte) 30);
        record.setBaseQualities(quals);
        record.setMappingQuality(60);
        return record;
    }

    private File writeNameGroupedBam(final SAMFileHeader header, final int fragments) throws IOException
    {
        File bam = mFolder.newFile("small.bam");
        try(SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, true, bam))
        {
            for(int i = 0; i < fragments; ++i)
            {
                String name = String.format("frag%04d", i);
                writer.addAlignment(read(header, name, 100 + i));   // R1
                writer.addAlignment(read(header, name, 300 + i));   // R2, same name
            }
        }
        return bam;
    }

    @Test
    public void smallInputHighShardCountNeverStartsAtEof() throws IOException
    {
        SAMFileHeader header = header();
        File bam = writeNameGroupedBam(header, 40);   // tiny: far fewer blocks than shardCount

        // A high shard count forces target offsets past the single data block, so probes hit the BGZF end marker.
        List<BamShardSplitter.ShardRange> ranges = BamShardSplitter.computeSplits(bam, header, 64);

        assertFalse(ranges.isEmpty());
        for(BamShardSplitter.ShardRange range : ranges)
            assertFalse("EOF used as a shard start", range.startVptr() == BamShardSplitter.EOF);

        // Read every shard back: pre-fix this threw when a range started at EOF. Assert all records are recovered
        // and no fragment is split across a boundary.
        List<String> names = new ArrayList<>();
        for(BamShardSplitter.ShardRange range : ranges)
        {
            try(ShardRecordIterator iterator = new ShardRecordIterator(bam, header, range))
            {
                while(iterator.hasNext())
                    names.add(iterator.next().getReadName());
            }
        }

        assertEquals(80, names.size());   // 40 fragments x 2 reads, none dropped or double-counted
        for(int i = 0; i < names.size(); i += 2)
            assertEquals("fragment split across a shard boundary", names.get(i), names.get(i + 1));
    }

    @Test
    public void singleShardSpansWholeFile() throws IOException
    {
        SAMFileHeader header = header();
        File bam = writeNameGroupedBam(header, 10);

        List<BamShardSplitter.ShardRange> ranges = BamShardSplitter.computeSplits(bam, header, 1);
        assertEquals(1, ranges.size());
        assertTrue(ranges.get(0).endVptr() == BamShardSplitter.EOF);
    }
}
