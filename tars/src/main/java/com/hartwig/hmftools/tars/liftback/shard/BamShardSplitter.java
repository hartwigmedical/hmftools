package com.hartwig.hmftools.tars.liftback.shard;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.BlockCompressedInputStream;

// Splits a name-grouped BAM into N byte ranges, each starting on a read-name boundary so a shard reads
// [start, end) with whole fragments and no cross-shard reconciliation. A split is a BGZF virtual pointer at a
// group's first record: from a target byte offset, find the next block (by its magic), decode until the name
// changes, and take that record's pointer. Shard 0 starts after the header; the last runs to EOF.
public final class BamShardSplitter
{
    private BamShardSplitter() {}

    // [startVptr, endVptr); endVptr == EOF means "read to end of file".
    public record ShardRange(long startVptr, long endVptr) {}

    public static final long EOF = Long.MAX_VALUE;

    private static final byte[] BGZF_MAGIC = { 0x1f, (byte) 0x8b, 0x08, 0x04 };
    private static final int SCAN_BUFFER = 64 * 1024;

    public static List<ShardRange> computeSplits(final File bam, final SAMFileHeader header, final int shardCount)
            throws IOException
    {
        final long firstRecordVptr = headerEndVptr(bam, header);
        if(shardCount <= 1)
            return List.of(new ShardRange(firstRecordVptr, EOF));

        final long compressedLength = bam.length();
        final List<Long> splits = new ArrayList<>();
        splits.add(firstRecordVptr);

        long lastBlockOffset = firstRecordVptr >>> 16;
        for(int k = 1; k < shardCount; ++k)
        {
            final long targetOffset = compressedLength * k / shardCount;
            if(targetOffset <= lastBlockOffset)
                continue; // shards collapsing onto the same block: fewer, larger shards is fine

            final long blockStart = findBlockStart(bam, targetOffset);
            if(blockStart < 0)
                break; // no further block boundary: the remaining file is one shard

            final long splitVptr = firstGroupBoundaryVptr(bam, header, blockStart << 16);
            if(splitVptr > splits.get(splits.size() - 1))
            {
                splits.add(splitVptr);
                lastBlockOffset = splitVptr >>> 16;
            }
        }

        final List<ShardRange> ranges = new ArrayList<>();
        for(int i = 0; i < splits.size(); ++i)
        {
            final long end = i + 1 < splits.size() ? splits.get(i + 1) : EOF;
            ranges.add(new ShardRange(splits.get(i), end));
        }
        return ranges;
    }

    // walk the BAM header binary on a BGZF stream and return the virtual pointer of the first alignment record.
    private static long headerEndVptr(final File bam, final SAMFileHeader header) throws IOException
    {
        try(BlockCompressedInputStream stream = new BlockCompressedInputStream(bam))
        {
            final BinaryCodec codec = new BinaryCodec(stream);
            codec.readBytes(new byte[4]);             // "BAM\1"
            final int textLength = codec.readInt();
            codec.readBytes(new byte[textLength]);    // header text
            final int refCount = codec.readInt();
            for(int i = 0; i < refCount; ++i)
            {
                final int nameLength = codec.readInt();
                codec.readBytes(new byte[nameLength]); // reference name
                codec.readInt();                      // reference length
            }
            return stream.getFilePointer();
        }
    }

    // scan raw bytes from targetOffset for the next valid BGZF block header, return its compressed offset (-1 none).
    private static long findBlockStart(final File bam, final long targetOffset) throws IOException
    {
        try(RandomAccessFile raf = new RandomAccessFile(bam, "r"))
        {
            final long length = raf.length();
            final byte[] buffer = new byte[SCAN_BUFFER];
            long offset = targetOffset;

            while(offset < length - 18)
            {
                raf.seek(offset);
                final int read = raf.read(buffer);
                if(read <= 0)
                    return -1;

                for(int i = 0; i + 18 <= read; ++i)
                {
                    if(buffer[i] == BGZF_MAGIC[0] && buffer[i + 1] == BGZF_MAGIC[1]
                            && buffer[i + 2] == BGZF_MAGIC[2] && buffer[i + 3] == BGZF_MAGIC[3]
                            && buffer[i + 12] == 0x42 && buffer[i + 13] == 0x43) // 'BC' subfield id
                        return offset + i;
                }

                offset += read - 17; // overlap so a header straddling the buffer edge is not missed
            }
        }
        return -1;
    }

    // decode records from blockVptr tracking names; return the virtual pointer of the first record whose name
    // differs from the first record seen -- i.e. the start of the second group at/after this block.
    private static long firstGroupBoundaryVptr(final File bam, final SAMFileHeader header, final long blockVptr)
            throws IOException
    {
        try(BlockCompressedInputStream stream = new BlockCompressedInputStream(bam))
        {
            stream.seek(blockVptr);
            final BAMRecordCodec codec = new BAMRecordCodec(header);
            codec.setInputStream(stream);

            String firstName = null;
            long recordVptr = stream.getFilePointer();
            SAMRecord record;
            while((record = codec.decode()) != null)
            {
                if(firstName == null)
                    firstName = record.getReadName();
                else if(!record.getReadName().equals(firstName))
                    return recordVptr;

                recordVptr = stream.getFilePointer();
            }
            return EOF; // no further boundary before EOF
        }
    }
}
