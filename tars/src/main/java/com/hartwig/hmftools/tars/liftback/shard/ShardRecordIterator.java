package com.hartwig.hmftools.tars.liftback.shard;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Iterator;
import java.util.NoSuchElementException;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.BlockCompressedInputStream;

// Decodes the BAM records in one shard range [startVptr, endVptr). Both ends are read-name-group boundaries
// (see BamShardSplitter), so stopping when the next record begins at/after endVptr leaves every group whole.
public class ShardRecordIterator implements Iterator<SAMRecord>, Closeable
{
    private final BlockCompressedInputStream mStream;
    private final BAMRecordCodec mCodec;
    private final long mEndVptr;
    private final long mStartOffset;
    private SAMRecord mNext;

    public ShardRecordIterator(final File bam, final SAMFileHeader header, final BamShardSplitter.ShardRange range)
    {
        try
        {
            mStream = new BlockCompressedInputStream(bam);
            mStream.seek(range.startVptr());
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }

        mCodec = new BAMRecordCodec(header);
        mCodec.setInputStream(mStream);
        mEndVptr = range.endVptr();
        mStartOffset = range.startVptr() >>> 16;
        mNext = readNext();
    }

    // compressed bytes this shard has consumed so far -- the progress monitor sums these across shards.
    public long consumedBytes()
    {
        return Math.max(0, (mStream.getFilePointer() >>> 16) - mStartOffset);
    }

    private SAMRecord readNext()
    {
        if(mEndVptr != BamShardSplitter.EOF && mStream.getFilePointer() >= mEndVptr)
        {
            return null;
        }

        return mCodec.decode(); // null at EOF
    }

    @Override
    public boolean hasNext()
    {
        return mNext != null;
    }

    @Override
    public SAMRecord next()
    {
        if(mNext == null)
        {
            throw new NoSuchElementException();
        }

        SAMRecord record = mNext;
        mNext = readNext();
        return record;
    }

    @Override
    public void close() throws IOException
    {
        mStream.close();
    }
}
