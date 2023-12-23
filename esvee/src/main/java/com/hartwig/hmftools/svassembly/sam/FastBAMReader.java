package com.hartwig.hmftools.svassembly.sam;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import com.hartwig.hmftools.svassembly.models.Record;

public class FastBAMReader implements SAMSource, AutoCloseable
{
    private final RandomAccessFile mFile;
    private final ByteBuffer mContents;
    private String mHeader;
    private Map<String, Integer> mChromosomeLengths;
    private final FastBAMIndex mIndex;

    public FastBAMReader(final File bamFile, final File bamIndex) throws IOException
    {
        mFile = new RandomAccessFile(bamFile, "r");
        final FileChannel channel = mFile.getChannel();
        final MappedByteBuffer mapped = channel.map(FileChannel.MapMode.READ_ONLY, 0, mFile.length());
        mapped.order(ByteOrder.LITTLE_ENDIAN);
        mContents = mapped.asReadOnlyBuffer();
        open();

        mIndex = new FastBAMIndex(bamIndex, new ArrayList<>(mChromosomeLengths.keySet()));
    }

    private synchronized void open()
    {
        final int magic = mContents.getInt();
        final int headTextLength = mContents.getInt();
        final byte[] headerBytes = new byte[headTextLength];
        mContents.get(headerBytes);
        final String headerText = new String(headerBytes);

        final int chromosomeCount = mContents.getInt();
        mChromosomeLengths = new LinkedHashMap<>();
        for (int i = 0; i < chromosomeCount; i++)
        {
            final int nameLength = mContents.getInt();
            final byte[] nameWithoutNullTerminator = new byte[nameLength - 1];
            mContents.get(nameWithoutNullTerminator);
            mContents.get(); // Consume null terminator
            final String chromosomeName = new String(nameWithoutNullTerminator);

            final int length = mContents.getInt();
            mChromosomeLengths.put(chromosomeName, length);
        }
    }

    public List<String> getChromosomeNames()
    {
        return new ArrayList<>(mChromosomeLengths.keySet());
    }

    @Override
    public Stream<Record> unmappedReads()
    {
        return null;
    }

    @Override
    public Stream<Record> streamReadsContaining(final String chromosome, final int startPosition, final int endPosition)
    {
        return null;
    }

    @Override
    public void close()
    {
        try
        {
            mFile.close();
        }
        catch(final IOException exception)
        {
            throw new RuntimeException(exception);
        }
    }
}
