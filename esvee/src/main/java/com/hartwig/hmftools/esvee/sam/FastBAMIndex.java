package com.hartwig.hmftools.esvee.sam;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class FastBAMIndex
{
    private final File mFile;
    private final List<String> mChromosomeNames;

    private final Map<String, List<IndexedPageRecord>> mLinearIndex = new HashMap<>();
    private IndexedPageRecord mUnmappedStart;
    private IndexedPageRecord mUnmappedEnd;

    public FastBAMIndex(final File file, final List<String> chromosomeNames)
    {
        mFile = file;
        mChromosomeNames = chromosomeNames;
    }

    public List<IndexedPageRecord> pagesFor(final String chromosome, final int startPosition, final int endPosition)
    {
        if (mLinearIndex.isEmpty())
        {
            try
            {
                open();
            }
            catch(final IOException exception)
            {
                throw new RuntimeException(exception);
            }
        }

        final List<IndexedPageRecord> linearIndex = mLinearIndex.get(chromosome);
        final int linearChunkIndexStart = startPosition / (16 * 1024);
        final int linearChunkIndexEnd = endPosition / (16 * 1024);

        final List<IndexedPageRecord> pages = new ArrayList<>();
        for (int i = linearChunkIndexStart; i <= linearChunkIndexEnd; i++)
        {
            final IndexedPageRecord chunk = linearIndex.get(i);
            if(pages.isEmpty())
                pages.add(chunk);
            else
            {
                final IndexedPageRecord lastPage = pages.get(pages.size() - 1);
                if(lastPage.PageFileOffset != chunk.PageFileOffset)
                    pages.add(new IndexedPageRecord(chunk.PageFileOffset, 0));
            }
        }
        return pages;
    }

    public IndexedPageRecord unmappedPagesStart()
    {
        return mUnmappedStart;
    }

    public IndexedPageRecord unmappedPagesEnd()
    {
        return mUnmappedEnd;
    }

    private void open() throws IOException
    {
        try(final RandomAccessFile randomAccessFile = new RandomAccessFile(mFile.getAbsolutePath(), "r"))
        {
            final FileChannel channel = randomAccessFile.getChannel();
            final MappedByteBuffer mapped = channel.map(FileChannel.MapMode.READ_ONLY, 0, randomAccessFile.length());
            mapped.order(ByteOrder.LITTLE_ENDIAN);
            final ByteBuffer buffer = mapped.asReadOnlyBuffer();

            final int magic = buffer.getInt();
            final int referenceSequenceCount = buffer.getInt();
            for(int chromosomeIndex = 0; chromosomeIndex < referenceSequenceCount; chromosomeIndex++)
            {
                final int binCount = buffer.getInt();
                for(int binIndex = 0; binIndex < binCount; binIndex++)
                {
                    final int bin = buffer.getInt();
                    final int chunkCount = buffer.getInt();

                    for(int chunkIndex = 0; chunkIndex < chunkCount; chunkIndex++)
                    {
                        final long chunkStart = buffer.getLong();
                        final long chunkEnd = buffer.getLong();

                        // We just ignore everything except the linear index
                    }
                }
                final int linearIndexSize = buffer.getInt();
                final List<IndexedPageRecord> linearIndex = new ArrayList<>(linearIndexSize);
                for(int i = 0; i < linearIndexSize; i++)
                    linearIndex.add(new IndexedPageRecord(buffer.getLong()));
                mLinearIndex.put(mChromosomeNames.get(chromosomeIndex), linearIndex);
            }
        }
    }

    public static class IndexedPageRecord
    {
        public final long PageFileOffset;
        public final long OffsetInPage;

        public IndexedPageRecord(final long virtualOffset)
        {
            this(virtualOffset >>> 16, virtualOffset & 0x00FFFF);
        }

        public IndexedPageRecord(final long pageFileOffset, final long offsetInPage)
        {
            PageFileOffset = pageFileOffset;
            OffsetInPage = offsetInPage;
        }
    }
}
