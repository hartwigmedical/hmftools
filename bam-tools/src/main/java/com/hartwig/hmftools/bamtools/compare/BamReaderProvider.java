package com.hartwig.hmftools.bamtools.compare;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

// interface to abstract away the bam reader. Main use is to
// allow thread local bam reader
public interface BamReaderProvider extends AutoCloseable
{
    SamReader getBamReader();

    // override to remove the checked exception
    @Override void close();

    // an implementation that simply wraps the input one
    static BamReaderProvider of(SamReader samReader)
    {
        return new SimpleBamReaderProvider(samReader);
    }

    // opens a sam reader per thread
    static BamReaderProvider makeThreadLocal(SamReaderFactory samReaderFactory, File bamFile)
    {
        return new ThreadLocalBamReader(samReaderFactory, bamFile);
    }

    class SimpleBamReaderProvider implements BamReaderProvider
    {
        private final SamReader mBamReader;

        public SimpleBamReaderProvider(final SamReader bamReader)
        {
            this.mBamReader = bamReader;
        }

        @Override
        public SamReader getBamReader() { return mBamReader; }

        @Override
        public void close()
        {
            try
            {
                mBamReader.close();
            }
            catch(IOException e)
            {
                throw new UncheckedIOException(e);
            }
        }
    }

    class ThreadLocalBamReader implements BamReaderProvider
    {
        private final ThreadLocal<SamReader> mThreadLocalBamReaders;

        private final List<SamReader> mSamReaderList = Collections.synchronizedList(new ArrayList<>());

        public ThreadLocalBamReader(SamReaderFactory samReaderFactory, File bamFile)
        {
            mThreadLocalBamReaders = ThreadLocal.withInitial(() -> {
                SamReader bamReader = samReaderFactory.open(bamFile);
                mSamReaderList.add(bamReader);
                return bamReader;
            });
        }

        @Override
        public SamReader getBamReader() { return mThreadLocalBamReaders.get(); }

        @Override
        public void close()
        {
            for(SamReader bamReader : mSamReaderList)
            {
                try
                {
                    bamReader.close();
                }
                catch(IOException e)
                {
                    throw new UncheckedIOException(e);
                }
            }
            mSamReaderList.clear();
        }
    }
}
