package com.hartwig.hmftools.bamtools.tofastq;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

// a simple class to manage thread local data
public class ThreadData implements AutoCloseable
{
    private final ThreadLocal<FastqWriter> mFastqWriters;
    private final ThreadLocal<PartitionReader> mPartitionReaders;
    private final List<PartitionReader> mPartitionReaderList = Collections.synchronizedList(new ArrayList<>());

    public List<PartitionReader> getPartitionReaders()
    {
        return mPartitionReaderList;
    }

    public ThreadData(final FastqConfig config, final FastqWriterCache writerCache, final RemoteReadHandler remoteReadHandler)
    {
        mFastqWriters = ThreadLocal.withInitial(() -> config.SplitMode == FileSplitMode.THREAD ? writerCache.createThreadedWriter() : null);

        // make a PartitionReader per thread
        mPartitionReaders = ThreadLocal.withInitial(() -> {
            PartitionReader partitionReader = new PartitionReader(config, writerCache, remoteReadHandler, getFastqWriter());
            mPartitionReaderList.add(partitionReader);
            return partitionReader;
        });
    }

    public PartitionReader getPartitionReader()
    {
        return mPartitionReaders.get();
    }

    public FastqWriter getFastqWriter()
    {
        return mFastqWriters.get();
    }

    @Override
    public void close()
    {
        for(PartitionReader partitionReader : mPartitionReaderList)
        {
            partitionReader.close();
        }
    }
}
