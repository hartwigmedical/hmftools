package com.hartwig.hmftools.bamtools.tofastq;

import java.util.HashMap;
import java.util.Map;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

public class FastqWriterCache implements AutoCloseable
{
    private final ToFastqConfig mConfig;
    private final String mThreadId;

    private long mNumReadsWritten;

    private final Map<String, FastqWriter> mReadGroupWriters;
    @Nullable final private FastqWriter mSingleWriter;

    public long numReadsWritten() { return mNumReadsWritten; }
    @Nullable public FastqWriter getFastqWriter() { return mSingleWriter; }

    @Nullable
    public FastqWriter getReadGroupFastqWriter(String readGroupId)
    {
        return mReadGroupWriters.get(readGroupId);
    }

    public FastqWriterCache(final ToFastqConfig config, String threadId)
    {
        mConfig = config;
        mThreadId = threadId;
        mReadGroupWriters = new HashMap<>();

        if(mConfig.SplitMode == FileSplitMode.READ_GROUP)
        {
            mSingleWriter = null;
            createReadGroupWriters();
        }
        else
        {
            mSingleWriter = new FastqWriter(mConfig.formFilePrefix(threadId, "", true));
        }
    }

    private void createReadGroupWriters()
    {
        for(SAMReadGroupRecord readGroup : ToFastqUtils.getReadGroups(mConfig))
        {
            FastqWriter fastqWriter = new FastqWriter(mConfig.formFilePrefix(mThreadId, readGroup.getId(), false));
            mReadGroupWriters.put(readGroup.getId(), fastqWriter);
        }
    }

    public void writeReadPair(final SAMRecord first, final SAMRecord second)
    {
        FastqWriter writer = getWriter(first);
        writer.writeReadPair(first, second);
        mNumReadsWritten += 2;
    }

    public void writeUnpairedRead(final SAMRecord read)
    {
        FastqWriter writer = getWriter(read);
        writer.writeUnpairedRead(read);
        mNumReadsWritten++;
    }

    private FastqWriter getWriter(final SAMRecord read)
    {
        return mSingleWriter == null ? mReadGroupWriters.get(read.getReadGroup().getReadGroupId()) : mSingleWriter;
    }

    @Override
    public void close()
    {
        if(mSingleWriter != null)
            mSingleWriter.close();
        mReadGroupWriters.values().forEach(FastqWriter::close);
    }
}
