package com.hartwig.hmftools.bamtools.tofastq;

import static java.lang.String.format;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class FastqWriterCache implements AutoCloseable
{
    private final FastqConfig mConfig;

    private final Map<String,FastqWriter> mReadGroupWriters;
    private final List<FastqWriter> mWriters;

    public FastqWriterCache(final FastqConfig config)
    {
        mConfig = config;

        mReadGroupWriters = Maps.newHashMap();
        mWriters = Lists.newArrayList();

        if(mConfig.SplitMode == FileSplitMode.READ_GROUP)
        {
            createReadGroupWriters();
        }
        else if(mConfig.SplitMode == FileSplitMode.NONE)
        {
            createThreadedWriter();
        }
    }

    private void createReadGroupWriters()
    {
        // create from BAM header
        final SAMFileHeader fileHeader;
        try(SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mConfig.BamFile)))
        {
            fileHeader = samReader.getFileHeader();
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }

        for(SAMReadGroupRecord readGroupRecord : fileHeader.getReadGroups())
        {
            String readGroupId = readGroupRecord.getId();
            String filePrefix = mConfig.OutputDir + readGroupId;
            FastqWriter fastqWriter = new FastqWriter(filePrefix, mConfig.WriteUnzipped);
            mReadGroupWriters.put(readGroupId, fastqWriter);
            mWriters.add(fastqWriter);
        }
    }

    public synchronized FastqWriter getThreadedWriter()
    {
        if(mWriters.isEmpty())
            return createThreadedWriter();

        return mWriters.get(0);
    }

    public synchronized FastqWriter createThreadedWriter()
    {
        String fileId = mConfig.Threads > 1 ? format("t%d", mWriters.size()) : "";
        String filePrefix = mConfig.formFilePrefix(fileId);
        FastqWriter fastqWriter = new FastqWriter(filePrefix, mConfig.WriteUnzipped);
        mWriters.add(fastqWriter);
        return fastqWriter;
    }

    public void processReadPair(final SAMRecord first, final SAMRecord second, @Nullable final FastqWriter fastqWriter)
    {
        if(fastqWriter != null)
        {
            fastqWriter.processReadPairNoSync(first, second);
            return;
        }

        FastqWriter writer = getWriter(first);
        writer.processReadPairSync(first, second);
    }

    public void processUnpairedRead(final SAMRecord read, @Nullable final FastqWriter fastqWriter)
    {
        if(fastqWriter != null)
        {
            fastqWriter.processUnpairedRead(read);
            return;
        }

        FastqWriter writer = getWriter(read);
        writer.processUnpairedRead(read);
    }

    public void writeUnpairedRead(final SAMRecord read)
    {
        FastqWriter writer = getWriter(read);
        writer.writeUnpairedRead(read);
    }

    private FastqWriter getWriter(final SAMRecord read)
    {
        return mConfig.SplitMode == FileSplitMode.READ_GROUP ?
                mReadGroupWriters.get(read.getReadGroup().getReadGroupId()) : mWriters.get(0);
    }

    public void writeUnpairedReads()
    {
        mWriters.forEach(x -> x.writeUnpairedReads());
    }

    @Override
    public void close()
    {
        mWriters.forEach(x -> x.close());
    }
}
