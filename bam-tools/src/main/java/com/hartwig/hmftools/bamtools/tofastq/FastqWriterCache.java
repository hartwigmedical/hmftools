package com.hartwig.hmftools.bamtools.tofastq;

import static java.lang.String.format;

import java.io.File;
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

public class FastqWriterCache
{
    private final FastqConfig mConfig;

    private final Map<String,FastqWriter> mReadGroupWriters;
    private final List<FastqWriter> mThreadedWriters;

    public FastqWriterCache(final FastqConfig config)
    {
        mConfig = config;

        mReadGroupWriters = Maps.newHashMap();
        mThreadedWriters = Lists.newArrayList();

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
        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mConfig.BamFile));

        final SAMFileHeader fileHeader = samReader.getFileHeader();

        for(SAMReadGroupRecord readGroupRecord : fileHeader.getReadGroups())
        {
            String readGroupId = readGroupRecord.getId();
            String filename = mConfig.formFilePrefix(readGroupId);
            FastqWriter fastqWriter = new FastqWriter(filename);
            mReadGroupWriters.put(readGroupId, fastqWriter);
        }
    }

    public FastqWriter createThreadedWriter()
    {
        String fileId = mConfig.Threads >  1 ? format("t%d", mThreadedWriters.size()) : "";
        String filename = mConfig.formFilePrefix(fileId);
        FastqWriter fastqWriter = new FastqWriter(filename);
        mThreadedWriters.add(fastqWriter);
        return fastqWriter;
    }

    public void processReadPair(final SAMRecord first, final SAMRecord second, @Nullable final FastqWriter fastqWriter)
    {
        if(fastqWriter != null)
        {
            fastqWriter.processReadPairNoSync(first, second);
            return;
        }

        if(mConfig.SplitMode == FileSplitMode.READ_GROUP)
        {
            FastqWriter writer = mReadGroupWriters.get(first.getReadGroup().getReadGroupId());

            if(writer != null)
                writer.processReadPairSync(first, second);
        }
        else
        {
            mThreadedWriters.get(0).processReadPairSync(first, second);
        }
    }

    public void processUnpairedRead(final SAMRecord read, @Nullable final FastqWriter fastqWriter)
    {
        if(fastqWriter != null)
        {
            fastqWriter.processUnpairedRead(read);
            return;
        }

        if(mConfig.SplitMode == FileSplitMode.READ_GROUP)
        {
            FastqWriter writer = mReadGroupWriters.get(read.getReadGroup().getReadGroupId());

            if(writer != null)
                writer.processUnpairedRead(read);
        }
        else
        {
            mThreadedWriters.get(0).processUnpairedRead(read);
        }
    }

    public void close()
    {
        if(!mThreadedWriters.isEmpty())
        {
            mThreadedWriters.forEach(x -> x.close());
        }
        else
        {
            mReadGroupWriters.values().forEach(x -> x.close());
        }
    }
}
