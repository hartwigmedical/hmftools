package com.hartwig.hmftools.virusdetect;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMRecord;

// One extraction worker's private FASTA shard, written lock-free by its owning thread and later concatenated
// byte-wise into the candidate FASTA. Not thread-safe.
class FastaPart
{
    private final Path mPath;
    private final BufferedWriter mWriter;
    private int mReadCount;
    private boolean mClosed;

    private static final Logger LOGGER = LogManager.getLogger(FastaPart.class);

    static FastaPart create(String outputFastaFile, int index)
    {
        String path = outputFastaFile + ".part" + index;
        try
        {
            return new FastaPart(path, createBufferedWriter(path));
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to create candidate FASTA part", e);
        }
    }

    private FastaPart(String path, BufferedWriter writer)
    {
        mPath = new File(path).toPath();
        mWriter = writer;
    }

    void add(SAMRecord record)
    {
        try
        {
            mWriter.write(">" + fastaLabel(record));
            mWriter.newLine();
            mWriter.write(record.getReadString());
            mWriter.newLine();
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to write candidate reads", e);
        }
        ++mReadCount;
    }

    Path path()
    {
        return mPath;
    }

    int readCount()
    {
        return mReadCount;
    }

    void close() throws IOException
    {
        if(!mClosed)
        {
            mClosed = true;
            mWriter.close();
        }
    }

    // Releases the part and removes its file, whether or not it was consumed. Never throws, so it is safe while
    // unwinding a failed extraction.
    void discard()
    {
        try
        {
            close();
            Files.deleteIfExists(mPath);
        }
        catch(IOException e)
        {
            LOGGER.warn("failed to discard candidate FASTA part {}: {}", mPath, e.getMessage());
        }
    }

    // Mates differ only by the suffix: the FASTA is single-ended, so a pair's two reads must stay distinguishable.
    private static String fastaLabel(SAMRecord record)
    {
        if(record.getReadPairedFlag())
        {
            return record.getReadName() + (record.getFirstOfPairFlag() ? "/1" : "/2");
        }
        return record.getReadName();
    }
}
