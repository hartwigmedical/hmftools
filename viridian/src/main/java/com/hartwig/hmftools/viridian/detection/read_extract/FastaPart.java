package com.hartwig.hmftools.viridian.detection.read_extract;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMRecord;

// One candidate read extraction worker's private FASTA shard.
// Written lock-free by its owning thread and later concatenated into the candidate FASTA.
// Not thread-safe.
public class FastaPart
{
    private final Path mPath;
    private final BufferedWriter mWriter;
    private int mReadCount;
    private boolean mClosed;

    private static final Logger LOGGER = LogManager.getLogger(FastaPart.class);

    static FastaPart create(String outputFastaFile, int index)
    {
        Path path = new File(outputFastaFile + ".part" + index).toPath();
        try
        {
            return new FastaPart(path, createBufferedWriter(path.toString()));
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to create candidate FASTA part", e);
        }
    }

    private FastaPart(Path path, BufferedWriter writer)
    {
        mPath = path;
        mWriter = writer;
    }

    public Path path() { return mPath; }

    public int readCount() { return mReadCount; }

    public void add(SAMRecord record)
    {
        try
        {
            // The FASTA is single-ended, so a pair's two reads must stay distinguishable.
            String suffix = record.getReadPairedFlag() ? (record.getFirstOfPairFlag() ? "/1" : "/2") : "";
            mWriter.write(">" + record.getReadName() + suffix);
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

    public void close() throws IOException
    {
        if(!mClosed)
        {
            mClosed = true;
            mWriter.close();
        }
    }

    // Releases the shard and removes its file, whether or not it was consumed. Never throws, so it is safe while
    // unwinding a failed extraction.
    public void discard()
    {
        try
        {
            close();
            Files.deleteIfExists(mPath);
        }
        catch(IOException e)
        {
            LOGGER.warn("Failed to discard candidate FASTA part {}: {}", mPath, e.getMessage());
        }
    }
}
