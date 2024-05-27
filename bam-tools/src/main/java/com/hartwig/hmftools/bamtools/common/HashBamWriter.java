package com.hartwig.hmftools.bamtools.common;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;

// a class to write reads into bam files hashed by the read names
public class HashBamWriter implements AutoCloseable
{
    private final int mNumHashBams;

    private final File[] mHashBams;
    private final SAMFileWriter[] mHashBamWriters;

    public int numHashBams() { return mNumHashBams; }

    public Map<Integer, File> getHashBams()
    {
        // convert to key / file map
        Map<Integer, File> hashBamMap = new HashMap<>();
        for(int i = 0; i < mHashBams.length; ++i)
        {
            hashBamMap.put(i, mHashBams[i]);
        }
        return hashBamMap;
    }

    public HashBamWriter(final SAMFileHeader samHeader, String tempDirPrefix, int numHashBams)
    {
        mNumHashBams = numHashBams;
        mHashBams = new File[mNumHashBams];
        mHashBamWriters = new SAMFileWriter[mNumHashBams];

        SAMFileHeader newSamHeader = samHeader.clone();
        newSamHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        try
        {
            // create a temp directory to house the files
            final String tempDir = Files.createTempDirectory(tempDirPrefix)
                    .toString();

            BT_LOGGER.info("creating {} hash bams in {}", mNumHashBams, tempDir);

            for(int i = 0; i < mNumHashBams; ++i)
            {
                File hashBamFile = new File(tempDir, String.format("%d.hash.bam", i));
                mHashBamWriters[i] = new SAMFileWriterFactory().makeBAMWriter(newSamHeader, false, hashBamFile);
                mHashBams[i] = hashBamFile;
            }

            // add a hook to delete the temp directory on exit
            Runtime.getRuntime().addShutdownHook(new Thread(() -> {
                Arrays.stream(mHashBams).forEach(File::delete);
                new File(tempDir).delete();
            } ));
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    private int hashBamIndex(final String readName)
    {
        return Math.abs(readName.hashCode()) % mNumHashBams;
    }

    public void writeToHashBam(SAMRecord read)
    {
        // we put the lock on individual bam writer
        SAMFileWriter bamWriter = mHashBamWriters[hashBamIndex(read.getReadName())];

        //noinspection SynchronizationOnLocalVariableOrMethodParameter
        synchronized(bamWriter)
        {
            bamWriter.addAlignment(read);
        }
    }

    @Override
    public void close()
    {
        Arrays.stream(mHashBamWriters).forEach(SAMFileWriter::close);
    }
}
