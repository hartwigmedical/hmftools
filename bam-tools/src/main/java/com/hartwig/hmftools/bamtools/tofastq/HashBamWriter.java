package com.hartwig.hmftools.bamtools.tofastq;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

// a class to write reads into bam files hashed by the read names
public class HashBamWriter implements AutoCloseable
{
    private static final int NUM_HASH_BAMS = 256;

    private final ToFastqConfig mConfig;

    private final File[] mHashBams = new File[NUM_HASH_BAMS];
    private final SAMFileWriter[] mHashBamWriters = new SAMFileWriter[NUM_HASH_BAMS];

    public int numHashBams() { return NUM_HASH_BAMS; }

    public File[] getHashBams()
    {
        return mHashBams;
    }

    public HashBamWriter(final ToFastqConfig config)
    {
        mConfig = config;
        try(SamReader samReader = ToFastqUtils.openSamReader(mConfig))
        {
            SAMFileHeader samHeader = samReader.getFileHeader();
            samHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            // create a temp directory to house the files
            String tempDir = Files.createTempDirectory(
                    Paths.get(mConfig.BamFile).getFileName().toString().replace('.', '_') + "_hashbams_")
                    .toString();

            BT_LOGGER.info("creating {} hash bams in {}", NUM_HASH_BAMS, tempDir);

            for(int i = 0; i < NUM_HASH_BAMS; ++i)
            {
                File hashBamFile = new File(tempDir, String.format("%d.hash.bam", i));
                mHashBamWriters[i] = new SAMFileWriterFactory().makeBAMWriter(samHeader, false, hashBamFile, 1);
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

    private static int hashBamIndex(final String readName)
    {
        return Math.abs(readName.hashCode()) % NUM_HASH_BAMS;
    }

    public void writeToHashBam(SAMRecord read)
    {
        // we put the lock on individual bam writer
        SAMFileWriter bamWriter = mHashBamWriters[hashBamIndex(read.getReadName())];
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
