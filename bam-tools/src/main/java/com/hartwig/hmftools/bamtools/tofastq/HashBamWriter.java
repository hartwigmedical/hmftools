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
import htsjdk.samtools.SamReaderFactory;

// a class to write reads into bam files hashed by the read names
public class HashBamWriter implements AutoCloseable
{
    private static final int NUM_HASH_BAMS = 256;

    private final FastqConfig mConfig;

    private final File[] mHashBams = new File[NUM_HASH_BAMS];
    private final SAMFileWriter[] mHashBamWriters = new SAMFileWriter[NUM_HASH_BAMS];

    public int numHashBams() { return NUM_HASH_BAMS; }

    public File[] getHashBams()
    {
        return mHashBams;
    }

    public HashBamWriter(final FastqConfig config)
    {
        mConfig = config;
        try(SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mConfig.BamFile)))
        {
            SAMFileHeader samHeader = samReader.getFileHeader();
            samHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            BT_LOGGER.info("creating {} hash bams", NUM_HASH_BAMS);

            // create a temp directory to house the files
            String tempDir = Files.createTempDirectory(Paths.get(mConfig.BamFile).getFileName() + "_hash_").toString();

            for(int i = 0; i < NUM_HASH_BAMS; ++i)
            {
                File tempBamFile = new File(tempDir, String.format("%d.hash.bam", i));
                tempBamFile.deleteOnExit();
                mHashBamWriters[i] = new SAMFileWriterFactory().makeBAMWriter(samHeader, false, tempBamFile);
                mHashBams[i] = tempBamFile;
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

    public void writeToHashBam(SAMRecord read)
    {
        // we put the lock on individual sam writer
        SAMFileWriter bamWriter = mHashBamWriters[getHashBamIndex(read.getReadName())];
        synchronized(bamWriter)
        {
            bamWriter.addAlignment(read);
        }
    }

    private static int getHashBamIndex(final String readName)
    {
        return Math.abs(readName.hashCode()) % NUM_HASH_BAMS;
    }

    @Override
    public void close()
    {
        Arrays.stream(mHashBamWriters).forEach(SAMFileWriter::close);
    }
}
