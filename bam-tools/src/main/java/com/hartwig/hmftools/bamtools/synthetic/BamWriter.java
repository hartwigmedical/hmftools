package com.hartwig.hmftools.bamtools.synthetic;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bamops.BamToolName.fromPath;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import com.hartwig.hmftools.bamtools.slice.IReadWriter;
import com.hartwig.hmftools.bamtools.slice.WriteType;
import com.hartwig.hmftools.common.bamops.BamOperations;
import com.hartwig.hmftools.common.bamops.BamToolName;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamWriter implements IReadWriter
{
    private final SyntheticConfig mConfig;

    private long mRecordWriteCount;
    private final SAMFileWriter mBamWriter;
    private String mUnsortedOutputBam;
    private String mFinalOutputBam;

    public BamWriter(final SyntheticConfig config)
    {
        mConfig = config;
        mRecordWriteCount = 0;
        mBamWriter = config.WriteBam ? initialiseBam() : null;
    }

    private SAMFileWriter initialiseBam()
    {
        SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));

        mFinalOutputBam = mConfig.formFilename(BAM_EXTENSION);

        int bamExtension = mFinalOutputBam.lastIndexOf(BAM_EXTENSION);
        mUnsortedOutputBam = mFinalOutputBam.substring(0, bamExtension) + ".unsorted.bam";

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();

        fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(mUnsortedOutputBam));
    }

    public synchronized void writeRead(final SAMRecord read) { doWriteRecord(read); }

    public synchronized long writeCount() { return mRecordWriteCount; }

    private void doWriteRecord(final SAMRecord read)
    {
        ++mRecordWriteCount;

        if(mBamWriter != null)
            mBamWriter.addAlignment(read);
    }

    public void close()
    {
        if(mBamWriter != null)
        {
            BT_LOGGER.info("{} records written to BAM", mRecordWriteCount);

            mBamWriter.close();

            BT_LOGGER.info("writing sorted BAM: {}", mFinalOutputBam);

            BamToolName toolName = fromPath(mConfig.BamToolPath);

            boolean success = BamOperations.sortBam(toolName, mConfig.BamToolPath, mUnsortedOutputBam, mFinalOutputBam, mConfig.Threads);

            if(!success)
            {
                BT_LOGGER.error("failed to sort BAM");
                System.exit(1);
            }

            if(toolName == BamToolName.SAMTOOLS)
            {
                success = BamOperations.indexBam(toolName, mConfig.BamToolPath, mFinalOutputBam, mConfig.Threads);
            }

            if(success)
            {
                try
                {
                    Files.deleteIfExists(Paths.get(mUnsortedOutputBam));
                }
                catch(IOException e)
                {
                    BT_LOGGER.error("error deleting interim file: {}", e.toString());
                }
            }
        }
    }
}
