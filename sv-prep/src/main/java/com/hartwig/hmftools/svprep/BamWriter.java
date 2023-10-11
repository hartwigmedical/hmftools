package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.WriteType.BAM;

import java.io.File;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamWriter
{
    private final SvConfig mConfig;

    private int mRecordWriteCount;
    private final SAMFileWriter mWriter;
    private String mOutputBam;

    public BamWriter(final SvConfig config)
    {
        mConfig = config;
        mRecordWriteCount = 0;
        mWriter = initialise();
    }

    private SAMFileWriter initialise()
    {
        if(!mConfig.WriteTypes.contains(BAM))
            return null;

        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));
        mOutputBam = mConfig.formFilename(BAM);

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();
        fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(mOutputBam));
    }

    public void writeRecord(final SAMRecord record)
    {
        if(mWriter == null)
            return;

        ++mRecordWriteCount;
        mWriter.addAlignment(record);
    }

    public void close()
    {
        if(mWriter != null)
        {
            SV_LOGGER.info("{} records written to BAM: {}", mRecordWriteCount, mOutputBam);
            mWriter.close();
        }
    }

}
