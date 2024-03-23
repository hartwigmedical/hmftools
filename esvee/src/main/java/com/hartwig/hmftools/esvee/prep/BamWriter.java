package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.CommonUtils.writeSortedBam;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.BAM;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.UNSORTED_BAM;

import java.io.File;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamWriter
{
    private final PrepConfig mConfig;

    private int mRecordWriteCount;
    private final SAMFileWriter mWriter;
    private String mUnsortedBam;

    public BamWriter(final PrepConfig config)
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
        mUnsortedBam = mConfig.formFilename(UNSORTED_BAM);

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();
        fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(mUnsortedBam));
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
        if(mWriter == null)
            return;

        SV_LOGGER.info("{} records written to BAM: {}", mRecordWriteCount, mUnsortedBam);
        mWriter.close();

        String sortedBam = mConfig.formFilename(BAM);
        writeSortedBam(mUnsortedBam, sortedBam, mConfig.BamToolPath, mConfig.Threads);
    }
}
