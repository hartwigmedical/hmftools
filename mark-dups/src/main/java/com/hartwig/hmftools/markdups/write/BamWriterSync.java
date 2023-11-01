package com.hartwig.hmftools.markdups.write;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import com.hartwig.hmftools.markdups.MarkDupsConfig;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

public class BamWriterSync extends BamWriter
{
    private long mWriteCount;

    public BamWriterSync(
            final String filename, final MarkDupsConfig config, final ReadDataWriter readDataWriter, final SAMFileWriter samFileWriter)
    {
        super(filename, config, readDataWriter, samFileWriter);
        mWriteCount = 0;
    }

    public boolean isSorted() { return false; }

    public void initialiseRegion(final String chromosome, int startPosition) {}
    public void setBoundaryPosition(int position, boolean isLower) {}
    public void onRegionComplete() {}

    @Override
    protected void writeRecord(final SAMRecord read) { writeRecordSync(read); }

    public synchronized void writeRecordSync(final SAMRecord read)
    {
        mSamFileWriter.addAlignment(read);
        ++mWriteCount;
    }

    @Override
    public void close()
    {
        if(mSamFileWriter != null)
        {
            MD_LOGGER.debug("unsorted-writer records written({}) to BAM({})", mWriteCount, filenamePart(mFilename));
            mSamFileWriter.close();
        }
    }
}
