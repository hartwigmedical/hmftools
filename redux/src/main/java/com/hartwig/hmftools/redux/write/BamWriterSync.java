package com.hartwig.hmftools.redux.write;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyser;
import com.hartwig.hmftools.redux.ReduxConfig;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

public class BamWriterSync extends BamWriter
{
    private long mWriteCount;

    public BamWriterSync(
            final String filename, final ReduxConfig config, final ReadDataWriter readDataWriter, final SAMFileWriter samFileWriter,
            @Nullable final JitterAnalyser jitterAnalyser)
    {
        super(filename, config, readDataWriter, samFileWriter, jitterAnalyser);
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
            RD_LOGGER.debug("unsorted-writer records written({}) to BAM({})", mWriteCount, filenamePart(mFilename));
            mSamFileWriter.close();
        }
    }
}
