package com.hartwig.hmftools.redux.write;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import com.hartwig.hmftools.redux.bqr.BaseQualRecalibration;
import com.hartwig.hmftools.redux.jitter.MsJitterAnalyser;
import com.hartwig.hmftools.redux.ReduxConfig;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

public class BamWriterSync extends BamWriter
{
    // a simple unsorted BAM writer, expected to be shared across threads so uses a synchronised write call
    private long mWriteCount;

    public BamWriterSync(
            final String filename, final ReduxConfig config, final ReadDataWriter readDataWriter, final SAMFileWriter samFileWriter,
            @Nullable final MsJitterAnalyser msJitterAnalyser, final BaseQualRecalibration bqr)
    {
        super(filename, config, readDataWriter, samFileWriter, msJitterAnalyser, bqr);
        mWriteCount = 0;
    }

    public boolean isSorted() { return false; }

    public void onRegionInitialised(final String chromosome, int startPosition) {}
    public void setBoundaryPosition(int position, boolean isLower) {}
    public void onRegionComplete() {}

    @Override
    protected void writeRecord(final SAMRecord read) { writeRecordSync(read); }

    @Override
    public long unsortedWriteCount() { return mWriteCount; }

    public synchronized void writeRecordSync(final SAMRecord read)
    {
        if(mSamFileWriter != null)
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
