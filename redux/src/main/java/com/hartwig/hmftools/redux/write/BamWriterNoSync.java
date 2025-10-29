package com.hartwig.hmftools.redux.write;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import com.hartwig.hmftools.redux.bqr.BaseQualRecalibration;
import com.hartwig.hmftools.redux.jitter.MsJitterAnalyser;
import com.hartwig.hmftools.redux.ReduxConfig;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

public class BamWriterNoSync extends BamWriter
{
    // writes reads in coordinate order using intelligent caching and without synchronised write calls
    private final SortedBamWriter mSortedBamWriter;

    private final SAMFileWriter mBamWriter; // when buffering and sorting are not required

    private int mUnsortedWriteCount;

    public BamWriterNoSync(
            final String filename, final ReduxConfig config, final ReadDataWriter readDataWriter, final SAMFileWriter samFileWriter,
            @Nullable final MsJitterAnalyser msJitterAnalyser, final BaseQualRecalibration bqr)
    {
        super(filename, config, readDataWriter, samFileWriter, msJitterAnalyser, bqr);

        if(config.SkipDuplicateMarking)
        {
            mBamWriter = samFileWriter;
            mSortedBamWriter = null;
        }
        else
        {
            mBamWriter = null;
            mSortedBamWriter = new SortedBamWriter(new SortedBamConfig(), samFileWriter);

            if(config.perfDebug())
                mSortedBamWriter.togglePerfDebug();
        }

        mUnsortedWriteCount = 0;
    }

    public boolean isSorted() { return true; }

    @Override
    public long unsortedWriteCount() { return mUnsortedWriteCount; }

    public void onRegionInitialised(final String chromosome, int startPosition)
    {
        if(mSortedBamWriter != null)
            mSortedBamWriter.initialiseStartPosition(chromosome, startPosition);
    }

    public void setBoundaryPosition(int position, boolean isLower)
    {
        if(mSortedBamWriter != null)
        {
            if(isLower)
                mSortedBamWriter.setUpperWritablePosition(position);
            else
                mSortedBamWriter.setUpperBoundPosition(position);
        }
    }

    public void onRegionComplete()
    {
        if(mSortedBamWriter != null)
            mSortedBamWriter.flush();
    }

    @Override
    protected void writeRecord(final SAMRecord read)
    {
        if(mSamFileWriter != null)
        {
            if(mSortedBamWriter != null)
            {
                if(mSortedBamWriter.canWriteRecord(read))
                {
                    mSortedBamWriter.addRecord(read);
                }
                else
                {
                    ++mUnsortedWriteCount;
                }
            }
            else
            {
                mBamWriter.addAlignment(read);
            }
        }
    }

    @Override
    public void close()
    {
        if(mSamFileWriter == null)
            return;

        if(mSortedBamWriter != null)
        {
            mSortedBamWriter.flush();

            RD_LOGGER.debug("sorted-writer records written({} writes={} avg={} max={}) maxCache({}) to BAM({})",
                    mSortedBamWriter.written(), mSortedBamWriter.writeCount(), mSortedBamWriter.avgWriteCount(),
                    mSortedBamWriter.maxWrite(), mSortedBamWriter.maxCache(), filenamePart(mFilename));
        }
        else if(mBamWriter != null)
        {
            mBamWriter.close();
        }

        mSamFileWriter.close();
    }
}
