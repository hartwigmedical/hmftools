package com.hartwig.hmftools.redux.write;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyser;
import com.hartwig.hmftools.redux.ReduxConfig;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

public class BamWriterNoSync extends BamWriter
{
    // writes reads in coordinate order using intelligent caching and without synchronised write calls
    private final SortedBamWriter mSortedBamWriter;

    private int mUnsortedWriteCount;

    public BamWriterNoSync(
            final String filename, final ReduxConfig config, final ReadDataWriter readDataWriter, final SAMFileWriter samFileWriter,
            @Nullable final JitterAnalyser jitterAnalyser)
    {
        super(filename, config, readDataWriter, samFileWriter, jitterAnalyser);

        mSortedBamWriter = new SortedBamWriter(new SortedBamConfig(), samFileWriter);

        if(config.perfDebug())
            mSortedBamWriter.togglePerfDebug();

        mUnsortedWriteCount = 0;
    }

    public boolean isSorted() { return true; }

    @Override
    public long unsortedWriteCount() { return mUnsortedWriteCount; }

    public void initialiseRegion(final String chromosome, int startPosition)
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
            if(mSortedBamWriter.canWriteRecord(read))
            {
                mSortedBamWriter.addRecord(read);
            }
            else
            {
                ++mUnsortedWriteCount;
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

        mSamFileWriter.close();
    }
}
