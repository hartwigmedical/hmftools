package com.hartwig.hmftools.markdups.write;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import com.hartwig.hmftools.markdups.MarkDupsConfig;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

public class BamWriterNoSync extends BamWriter
{
    private final SortedBamWriter mSortedBamWriter;
    private final BamWriter mSharedUnsortedWriter;
    private boolean mWriteSorted;

    public BamWriterNoSync(
            final String filename, final MarkDupsConfig config, final ReadDataWriter readDataWriter, final SAMFileWriter samFileWriter,
            boolean writeSorted, final BamWriter sharedUnsortedWriter)
    {
        super(filename, config, readDataWriter, samFileWriter);

        mWriteSorted = writeSorted;

        if(mWriteSorted)
        {
            mSortedBamWriter = new SortedBamWriter(new SortedBamConfig(), samFileWriter);

            if(config.PerfDebug)
                mSortedBamWriter.togglePerfDebug();

            mSharedUnsortedWriter = sharedUnsortedWriter;
        }
        else
        {
            mSortedBamWriter = null;
            mSharedUnsortedWriter = null;
        }
    }

    public boolean isSorted() { return mWriteSorted; }

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
            if(mWriteSorted)
            {
                if(mSortedBamWriter.canWriteRecord(read))
                {
                    mSortedBamWriter.addRecord(read);
                }
                else
                {
                    mSharedUnsortedWriter.writeRecord(read);
                }
            }
            else
            {
                mSamFileWriter.addAlignment(read);
            }
        }
    }

    @Override
    public void close()
    {
        if(mSortedBamWriter != null)
        {
            mSortedBamWriter.flush();

            MD_LOGGER.debug("sorted-writer records written({} writes={} avg={} max={}) maxCache({}) to BAM({})",
                    mSortedBamWriter.written(), mSortedBamWriter.writeCount(), mSortedBamWriter.avgWriteCount(),
                    mSortedBamWriter.maxWrite(), mSortedBamWriter.maxCache(), filenamePart(mFilename));
        }

        if(mSamFileWriter != null)
            mSamFileWriter.close();
    }
}
