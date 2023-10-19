package com.hartwig.hmftools.sieve.unmap;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.io.File;
import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

// TODO(m_cooper): Duplicate code?
public class RecordWriter
{
    private final UnmapperConfig mConfig;
    private final SAMFileWriter mBamWriter;
    private int mWriteCount;

    public RecordWriter(final UnmapperConfig config)
    {
        mConfig = config;
        mBamWriter = mConfig.OutputBamFile != null ? initialiseBam() : null;
        mWriteCount = 0;
        MD_LOGGER.info("writing new BAM file: {}", mConfig.OutputBamFile);
    }

    private SAMFileWriter initialiseBam()
    {
        final SamReader samReader =
                SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenome)).open(new File(mConfig.BamFile));

        final SAMFileHeader fileHeader = samReader.getFileHeader().clone();
        fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(mConfig.OutputBamFile));
    }

    public synchronized void writeRead(final SAMRecord read)
    {
        ++mWriteCount;
        if(mBamWriter != null)
        {
            mBamWriter.addAlignment(read);
        }
    }

    public synchronized void writeReads(final List<SAMRecord> reads)
    {
        for(final SAMRecord read : reads)
        {
            ++mWriteCount;
            if(mBamWriter != null)
            {
                mBamWriter.addAlignment(read);
            }
        }
    }

    public void close()
    {
        MD_LOGGER.info("{} records written to BAM", mWriteCount);
        if(mBamWriter != null)
        {
            mBamWriter.close();
        }
    }
}
