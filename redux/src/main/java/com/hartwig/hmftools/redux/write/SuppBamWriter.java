package com.hartwig.hmftools.redux.write;

import static com.hartwig.hmftools.common.bamops.BamMerger.buildCombinedHeader;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.redux.write.FileWriterCache.BAM_FILE_ID;

import java.io.File;

import com.hartwig.hmftools.redux.ReduxConfig;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class SuppBamWriter
{
    private final String mFilename;
    private final SAMFileWriter mBamWriter;
    private int mReadCount;

    public SuppBamWriter(final String filename, final ReduxConfig config)
    {
        mFilename = filename;

        SAMFileHeader fileHeader = buildCombinedHeader(config.BamFiles, config.RefGenomeFile);

        fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        mBamWriter = new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(filename));
        mReadCount = 0;
    }

    public String filename() { return mFilename; }
    public int readCount() { return mReadCount; }

    public static String formBamFilename(final ReduxConfig config, final String multiId)
    {
        String filename = config.OutputDir + config.SampleId + "." + BAM_FILE_ID;

        if(config.OutputId != null)
            filename += "." + config.OutputId;

        filename += "." + multiId + ".supps" + BAM_EXTENSION;

        return filename;
    }

    public void writeRecord(final SAMRecord read)
    {
        mBamWriter.addAlignment(read);
        ++mReadCount;
    }

    public void close() { mBamWriter.close(); }
}
