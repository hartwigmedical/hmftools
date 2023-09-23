package com.hartwig.hmftools.markdups;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.io.File;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class FileWriterCache
{
    private final MarkDupsConfig mConfig;
    private final ReadDataWriter mReadDataWriter;

    private final BamWriter mSharedBamWriter;
    private final List<BamWriter> mBamWriters;

    private static final String BAM_FILE_ID = "mark_dups";

    public FileWriterCache(final MarkDupsConfig config)
    {
        mConfig = config;

        mReadDataWriter = new ReadDataWriter(mConfig);

        mBamWriters = Lists.newArrayList();
        mSharedBamWriter = mConfig.MultiBam ? null : createBamWriter(null);
    }

    public BamWriter getBamWriter(final String chromosome)
    {
        if(mSharedBamWriter != null)
            return mSharedBamWriter;

        return createBamWriter(chromosome);
    }

    public void close()
    {
        mReadDataWriter.close();
        mBamWriters.forEach(x -> x.close());
    }

    public BamWriter createBamWriter(@Nullable final String multiId)
    {
        SAMFileWriter samFileWriter = null;

        if(mConfig.WriteBam)
        {
            String filename = formBamFilename(multiId);

            if(multiId == null)
            {
                MD_LOGGER.info("writing BAM file: {}", filename);

            }
            else
            {
                MD_LOGGER.debug("writing tmp BAM file: {}", filename);
            }

            samFileWriter = initialiseSamFileWriter(filename);
        }

        BamWriter bamWriter = new BamWriter(mConfig, mReadDataWriter, samFileWriter);
        mBamWriters.add(bamWriter);
        return bamWriter;
    }

    private SAMFileWriter initialiseSamFileWriter(final String filename)
    {
        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();

        if(mConfig.SortedBam)
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        else
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(filename));
    }

    private String formBamFilename(@Nullable final String multiId)
    {
        String filename = mConfig.OutputDir + mConfig.SampleId + "." + BAM_FILE_ID;

        if(mConfig.OutputId != null)
            filename += "." + mConfig.OutputId;

        if(multiId != null)
            filename += "." + multiId;

        filename += ".bam";
        return filename;
    }

}
