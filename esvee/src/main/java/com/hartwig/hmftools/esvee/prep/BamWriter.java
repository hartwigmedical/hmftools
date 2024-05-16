package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.writeSortedBam;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.BAM_RECORD_SAMPLE_ID_TAG;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.BAM;
import static com.hartwig.hmftools.esvee.prep.types.WriteType.UNSORTED_BAM;

import java.io.File;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

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

    private final Map<String,SAMFileWriter> mWriters;
    private final List<String> mUnsortedBamFiles;
    private final SAMFileWriter mSingleWriter;

    public BamWriter(final PrepConfig config)
    {
        mConfig = config;
        mRecordWriteCount = 0;
        mWriters = Maps.newHashMap();
        mUnsortedBamFiles = Lists.newArrayList();

        initialiseWriters();
        mSingleWriter = mWriters.size() > 1 || mWriters.isEmpty() ? null : mWriters.get(mConfig.sampleId());
    }

    private void initialiseWriters()
    {
        if(!mConfig.WriteTypes.contains(BAM))
            return;

        for(int i = 0; i < mConfig.SampleIds.size(); ++i)
        {
            String sampleId = mConfig.SampleIds.get(i);
            String bamFile = mConfig.BamFiles.get(i);
            SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(bamFile));

            String unsortedBamFile = mConfig.formFilename(UNSORTED_BAM, sampleId);
            mUnsortedBamFiles.add(unsortedBamFile);

            SAMFileHeader fileHeader = samReader.getFileHeader().clone();
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            mWriters.put(sampleId, new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(unsortedBamFile)));
        }
    }

    public void writeRecord(final SAMRecord record)
    {
        if(mWriters.isEmpty())
            return;

        ++mRecordWriteCount;

        String sampleId = record.getStringAttribute(BAM_RECORD_SAMPLE_ID_TAG);
        record.setAttribute(BAM_RECORD_SAMPLE_ID_TAG, null); // remove since not required downstream

        if(mSingleWriter != null)
            mSingleWriter.addAlignment(record);
        else
            mWriters.get(sampleId).addAlignment(record);
    }

    public void close()
    {
        if(mWriters.isEmpty())
            return;

        SV_LOGGER.info("{} records written to BAM", mRecordWriteCount);

        for(int i = 0; i < mConfig.SampleIds.size(); ++i)
        {
            String sampleId = mConfig.SampleIds.get(i);
            mWriters.get(sampleId).close();

            String unsortedBam = mUnsortedBamFiles.get(i);
            String sortedBam = mConfig.formFilename(BAM, sampleId);

            writeSortedBam(unsortedBam, sortedBam, mConfig.BamToolPath, mConfig.Threads);
        }
    }
}
