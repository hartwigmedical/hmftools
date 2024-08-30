package com.hartwig.hmftools.common.bamops;

import static com.hartwig.hmftools.common.bamops.BamMerger.UNMAPPED_READS;
import static com.hartwig.hmftools.common.bamops.BamMerger.buildCombinedHeader;
import static com.hartwig.hmftools.common.bamops.BamMerger.formBamFilename;

import java.io.File;
import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class UnmappedMergeTask extends Thread
{
    private final List<String> mInputBams;
    private final String mOutputBamPrefix;

    private final String mRefGenomeFile;

    public UnmappedMergeTask(final List<String> inputBams, final String refGenomeFile, final String outputBamPrefix)
    {
        mOutputBamPrefix = outputBamPrefix;
        mInputBams = inputBams;
        mRefGenomeFile = refGenomeFile;

        start();
    }

    public void run()
    {
        String unmappedBam = formBamFilename(mOutputBamPrefix, UNMAPPED_READS);

        SAMFileHeader fileHeader = buildCombinedHeader(mInputBams, mRefGenomeFile);

        SAMFileWriter bamWriter = new SAMFileWriterFactory().makeBAMWriter(fileHeader, true, new File(unmappedBam));

        for(String inputBam : mInputBams)
        {
            SamReader bamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mRefGenomeFile)).open(new File(inputBam));

            SAMRecordIterator iterator = bamReader.queryUnmapped();
            while(iterator.hasNext())
            {
                bamWriter.addAlignment(iterator.next());
            }
        }

        bamWriter.close();
    }
}
