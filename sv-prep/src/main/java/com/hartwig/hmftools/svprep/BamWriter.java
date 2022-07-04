package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.svprep.WriteType.BAM;

import java.io.File;
import java.util.List;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamWriter
{
    private final SvConfig mConfig;

    private SAMFileWriter mWriter;

    public BamWriter(final SvConfig config)
    {
        mConfig = config;
        mWriter = initialise();
    }

    private SAMFileWriter initialise()
    {
        if(!mConfig.WriteTypes.contains(BAM))
            return null;

        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));
        String outputBam = mConfig.formFilename(BAM);

        return new SAMFileWriterFactory().makeBAMWriter(samReader.getFileHeader(), false, new File(outputBam));
    }

    public synchronized void writeRecords(final List<ReadRecord> reads)
    {
        if(mWriter == null)
            return;

        reads.forEach(x -> mWriter.addAlignment(x.record()));
    }

    public void close()
    {
        if(mWriter != null)
            mWriter.close();
    }

}
