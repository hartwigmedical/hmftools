package com.hartwig.hmftools.bamtools.biomodalcollapse;

import java.io.BufferedWriter;
import java.io.IOException;

import htsjdk.samtools.fastq.FastqRecord;

// TODO: Remove
public class SynchronizedPairedFastqWriter
{
    private final BufferedWriter mFastq1Writer;
    private final BufferedWriter mFastq2Writer;

    public SynchronizedPairedFastqWriter(final BufferedWriter fastq1Writer, final BufferedWriter fastq2Writer)
    {
        mFastq1Writer = fastq1Writer;
        mFastq2Writer = fastq2Writer;
    }

    public synchronized void write(final FastqRecord fastq1, final FastqRecord fastq2)
    {
        try
        {
            mFastq1Writer.write('@');
            mFastq1Writer.write(fastq1.getReadName());
            mFastq1Writer.newLine();
            mFastq1Writer.write(fastq1.toString());
            mFastq1Writer.newLine();
            mFastq1Writer.write(fastq1.getBaseQualityHeader());
            mFastq1Writer.newLine();
            mFastq1Writer.write(fastq1.getBaseQualityString());
            mFastq1Writer.newLine();

            mFastq2Writer.write('@');
            mFastq2Writer.write(fastq2.getReadName());
            mFastq2Writer.newLine();
            mFastq2Writer.write(fastq2.toString());
            mFastq2Writer.newLine();
            mFastq2Writer.write(fastq2.getBaseQualityHeader());
            mFastq2Writer.newLine();
            mFastq2Writer.write(fastq2.getBaseQualityString());
            mFastq2Writer.newLine();
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }
}
