package com.hartwig.hmftools.bamtools.tofastq;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.readToString;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;

import htsjdk.samtools.SAMRecord;

public class FastqWriter
{
    private final String mFilePrefix;

    private final BufferedWriter mWriterR1;
    private final BufferedWriter mWriterR2;

    private final List<SAMRecord> mUnpairedReads;

    private static final String ZIPPED_SUFFIX = ".fastq.gz";
    private static final String UNZIPPED_SUFFIX = ".fastq";

    public FastqWriter(final String filePrefix, boolean writeUnzipped)
    {
        mFilePrefix = filePrefix;
        mUnpairedReads = new ArrayList<>();

        String fastqR1 = formFilename(mFilePrefix, true, writeUnzipped);
        mWriterR1 = initialise(fastqR1);

        String fastqR2 = formFilename(mFilePrefix, false, writeUnzipped);
        mWriterR2 = initialise(fastqR2);
    }

    private String formFilename(final String filePrefix, final boolean isR1, boolean writeUnzipped)
    {
        String filename = filePrefix + ".";
        filename += isR1 ? "R1" : "R2";
        filename += writeUnzipped ? UNZIPPED_SUFFIX : ZIPPED_SUFFIX;
        return filename;
    }

    private BufferedWriter initialise(final String filename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename);

            return writer;
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to initialise fastq file({}): {}", filename, e.toString());
            System.exit(1);
        }

        return null;
    }

    public synchronized void processReadPairSync(final SAMRecord first, final SAMRecord second)
    {
        processReadPair(first, second);
    }

    public void processReadPairNoSync(final SAMRecord first, final SAMRecord second)
    {
        processReadPair(first, second);
    }

    private void processReadPair(final SAMRecord first, final SAMRecord second)
    {
        if(first.getFirstOfPairFlag())
        {
            writeFastqRecord(first, mWriterR1);
            writeFastqRecord(second, mWriterR2);
        }
        else
        {
            writeFastqRecord(first, mWriterR2);
            writeFastqRecord(second, mWriterR1);
        }
    }

    public synchronized void processUnpairedRead(final SAMRecord read)
    {
        // cache until all paired reads have been written
        mUnpairedReads.add(read);
    }

    public void writeUnpairedRead(final SAMRecord read)
    {
        writeFastqRecord(read, mWriterR1);
    }

    public void writeUnpairedReads()
    {
        if(!mUnpairedReads.isEmpty())
        {
            for(SAMRecord read : mUnpairedReads)
            {
                if(read.getReadPairedFlag())
                {
                    BT_LOGGER.error("mate not found for paired read: {}", read);
                }

                // write these to the first fastq file only
                writeFastqRecord(read, mWriterR1);
            }

            mUnpairedReads.clear();
        }
    }

    private void writeFastqRecord(final SAMRecord read, BufferedWriter writer)
    {
        try
        {
            writer.write('@');
            writer.write(read.getReadName());
            writer.write('\n'); // must use this instead of newline, otherwise would write \r\n in windows

            String readBases = read.getReadString();
            if(read.getReadNegativeStrandFlag())
            {
                readBases = Nucleotides.reverseComplementBases(readBases);
            }

            writer.write(readBases);
            writer.write('\n');

            writer.write('+');
            writer.write('\n');

            String baseQuals = read.getBaseQualityString();
            if(read.getReadNegativeStrandFlag())
            {
                StringBuilder baseQualsBuilder = new StringBuilder(baseQuals);
                baseQualsBuilder.reverse();
                baseQuals = baseQualsBuilder.toString();
            }

            writer.write(baseQuals);
            writer.write('\n');
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write read({}): {}", readToString(read), e.toString());
            throw new UncheckedIOException(e);
        }
    }

    public void close()
    {
        closeBufferedWriter(mWriterR1);
        closeBufferedWriter(mWriterR2);
    }
}
