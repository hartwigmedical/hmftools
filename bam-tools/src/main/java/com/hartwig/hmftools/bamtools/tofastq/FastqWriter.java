package com.hartwig.hmftools.bamtools.tofastq;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.tofastq.ToFastqUtils.R1;
import static com.hartwig.hmftools.bamtools.tofastq.ToFastqUtils.R2;
import static com.hartwig.hmftools.bamtools.tofastq.ToFastqUtils.UNPAIRED;
import static com.hartwig.hmftools.bamtools.tofastq.ToFastqUtils.formFilename;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.readToString;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;

import com.hartwig.hmftools.common.codon.Nucleotides;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public class FastqWriter
{
    private final String mFastqR1;
    private final String mFastqR2;

    private final String mFastqUnpaired;

    private final BufferedWriter mWriterR1;
    private final BufferedWriter mWriterR2;

    @Nullable private BufferedWriter mWriterUnpaired;

    public String getFastqR1() { return mFastqR1; }

    public String getFastqR2() { return mFastqR2; }

    @Nullable
    public String getFastqUnpaired() { return mFastqUnpaired; }

    public FastqWriter(final String filePrefix)
    {
        mFastqR1 = formFilename(filePrefix, R1);
        mWriterR1 = initialise(mFastqR1);

        mFastqR2 = formFilename(filePrefix, R2);
        mWriterR2 = initialise(mFastqR2);

        mFastqUnpaired = formFilename(filePrefix, UNPAIRED);
    }

    private BufferedWriter initialise(final String filename)
    {
        try
        {
            return createBufferedWriter(filename);
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to initialise fastq file({}): {}", filename, e);
            throw new UncheckedIOException(e);
        }
    }

    public void writeReadPair(final SAMRecord first, final SAMRecord second)
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

    public void writeUnpairedRead(final SAMRecord read)
    {
        if(mWriterUnpaired == null)
        {
            // do not create the unpaired one until needed
            mWriterUnpaired = initialise(mFastqUnpaired);
        }
        if(read.getReadPairedFlag())
        {
            BT_LOGGER.error("mate not found for paired read: {}", read);
            throw new RuntimeException("mate not found for paired read");
        }

        // write to the unpaired fastq
        writeFastqRecord(read, mWriterUnpaired);
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
        if(mWriterUnpaired != null)
            closeBufferedWriter(mWriterUnpaired);
    }
}
