package com.hartwig.hmftools.bamtools.btofmc.writer;

import static com.google.common.io.CharStreams.nullWriter;
import static com.hartwig.hmftools.bamtools.btofmc.BamToFastqConfig.BFQ_LOGGER;
import static com.hartwig.hmftools.bamtools.btofmc.BindingAnnotations.Fastq1;
import static com.hartwig.hmftools.bamtools.btofmc.BindingAnnotations.Fastq2;
import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;

import java.io.BufferedWriter;
import java.io.IOException;

import javax.annotation.Nullable;

import com.google.common.io.CharSink;
import com.google.common.io.Closer;
import com.google.inject.Inject;
import com.hartwig.hmftools.bamtools.btofmc.BamToFastqConfig;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import org.apache.commons.lang3.tuple.Pair;

import htsjdk.samtools.SAMRecord;

// TODO NEXT: TEST
public class NoSyncPairedFastqWriter extends PairedFastqWriterInterface
{
    private final BamToFastqConfig mConfig;
    private final BufferedWriter mFastq1Writer;
    private final BufferedWriter mFastq2Writer;
    private final Closer mCloser;

    private final PerformanceCounter mWritePc;

    private int mPairedRecordCount;

    @Inject
    public NoSyncPairedFastqWriter(final BamToFastqConfig config, @Fastq1 @Nullable final CharSink fastq1,
            @Fastq2 @Nullable final CharSink fastq2)
    {
        mConfig = config;

        mCloser = Closer.create();
        Pair<BufferedWriter, BufferedWriter> fastqWriters = openFastqWriters(fastq1, fastq2);
        mFastq1Writer = fastqWriters.getLeft();
        mFastq2Writer = fastqWriters.getRight();

        mWritePc = new PerformanceCounter("Total write", false);

        mPairedRecordCount = 0;
    }

    private Pair<BufferedWriter, BufferedWriter> openFastqWriters(@Nullable final CharSink fastq1, @Nullable final CharSink fastq2)
    {
        try
        {
            BufferedWriter fastq1Writer = mCloser.register(new BufferedWriter(fastq1 == null ? nullWriter() : fastq1.openStream()));
            BufferedWriter fastq2Writer = mCloser.register(new BufferedWriter(fastq2 == null ? nullWriter() : fastq2.openStream()));
            return Pair.of(fastq1Writer, fastq2Writer);
        }
        catch(IOException e)
        {
            try
            {
                mCloser.close();
            }
            catch(IOException e2)
            {
            }

            throw new RuntimeException("Failed to open writers to fastq output files", e);
        }
    }

    @Override
    public void writePairedFastqRecord(final SAMRecord read1, final SAMRecord read2)
    {
        mWritePc.start();
        writeFastqRecord(read1);
        writeFastqRecord(read2);
        mWritePc.stop();

        ++mPairedRecordCount;
    }

    private void writeFastqRecord(final SAMRecord read)
    {
        try
        {
            BufferedWriter writer = read.getFirstOfPairFlag() ? mFastq1Writer : mFastq2Writer;

            writer.write('@');
            writer.write(read.getReadName());
            writer.newLine();

            String readBases = read.getReadString();
            if(read.getReadNegativeStrandFlag())
            {
                StringBuilder readBasesBuilder = new StringBuilder(readBases);

                int left = 0;
                int right = readBasesBuilder.length() - 1;
                while(left < right)
                {
                    char leftCh = swapDnaBase(readBasesBuilder.charAt(right));
                    char rightCh = swapDnaBase(readBasesBuilder.charAt(left));
                    readBasesBuilder.setCharAt(left, leftCh);
                    readBasesBuilder.setCharAt(right, rightCh);
                    ++left;
                    --right;
                }

                if(left == right)
                {
                    readBasesBuilder.setCharAt(left, swapDnaBase((readBasesBuilder.charAt(left))));
                }

                readBases = readBasesBuilder.toString();
            }

            writer.write(readBases);
            writer.newLine();

            writer.write('+');
            writer.newLine();

            String baseQuals = read.getBaseQualityString();
            if(read.getReadNegativeStrandFlag())
            {
                StringBuilder baseQualsBuilder = new StringBuilder(baseQuals);
                baseQualsBuilder.reverse();
                baseQuals = baseQualsBuilder.toString();
            }

            writer.write(baseQuals);
            writer.newLine();
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to write fastq record to fastq output", e);
        }
    }

    public void mergeStats(final NoSyncPairedFastqWriter other)
    {
        mWritePc.merge(other.mWritePc);
        mPairedRecordCount += other.mPairedRecordCount;
    }

    public void logStats()
    {
        BFQ_LOGGER.info("NoSyncPairedFastqWriter stats: pairedRecordCount({})", mPairedRecordCount);
        if(mConfig.PerfDebug)
        {
            mWritePc.logStats();
        }
    }

    @Override
    public void close() throws IOException
    {
        mCloser.close();
    }
}
