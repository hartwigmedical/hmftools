package com.hartwig.hmftools.bamtools.compare;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.DELIM;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.GENOME_TERRITORY_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.MAD_COVERAGE_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.MEAN_COVERAGE_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.MEDIAN_COVERAGE_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.PCT_EXC_ADAPTER_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.PCT_EXC_TOTAL_COLUMN;
import static com.hartwig.hmftools.common.metrics.WGSMetricsFile.SD_COVERAGE_COLUMN;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bamtools.metrics.FilterType;
import com.hartwig.hmftools.bamtools.metrics.Metrics;
import com.hartwig.hmftools.bamtools.metrics.MetricsConfig;
import com.hartwig.hmftools.bamtools.metrics.Statistics;
import com.hartwig.hmftools.common.metrics.WGSMetricsFile;
import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public class ReadWriter
{
    private final BufferedWriter mWriter;

    public ReadWriter(final CompareConfig config)
    {
        mWriter = initialiseWriter(config.OutputFile);
    }

    public boolean initialised() { return mWriter != null; }

    private BufferedWriter initialiseWriter(final String filename)
    {
        try
        {
            // write summary metrics
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("ReadId,Chromosome,PosStart,MismatchType,Diff");
            writer.write(",Cigar,Flags,Paired,IsFirst,NegStrand,Duplicate,IsSupp,HasSupp");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to initialise reads comparison file: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeComparison(
            final SAMRecord read, final MismatchType mismatchType, final String diffDetails)
    {
        try
        {
            mWriter.write(format("%s,%s,%d,%s,%s",
                    read.getReadName(), read.getReferenceName(), read.getAlignmentStart(), mismatchType, diffDetails));

            mWriter.write(format(",%s,%d,%s,%s,%s,%s,%s,%s",
                    read.getCigarString(), read.getFlags(), read.getReadPairedFlag(), read.getFirstOfPairFlag(),
                    read.getReadNegativeStrandFlag(), read.getDuplicateReadFlag(), read.getSupplementaryAlignmentFlag(),
                    read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE)));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write coverage frequency file: {}", e.toString());
        }
    }

    public void close() { closeBufferedWriter(mWriter); }

    public static String readDetails(final SAMRecord read)
    {
        boolean unmapped = read.getReadUnmappedFlag();

        return format("%s_%s_%d_%s_%s_%s",
                read.getReadName(),
                unmapped ? "unmapped" : read.getReferenceName(),
                unmapped ? 0 : read.getAlignmentStart(),
                read.getReadNegativeStrandFlag() ? "fwd" : "rev",
                read.getFirstOfPairFlag() ? "R1" : "R2",
                read.getSupplementaryAlignmentFlag() ? "supp" : "prim");

    }
}
