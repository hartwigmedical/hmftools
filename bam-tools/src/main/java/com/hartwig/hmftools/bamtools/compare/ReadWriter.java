package com.hartwig.hmftools.bamtools.compare;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.samtools.SupplementaryReadData;

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

            writer.write("ReadId\tChromosome\tPosStart\tMismatchType\tDiff\tMateChr\tMatePos");
            writer.write("\tCigar\tFlags\tMapQual\tPaired\tIsFirst\tNegStrand\tDuplicate\tIsSupp\tSuppData");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to initialise BAM comparison file: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeComparison(final SAMRecord read, final MismatchType mismatchType, final List<String> diffList)
    {
        try
        {
            String diffDetails = diffList != null ? diffList.stream().collect(Collectors.joining(";")) : "";
            mWriter.write(format("%s\t%s\t%d\t%s\t%s\t%s\t%d",
                    read.getReadName(), read.getReferenceName(), read.getAlignmentStart(), mismatchType, diffDetails,
                    read.getMateReferenceName(), read.getMateAlignmentStart()));

            mWriter.write(format("\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s",
                    read.getCigarString(), read.getFlags(), read.getMappingQuality(), read.getReadPairedFlag(), read.getFirstOfPairFlag(),
                    read.getReadNegativeStrandFlag(), read.getDuplicateReadFlag(), read.getSupplementaryAlignmentFlag(),
                    read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE) ? SupplementaryReadData.extractAlignment(read).asCsv() : "N/A"));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write BAM comparison file: {}", e.toString());
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
