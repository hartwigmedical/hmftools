package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

public class ReadWriter implements AutoCloseable
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
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add("ReadId").add("Chromosome").add("PosStart").add("MismatchType").add("Diff").add("MateChr").add("MatePos");
            sj.add("Cigar").add("Flags").add("MapQual").add("IsFirst").add("NegStrand").add("Duplicate").add("IsSupp").add("SuppData");
            writer.write(sj.toString());
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
            String diffDetails = diffList != null ? diffList.stream().collect(Collectors.joining(ITEM_DELIM)) : "";

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(read.getReadName());
            sj.add(read.getReferenceName());
            sj.add(String.valueOf(read.getAlignmentStart()));
            sj.add(String.valueOf(mismatchType));
            sj.add(diffDetails);
            sj.add(read.getMateReferenceName());
            sj.add(String.valueOf(read.getMateAlignmentStart()));
            sj.add(read.getCigarString());
            sj.add(String.valueOf(read.getFlags()));
            sj.add(String.valueOf(read.getMappingQuality()));
            sj.add(String.valueOf(!read.getReadPairedFlag() || read.getFirstOfPairFlag()));
            sj.add(String.valueOf(read.getReadNegativeStrandFlag()));
            sj.add(String.valueOf(read.getDuplicateReadFlag()));
            sj.add(String.valueOf(read.getSupplementaryAlignmentFlag()));
            sj.add(read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE) ? SupplementaryReadData.extractAlignment(read).asCsv() : "N/A");
            mWriter.write(sj.toString());
            mWriter.newLine();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write BAM comparison file: {}", e.toString());
        }
    }

    @Override
    public void close() { closeBufferedWriter(mWriter); }
}
