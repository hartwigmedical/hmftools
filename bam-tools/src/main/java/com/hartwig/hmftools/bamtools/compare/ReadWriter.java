package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

public class ReadWriter implements AutoCloseable
{
    static final String[] MATE_DIFF_COLUMN_NAMES = { "mateChr", "matePos", "mateNegStrand" };

    private static final String[] HEADERS = {
            "ReadId", "Chromosome", "PosStart", "MismatchType", "DiffBucket", "Diff",
            "MateChr", "MatePos", "Cigar", "Flags", "MapQual", "IsFirst",
            "NegStrand", "Duplicate", "Supplementary", "SuppData", "Secondary",
            "MateChrDiff", "MatePosDiff", "MateNegStrandDiff"
    };

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
            writer.write(String.join(TSV_DELIM, HEADERS));
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to initialise BAM comparison file: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeComparison(
            final SAMRecord read, final MismatchType type, final List<String> diffs, final DiffBucket bucket)
    {
        Map<String, String> mateColumnValues = new HashMap<>();
        String packedDiff = packDiff(diffs, mateColumnValues);
        String suppData = read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE)
                ? SupplementaryReadData.extractAlignment(read).asDelimStr() : "N/A";
        boolean isFirst = !read.getReadPairedFlag() || read.getFirstOfPairFlag();

        String[] row = {
                read.getReadName(),
                read.getReferenceName(),
                String.valueOf(read.getAlignmentStart()),
                String.valueOf(type),
                bucket != null ? bucket.name() : "",
                packedDiff,
                read.getMateReferenceName(),
                String.valueOf(read.getMateAlignmentStart()),
                read.getCigarString(),
                String.valueOf(read.getFlags()),
                String.valueOf(read.getMappingQuality()),
                String.valueOf(isFirst),
                String.valueOf(read.getReadNegativeStrandFlag()),
                String.valueOf(read.getDuplicateReadFlag()),
                String.valueOf(read.getSupplementaryAlignmentFlag()),
                suppData,
                String.valueOf(read.isSecondaryAlignment()),
                mateColumnValues.getOrDefault("mateChr", ""),
                mateColumnValues.getOrDefault("matePos", ""),
                mateColumnValues.getOrDefault("mateNegStrand", "")
        };

        try
        {
            mWriter.write(String.join(TSV_DELIM, row));
            mWriter.newLine();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to write BAM comparison file: {}", e.toString());
        }
    }

    private static String packDiff(final List<String> diffs, final Map<String, String> mateColumnValues)
    {
        if(diffs == null)
            return "";

        StringJoiner rest = new StringJoiner(ITEM_DELIM);
        for(String entry : diffs)
        {
            int paren = entry.indexOf('(');
            String name = paren > 0 ? entry.substring(0, paren) : null;
            if(name != null && isMateColumn(name) && entry.endsWith(")"))
                mateColumnValues.put(name, entry.substring(paren + 1, entry.length() - 1));
            else
                rest.add(entry);
        }
        return rest.toString();
    }

    private static boolean isMateColumn(final String name)
    {
        for(String col : MATE_DIFF_COLUMN_NAMES)
        {
            if(col.equals(name))
                return true;
        }
        return false;
    }

    @Override
    public void close() { closeBufferedWriter(mWriter); }
}
