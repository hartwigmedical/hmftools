package com.hartwig.hmftools.bamtools.compare;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

import java.util.List;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.commons.lang3.tuple.Triple;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public class ReadWriter implements AutoCloseable
{
    enum Column
    {
        readId, chromosome, posStart, mismatchType, diff, mateChr, matePos,
        cigar, flags, mapQual, paired, isFirst, negStrand, duplicate, isSupp, suppData
    }

    private final DelimFileWriter<Triple<SAMRecord, MismatchType, List<String>>> mWriter;

    public ReadWriter(final CompareConfig config)
    {
        mWriter = new DelimFileWriter<>(config.OutputFile, Column.values(),
            (t, row) ->
            {
                SAMRecord read = t.getLeft();
                @Nullable List<String> diffList = t.getRight();
                String diffDetails = diffList != null ? String.join(";", diffList) : "";
                row.set(Column.readId, read.getReadName());
                row.set(Column.chromosome, read.getReferenceName());
                row.set(Column.posStart, read.getAlignmentStart());
                row.set(Column.mismatchType, t.getMiddle().name());
                row.set(Column.diff, diffDetails);
                row.set(Column.mateChr, read.getMateReferenceName());
                row.set(Column.matePos, read.getMateAlignmentStart());
                row.set(Column.cigar, read.getCigarString());
                row.set(Column.flags, read.getFlags());
                row.set(Column.mapQual, read.getMappingQuality());
                row.set(Column.paired, Boolean.toString(read.getReadPairedFlag()));
                row.set(Column.isFirst, Boolean.toString(read.getFirstOfPairFlag()));
                row.set(Column.negStrand, Boolean.toString(read.getReadNegativeStrandFlag()));
                row.set(Column.duplicate, Boolean.toString(read.getDuplicateReadFlag()));
                row.set(Column.isSupp, Boolean.toString(read.getSupplementaryAlignmentFlag()));
                row.set(Column.suppData, read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE) ? SupplementaryReadData.extractAlignment(read).asCsv() : "N/A");
            });
    }

    public boolean initialised() { return mWriter != null; }

    public void writeComparison(final SAMRecord read, final MismatchType mismatchType, @Nullable final List<String> diffList)
    {
        mWriter.writeRow(Triple.of(read, mismatchType, diffList));
    }

    @Override
    public void close() { mWriter.close(); }

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
