package com.hartwig.hmftools.common.cider;

import java.io.File;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;

public class Cdr3SequenceFile
{
    public enum Column
    {
        cdr3Seq,
        cdr3AA,
        locus,
        filter,
        alignmentStatus,
        minHighQualBaseReads,
        assignedReads,
        inFrame,
        containsStop
    }

    private static final String OLD_BLASTN_STATUS_COLUMN_NAME = "blastnStatus";
    private static final String FILE_EXTENSION = ".cider.vdj.tsv.gz";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static List<Cdr3Sequence> read(final String filename)
    {
        try (DelimFileReader reader = new DelimFileReader(filename))
        {
            return reader.stream().map(row -> readRow(row)).collect(Collectors.toList());
        }
    }

    private static ImmutableCdr3Sequence readRow(final DelimFileReader.Row row)
    {
        String blastnStatus;
        try
        {
            // Try old name first such that error message refers to new name if both fail
            blastnStatus = row.get(OLD_BLASTN_STATUS_COLUMN_NAME);
        }
        catch(NoSuchElementException exception)
        {
            blastnStatus = row.get(Column.alignmentStatus);
        }

        return ImmutableCdr3Sequence.builder()
                .cdr3Seq(row.get(Column.cdr3Seq))
                .cdr3AA(row.get(Column.cdr3AA))
                .locus(row.get(Column.locus))
                .filter(row.get(Column.filter))
                .blastnStatus(blastnStatus)
                .minHighQualBaseReads(row.getInt(Column.minHighQualBaseReads))
                .assignedReads(row.getInt(Column.assignedReads))
                .inFrame(Boolean.parseBoolean(row.get(Column.inFrame)))
                .containsStop(Boolean.parseBoolean(row.get(Column.containsStop)))
                .build();
    }
}
