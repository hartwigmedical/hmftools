package com.hartwig.hmftools.common.cider;

import java.io.File;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;

public class Cdr3SequenceFile
{
    enum Column
    {
        cdr3Seq,
        cdr3AA,
        locus,
        filter,
        blastnStatus,
        minHighQualBaseReads,
        assignedReads,
        inFrame,
        containsStop
    }

    private static final String FILE_EXTENSION = ".cider.vdj.tsv.gz";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static List<Cdr3Sequence> read(final String filename)
    {
        try (DelimFileReader reader = new DelimFileReader(filename))
        {
            return reader.stream().map(row -> ImmutableCdr3Sequence.builder()
                        .cdr3Seq(Objects.requireNonNull(row.get(Column.cdr3Seq)))
                        .cdr3AA(Objects.requireNonNull(row.get(Column.cdr3AA)))
                        .locus(Objects.requireNonNull(row.get(Column.locus)))
                        .filter(Objects.requireNonNull(row.get(Column.filter)))
                        .blastnStatus(Objects.requireNonNull(row.get(Column.blastnStatus)))
                        .minHighQualBaseReads(Objects.requireNonNull(row.getInt(Column.minHighQualBaseReads)))
                        .assignedReads(Objects.requireNonNull(row.getInt(Column.assignedReads)))
                        .inFrame(Boolean.parseBoolean(Objects.requireNonNull(row.get(Column.inFrame))))
                        .containsStop(Boolean.parseBoolean(Objects.requireNonNull(row.get(Column.containsStop))))
                    .build()).collect(Collectors.toList());
        }
    }
}
