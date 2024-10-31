package com.hartwig.hmftools.common.cider;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;

public class Cdr3LocusSummaryFile
{
    public enum Column
    {
        locus,
        readsUsed,
        readsTotal,
        downSampled,
        sequences,
        passSequences
    }

    private static final String FILE_EXTENSION = ".cider.locus_stats.tsv";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static List<Cdr3LocusSummary> read(final String filename)
    {
        try (DelimFileReader reader = new DelimFileReader(filename))
        {
            return reader.stream().map(row -> ImmutableCdr3LocusSummary.builder()
                    .locus(row.get(Column.locus))
                    .readsUsed(row.getInt(Column.readsUsed))
                    .readsTotal(row.getInt(Column.readsTotal))
                    .downSampled(Boolean.parseBoolean(row.get(Column.downSampled)))
                    .sequences(row.getInt(Column.sequences))
                    .passSequences(row.getInt(Column.passSequences))
                .build()).collect(Collectors.toList());
        }
    }
}
