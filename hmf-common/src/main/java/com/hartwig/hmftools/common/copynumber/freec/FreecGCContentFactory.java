package com.hartwig.hmftools.common.copynumber.freec;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.io.path.PathRegexFinder;
import com.hartwig.hmftools.common.io.reader.LineReader;

import org.jetbrains.annotations.NotNull;

public enum FreecGCContentFactory {
    ;
    private static final String GC_FILE = "GC_profile.1000bp.cnp";

    private static final String RATIO_COLUMN_SEPARATOR = "\t";
    private static final int CHROMOSOME_COLUMN = 0;
    private static final int START_FIELD_COLUMN = 1;
    private static final int GC_CONTENT_COLUMN = 2;
    private static final int NON_N_PERCENTAGE_COLUMN = 3;
    private static final int MAPPABLE_PERCENTAGE_COLUMN = 4;

    @NotNull
    public static Multimap<String, FreecGCContent> loadGCContent(@NotNull final String basePath) throws IOException, HartwigException {
        final Path path = PathRegexFinder.build().findPath(basePath, GC_FILE);
        return loadGCContent(LineReader.build().readLines(path, x -> true));
    }

    @NotNull
    private static Multimap<String, FreecGCContent> loadGCContent(@NotNull final List<String> lines) {

        final Multimap<String, FreecGCContent> result = ArrayListMultimap.create();
        for (String line : lines) {
            final FreecGCContent gcContent = fromLine(line);
            result.put(gcContent.chromosome(), gcContent);
        }

        return result;
    }

    @NotNull
    @VisibleForTesting
    static FreecGCContent fromLine(@NotNull final String ratioLine) {
        final String[] values = ratioLine.split(RATIO_COLUMN_SEPARATOR);

        final String chromosome = values[CHROMOSOME_COLUMN].trim();
        final long position = Long.valueOf(values[START_FIELD_COLUMN].trim());
        final double gcContent = Double.valueOf(values[GC_CONTENT_COLUMN].trim());
        final double nonNPercentage = Double.valueOf(values[NON_N_PERCENTAGE_COLUMN].trim());
        final double mappablePercentage = Double.valueOf(values[MAPPABLE_PERCENTAGE_COLUMN].trim());

        return ImmutableFreecGCContent.builder()
                .chromosome(chromosome)
                .position(position)
                .gcContent(gcContent)
                .nonNPercentage(nonNPercentage)
                .mappablePercentage(mappablePercentage)
                .build();
    }
}
