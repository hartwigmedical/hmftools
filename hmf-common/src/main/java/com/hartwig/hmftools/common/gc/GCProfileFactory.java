package com.hartwig.hmftools.common.gc;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.io.reader.LineReader;

import org.jetbrains.annotations.NotNull;

public enum GCProfileFactory {
    ;

    private static final String RATIO_COLUMN_SEPARATOR = "\t";
    private static final int CHROMOSOME_COLUMN = 0;
    private static final int START_FIELD_COLUMN = 1;
    private static final int GC_CONTENT_COLUMN = 2;
    private static final int NON_N_PERCENTAGE_COLUMN = 3;
    private static final int MAPPABLE_PERCENTAGE_COLUMN = 4;

    @Deprecated
    @NotNull
    public static Multimap<String, GCProfile> loadGCContentOld(int windowSize, @NotNull final String fileName) throws IOException {
        return loadGCContentOld(windowSize, LineReader.build().readLines(new File(fileName).toPath(), x -> true));
    }

    @Deprecated
    @NotNull
    private static Multimap<String, GCProfile> loadGCContentOld(int windowSize, @NotNull final List<String> lines) {
        final Multimap<String, GCProfile> result = ArrayListMultimap.create();
        for (String line : lines) {
            final GCProfile gcProfile = fromLine(windowSize, line);
            result.put(gcProfile.chromosome(), gcProfile);
        }

        return result;
    }

    @NotNull
    public static Multimap<Chromosome, GCProfile> loadGCContent(int windowSize, @NotNull final String fileName) throws IOException {
        return loadGCContent(windowSize, LineReader.build().readLines(new File(fileName).toPath(), x -> true));
    }

    @NotNull
    private static Multimap<Chromosome, GCProfile> loadGCContent(int windowSize, @NotNull final List<String> lines) {
        final Multimap<Chromosome, GCProfile> result = ArrayListMultimap.create();
        for (String line : lines) {
            final GCProfile gcProfile = fromLine(windowSize, line);
            if (HumanChromosome.contains(gcProfile.chromosome())) {
                result.put(HumanChromosome.fromString(gcProfile.chromosome()), gcProfile);
            }
        }

        return result;
    }

    @NotNull
    @VisibleForTesting
    static GCProfile fromLine(int windowSize, @NotNull final String ratioLine) {
        final String[] values = ratioLine.split(RATIO_COLUMN_SEPARATOR);

        final String chromosome = values[CHROMOSOME_COLUMN].trim();
        final long position = Long.valueOf(values[START_FIELD_COLUMN].trim());
        final double gcContent = Double.valueOf(values[GC_CONTENT_COLUMN].trim());
        final double nonNPercentage = Double.valueOf(values[NON_N_PERCENTAGE_COLUMN].trim());
        final double mappablePercentage = Double.valueOf(values[MAPPABLE_PERCENTAGE_COLUMN].trim());

        return ImmutableGCProfile.builder()
                .chromosome(chromosome)
                .start(position + 1) // GCProfile is zero-indexed
                .end(position + windowSize)
                .gcContent(gcContent)
                .nonNPercentage(nonNPercentage)
                .mappablePercentage(mappablePercentage)
                .build();
    }
}
