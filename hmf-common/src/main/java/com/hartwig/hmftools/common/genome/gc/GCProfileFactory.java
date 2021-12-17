package com.hartwig.hmftools.common.genome.gc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

public final class GCProfileFactory {

    private static final String RATIO_COLUMN_SEPARATOR = "\t";
    private static final int CHROMOSOME_COLUMN = 0;
    private static final int START_FIELD_COLUMN = 1;
    private static final int GC_CONTENT_COLUMN = 2;
    private static final int NON_N_PERCENTAGE_COLUMN = 3;
    private static final int MAPPABLE_PERCENTAGE_COLUMN = 4;

    private GCProfileFactory() {
    }

    @NotNull
    public static Multimap<Chromosome, GCProfile> loadGCContent(int windowSize, @NotNull final String fileName) throws IOException {
        return loadGCContent(windowSize, Files.readAllLines(new File(fileName).toPath()));
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
    private static GCProfile fromLine(int windowSize, @NotNull final String ratioLine) {
        final String[] values = ratioLine.split(RATIO_COLUMN_SEPARATOR);

        final String chromosome = values[CHROMOSOME_COLUMN].trim();
        final int position = Integer.parseInt(values[START_FIELD_COLUMN].trim());
        final double gcContent = Double.parseDouble(values[GC_CONTENT_COLUMN].trim());
        final double nonNPercentage = Double.parseDouble(values[NON_N_PERCENTAGE_COLUMN].trim());
        final double mappablePercentage = Double.parseDouble(values[MAPPABLE_PERCENTAGE_COLUMN].trim());

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
