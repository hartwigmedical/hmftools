package com.hartwig.hmftools.svvisualise.data;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class Exons {

    private static final String HEADER = "SampleId";
    private static final String COMMENT = "#";
    private static final String DELIMITER = ",";

    @NotNull
    public static List<Exon> readExons(@NotNull final String fileName) throws IOException {
        return fromString(Files.readAllLines(new File(fileName).toPath()));
    }

    @VisibleForTesting
    @NotNull
    static List<Exon> fromString(@NotNull final List<String> lines) {
        final List<Exon> result = Lists.newArrayList();

        for (final String line : lines) {

            if (!line.startsWith(COMMENT) && !line.startsWith(HEADER)) {
                String[] values = line.split(DELIMITER);

                final Exon newExon = ImmutableExon.builder()
                        .sampleId(values[0])
                        .clusterId(Integer.valueOf(values[1]))
                        .gene(values[2])
                        .chromosome(values[4])
                        .rank(Integer.valueOf(values[6]))
                        .start(Long.valueOf(values[7]))
                        .end(Long.valueOf(values[8]))
                        .build();

                result.add(newExon);

            }
        }

        return result;
    }

}
