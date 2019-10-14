package com.hartwig.hmftools.common.linx;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public final class HomozygousDisruptionFile {
    private static final String DELIMITER = "\t";

    private HomozygousDisruptionFile(){

    }

    @NotNull
    public static List<HomozygousDisruption> read(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    private static List<HomozygousDisruption> fromLines(@NotNull List<String> lines) {
        return lines.stream()
                .skip(1)
                .map(HomozygousDisruptionFile::fromString)
                .collect(toList());
    }

    @NotNull
    private static HomozygousDisruption fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);

        return ImmutableHomozygousDisruption.builder()
                .clusterId(Integer.parseInt(values[0]))
                .gene(values[1])
                .eventType(values[2])
                .build();
    }

}
