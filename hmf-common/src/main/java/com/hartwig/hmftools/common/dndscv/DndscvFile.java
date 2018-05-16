package com.hartwig.hmftools.common.dndscv;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

final class DndscvFile {

    private static final String DELIMITER = "\t";

    @NotNull
    public static List<Dndscv> read(@NotNull final String fileName) throws IOException {
        return Files.readAllLines(new File(fileName).toPath()).stream().skip(1).map(DndscvFile::fromString).collect(Collectors.toList());
    }

    @NotNull
    private static Dndscv fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableDndscv.builder()
                .gene(values[0])
                .nSynonymous(Integer.valueOf(values[1]))
                .nMissense(Integer.valueOf(values[2]))
                .nNonSynonymous(Integer.valueOf(values[3]))
                .nSplice(Integer.valueOf(values[4]))
                .nIndel(Integer.valueOf(values[5]))
                .wMissense(Double.valueOf(values[6]))
                .wNonSynonymous(Double.valueOf(values[7]))
                .wSplice(Double.valueOf(values[8]))
                .wIndel(Double.valueOf(values[9]))
                .pScore(Double.valueOf(values[values.length - 2]))
                .qScore(Double.valueOf(values[values.length - 1]))
                .build();
    }
}
