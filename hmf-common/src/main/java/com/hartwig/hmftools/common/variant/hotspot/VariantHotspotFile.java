package com.hartwig.hmftools.common.variant.hotspot;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;

public class VariantHotspotFile {

    private static final String DELIMITER = "\t";

    @NotNull
    public static List<VariantHotspot> read(@NotNull final String fileName) throws IOException {
        return Files.readAllLines(new File(fileName).toPath()).stream().map(VariantHotspotFile::fromString).collect(Collectors.toList());
    }

    @NotNull
    private static VariantHotspot fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableVariantHotspot.builder()
                .chromosome(values[0])
                .position(Long.valueOf(values[1]))
                .ref(values[2])
                .alt(values[3])
                .build();
    }
}
