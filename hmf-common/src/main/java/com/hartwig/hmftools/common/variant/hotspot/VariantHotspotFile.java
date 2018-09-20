package com.hartwig.hmftools.common.variant.hotspot;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

import org.jetbrains.annotations.NotNull;

public class VariantHotspotFile {

    private static final String DELIMITER = "\t";

    @NotNull
    public static Multimap<String, VariantHotspot> read(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    private static Multimap<String, VariantHotspot> fromLines(@NotNull List<String> lines) {
        Multimap<String, VariantHotspot> result = ArrayListMultimap.create();
        for (String line : lines) {
            VariantHotspot position = fromString(line);
            result.put(position.chromosome(), position);
        }

        return result;
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
