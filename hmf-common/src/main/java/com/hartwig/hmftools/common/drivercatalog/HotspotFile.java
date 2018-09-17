package com.hartwig.hmftools.common.drivercatalog;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;

import org.jetbrains.annotations.NotNull;

public class HotspotFile {

    private static final String DELIMITER = "\t";

    @NotNull
    public static Multimap<String, GenomePosition> read(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    private static Multimap<String, GenomePosition> fromLines(@NotNull List<String> lines) {
        Multimap<String, GenomePosition> result = ArrayListMultimap.create();
        for (String line : lines) {
            GenomePosition position = fromString(line);
            result.put(position.chromosome(), position);
        }

        return result;
    }

    @NotNull
    private static GenomePosition fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return GenomePositions.create(values[0], Long.valueOf(values[1]));
    }

}
