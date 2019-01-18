package com.hartwig.hmftools.svanalysis.visualisation;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;

import org.jetbrains.annotations.NotNull;

public class SvRegionFile {

    private static final String COMMENT = "#";
    private static final String DELIMITER = "\t";

    @NotNull
    public static List<GenomeRegion> readLinks(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    private static List<GenomeRegion> fromLines(@NotNull List<String> lines) {
        final List<GenomeRegion> results = Lists.newArrayList();
        for (final String line : lines) {
            if (!line.startsWith(COMMENT)) {
                results.add(fromString(line));
            }
        }

        return results;

    }

    @NotNull
    private static GenomeRegion fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return GenomeRegionFactory.create(values[0], Long.valueOf(values[1]), Long.valueOf(values[2]));
    }

}
