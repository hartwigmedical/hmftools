package com.hartwig.hmftools.serve.sources.compassionateuse;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class CompassionateUseProgramFile {

    private static final String DELIMITER = "\t";

    private CompassionateUseProgramFile() {
    }

    @NotNull
    public static List<CompassionateUseProgram> read(@NotNull String compassionateUseProgramTsv) throws IOException {
        return fromLines(Files.readAllLines(new File(compassionateUseProgramTsv).toPath()));
    }

    @NotNull
    private static List<CompassionateUseProgram> fromLines(@NotNull List<String> lines) {
        List<CompassionateUseProgram> compassionateUsePrograms = Lists.newArrayList();
        // Skip header line
        for (String line : lines.subList(1, lines.size())) {
            compassionateUsePrograms.add(fromLine(line));
        }
        return compassionateUsePrograms;
    }

    @NotNull
    private static CompassionateUseProgram fromLine(@NotNull String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableCompassionateUseProgram.builder()
                .trialAcronymSite(values[0])
                .variants(values[1])
                .source(values[2])
                .drug(values[3])
                .drugType(values[4])
                .cancerType(values[5])
                .level(values[6])
                .direction(values[7])
                .link(values[8])
                .build();
    }
}
