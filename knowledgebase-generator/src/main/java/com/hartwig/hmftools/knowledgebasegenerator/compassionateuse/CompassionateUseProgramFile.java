package com.hartwig.hmftools.knowledgebasegenerator.compassionateuse;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class CompassionateUseProgramFile {

    private static final String DELIMITER = "\t";

    @NotNull
    public static List<CompassionateUseProgram> read(@NotNull String filePath) throws IOException {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));


    }

    @NotNull
    private static List<CompassionateUseProgram> fromLines (@NotNull List<String> lines) {
        List<CompassionateUseProgram> compassionateUsePrograms = Lists.newArrayList();
        String[] headers = lines.get(0).split(DELIMITER);
        // Skip header line
        for (String line : lines.subList(1, lines.size())) {
            compassionateUsePrograms.add(fromLine(line));
        }
        return compassionateUsePrograms;
    }


    @NotNull
    @VisibleForTesting
    private static CompassionateUseProgram fromLine(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        ImmutableCompassionateUseProgram.Builder builder = ImmutableCompassionateUseProgram.builder()
                .trialAcronymSite(values[0])
                .variants(values[1])
                .source(values[2])
                .drug(values[3])
                .drugType(values[4])
                .cancerType(values[5])
                .level(values[6])
                .direction(values[7])
                .link(values[8]);
        return builder.build();
    }
}
