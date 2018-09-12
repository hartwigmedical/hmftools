package com.hartwig.hmftools.actionability.variants;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;


import org.apache.logging.log4j.LogManager;
import org.jetbrains.annotations.NotNull;

public class ActionabilityAnalyzer {
    private static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(ActionabilityAnalyzer.class);
    static final String DELIMITER = "\t";

    @NotNull
    public static ActionabilityAnalyzer loadFromFile(String file) throws IOException {
        final String line =  Files.readAllLines(new File(file).toPath()).get(1);
        fromLine(line);
        LOGGER.info(fromLine(line));
        return null;
    }

    @NotNull
    static ActionabilityVariantsSOC fromLine(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionabilityVariantsSOC.builder()
                .gene(values[0])
                .chromosome(values[1])
                .position(values[2])
                .ref(values[3])
                .alt(values[4])
                .source(values[5])
                .reference(values[6])
                .drug(values[7])
                .drugsType(values[8])
                .cancerType(values[9])
                .levelSource(values[10])
                .levelHmf(values[11])
                .evidenceType(values[12])
                .significanceSource(values[13])
                .hmfResponse(values[14])
                .build();
    }
}
