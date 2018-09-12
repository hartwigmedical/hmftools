package com.hartwig.hmftools.actionability.variants;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.jetbrains.annotations.NotNull;

public class ActionabilityAnalyzer {
    private static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(ActionabilityAnalyzer.class);
    static final String DELIMITER = "\t";

    @NotNull
    public static ActionabilityAnalyzer loadFromFile(String file) throws IOException {
        final List<String> line =  Files.readAllLines(new File(file).toPath());
        for (int i = 1; i< line.size(); i++) {
            fromLine(line.get(i));
            LOGGER.info(fromLine(line.get(i)));
        }
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
