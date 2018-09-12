package com.hartwig.hmftools.actionability.variants;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ActionabilityVariantsAnalyzer {
    private static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(ActionabilityVariantsAnalyzer.class);
    static final String DELIMITER = "\t";

    @NotNull
    public static ActionabilityVariantsAnalyzer loadFromFile(String file) throws IOException {
        final List<ActionabilityVariantsSOC> VariantsFile = new ArrayList<>();
        final List<String> line =  Files.readAllLines(new File(file).toPath());
        for (int i = 1; i< line.size(); i++) {
            fromLineVariants(line.get(i));
            VariantsFile.add(fromLineVariants(line.get(i)));
            LOGGER.info(fromLineVariants(line.get(i)));
        }
        return null;
    }

    @NotNull
    static ActionabilityVariantsSOC fromLineVariants(@NotNull String line) {
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

    @NotNull
    static ActionabilityRanges fromLineRanges(@NotNull String line) {
        final String [] values = line.split(DELIMITER);
        return ImmutableActionabilityRanges.builder()
                .mutationTranscript(values[0])
                .chromosome(values[1])
                .start(values[2])
                .stop(values[3])
                .geneTranscript(values[4])
                .source(values[5])
                .reference(values[6])
                .drugsName(values[7])
                .drugsType(values[8])
                .cancerType(values[9])
                .levelSource(values[10])
                .hmfLevel(values[11])
                .evidenceType(values[12])
                .significanceSource(values[13])
                .hmfResponse(values[14])
                .build();
    }
}
