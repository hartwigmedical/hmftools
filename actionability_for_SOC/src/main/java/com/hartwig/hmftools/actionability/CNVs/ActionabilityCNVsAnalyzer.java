package com.hartwig.hmftools.actionability.CNVs;

import org.apache.logging.log4j.LogManager;
import org.jetbrains.annotations.NotNull;

public class ActionabilityCNVsAnalyzer {
    private static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(ActionabilityCNVsAnalyzer.class);
    static final String DELIMITER = "\t";

    @NotNull
    static ActionabilityCNVs fromLineCNVs(@NotNull String line){
        final String[] values = line.split(DELIMITER);
        return ImmutableActionabilityCNVs.builder()
                .gene(values[0])
                .cnvType(values[1])
                .source(values[2])
                .reference(values[3])
                .drugsName(values[4])
                .drugsType(values[5])
                .cancerType(values[6])
                .levelSource(values[7])
                .hmfLevel(values[8])
                .evidenceType(values[9])
                .significanceSource(values[10])
                .hmfResponse(values[11])
                .build();

    }
}
