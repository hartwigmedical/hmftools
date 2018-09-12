package com.hartwig.hmftools.actionability.fusions;

import org.apache.logging.log4j.LogManager;
import org.jetbrains.annotations.NotNull;

public class ActionabilityFusionAnalyzer {
    private static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(ActionabilityFusionAnalyzer.class);
    static final String DELIMITER = "\t";

    @NotNull
    static ActionabilityPromiscuosThree fromLinePromiscuousThree(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionabilityPromiscuosThree.builder()
                .gene(values[0])
                .source(values[1])
                .reference(values[2])
                .drugsName(values[3])
                .drugsType(values[4])
                .cancerType(values[5])
                .level(values[6])
                .hmfLevel(values[7])
                .evidenceType(values[8])
                .significanceSource(values[9])
                .hmfResponse(values[10])
                .build();
    }

    @NotNull
    static ActionabilityPromiscuousFive fromLinePromiscuousFive(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionabilityPromiscuousFive.builder()
                .gene(values[0])
                .source(values[1])
                .reference(values[2])
                .drugsName(values[3])
                .drugsType(values[4])
                .cancerType(values[5])
                .levelSource(values[6])
                .hmfLevel(values[7])
                .evidenceType(values[8])
                .significanceSource(values[9])
                .hmfResponse(values[10])
                .build();
    }

    @NotNull
    static ActionabilityFusionPairs fromLineFusionPairs(@NotNull String line){
        final String[] values = line.split(DELIMITER);
        return ImmutableActionabilityFusionPairs.builder()
                .fiveGene(values[0])
                .threeGene(values[1])
                .source(values[2])
                .reference(values[3])
                .drugsName(values[4])
                .drugsName(values[5])
                .cancerType(values[6])
                .levelSource(values[7])
                .hmfLevel(values[8])
                .evidenceType(values[9])
                .significanceSource(values[10])
                .hmfResponse(values[11])
                .build();
    }

}
