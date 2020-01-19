package com.hartwig.hmftools.protect.actionability.fusion;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.protect.actionability.util.MultiDrugCurator;

import org.jetbrains.annotations.NotNull;

public final class FusionEvidenceAnalyzerFactory {

    private static final String DELIMITER = "\t";

    private FusionEvidenceAnalyzerFactory() {
    }

    @NotNull
    public static FusionEvidenceAnalyzer loadFromFileFusions(@NotNull String actionableFusionPairsTsv,
            @NotNull String actionablePromiscuousFiveTsv, @NotNull String actionablePromiscuousThreeTsv) throws IOException {
        final List<ActionableFusion> fusionPairs = Lists.newArrayList();
        final List<ActionablePromiscuous> promiscuousFive = Lists.newArrayList();
        final List<ActionablePromiscuous> promiscuousThree = Lists.newArrayList();

        final List<String> lineFusionPairs = Files.readAllLines(new File(actionableFusionPairsTsv).toPath());
        final List<String> linePromiscuousFives = Files.readAllLines(new File(actionablePromiscuousFiveTsv).toPath());
        final List<String> linePromiscuousThrees = Files.readAllLines(new File(actionablePromiscuousThreeTsv).toPath());

        for (String lineFusionPair : lineFusionPairs.subList(1, lineFusionPairs.size())) {
            fusionPairs.add(fromLineFusionPairs(lineFusionPair));
        }

        for (String linePromiscuousFive : linePromiscuousFives.subList(1, linePromiscuousFives.size())) {
            promiscuousFive.add(fromLinePromiscuous(linePromiscuousFive));
        }

        for (String linePromiscuousThree : linePromiscuousThrees.subList(1, linePromiscuousThrees.size())) {
            promiscuousThree.add(fromLinePromiscuous(linePromiscuousThree));
        }

        return new FusionEvidenceAnalyzer(fusionPairs, promiscuousFive, promiscuousThree);
    }

    @NotNull
    private static ActionableFusion fromLineFusionPairs(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableFusion.builder()
                .fiveGene(values[0])
                .threeGene(values[1])
                .source(values[2])
                .reference(values[3])
                .drug(MultiDrugCurator.reformat(values[4]))
                .drugsType(values[5])
                .cancerType(values[6])
                .level(values[8])
                .response(values[11])
                .build();
    }

    @NotNull
    private static ActionablePromiscuous fromLinePromiscuous(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionablePromiscuous.builder()
                .gene(values[0])
                .source(values[1])
                .reference(values[2])
                .drug(MultiDrugCurator.reformat(values[3]))
                .drugsType(values[4])
                .cancerType(values[5])
                .level(values[7])
                .response(values[10])
                .build();
    }
}

