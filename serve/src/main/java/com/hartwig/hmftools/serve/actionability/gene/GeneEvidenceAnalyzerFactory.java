package com.hartwig.hmftools.serve.actionability.gene;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class GeneEvidenceAnalyzerFactory {

    private static final String DELIMITER = "\t";

    private GeneEvidenceAnalyzerFactory() {
    }

    @NotNull
    public static GeneEvidenceAnalyzer loadFromActionableGeneTsv(@NotNull String actionableGeneTsv) throws IOException {
        List<ActionableGene> actionableGenes = Lists.newArrayList();
        List<String> lines = Files.readAllLines(new File(actionableGeneTsv).toPath());

        // Skip header line for gene
        for (String line : lines.subList(1, lines.size())) {
            actionableGenes.add(fromLine(line));
        }
        return new GeneEvidenceAnalyzer(actionableGenes);
    }

    @NotNull
    private static ActionableGene fromLine(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableGene.builder()
                .gene(values[0])
                .type(values[1])
                .source(values[2])
                .drug(values[3])
                .cancerType(values[4])
                .doid(values[5])
                .level(values[6])
                .direction(values[7])
                .build();
    }
}
