package com.hartwig.hmftools.serve.actionability.gene;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class ActionableGeneFile {

    private static final String DELIMITER = "\t";
    private static final String ACTIONABLE_GENE_TSV = "actionableGenes.tsv";

    private ActionableGeneFile() {
    }

    @NotNull
    public static String actionableGeneTsvFilePath(@NotNull String serveActionabilityDir) {
        return serveActionabilityDir + File.separator + ACTIONABLE_GENE_TSV;
    }

    public static void writeToActionableFusionTsv(@NotNull List<ActionableGene> actionableGenes) {
        // TODO Implement
    }

    @NotNull
    public static List<ActionableGene> loadFromActionableGeneTsv(@NotNull String actionableGeneTsv) throws IOException {
        List<ActionableGene> actionableGenes = Lists.newArrayList();
        List<String> lines = Files.readAllLines(new File(actionableGeneTsv).toPath());

        // Skip header line for gene
        for (String line : lines.subList(1, lines.size())) {
            actionableGenes.add(fromLine(line));
        }
        return actionableGenes;
    }

    @NotNull
    private static ActionableGene fromLine(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableGene.builder()
                .gene(values[0])
                .type(values[1])
                .source(values[2])
                .treatment(values[3])
                .cancerType(values[4])
                .doid(values[5])
                .level(values[6])
                .direction(values[7])
                .build();
    }
}
