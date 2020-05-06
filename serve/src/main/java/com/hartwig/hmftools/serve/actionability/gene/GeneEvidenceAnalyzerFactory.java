package com.hartwig.hmftools.serve.actionability.gene;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class GeneEvidenceAnalyzerFactory {

    private static final String DELIMITER = "\t";

    private GeneEvidenceAnalyzerFactory(){
    }

    @NotNull
    public static GeneEvidenceAnalyzer loadFromFileGene(@NotNull String actionableGeneTsv) throws IOException {
        final List<ActionableGene> genes = Lists.newArrayList();
        final List<String> lineGene = Files.readAllLines(new File(actionableGeneTsv).toPath());

        // Skip header line for gene
        for (String lineGenes : lineGene.subList(1, lineGene.size())) {
            genes.add(fromLineGene(lineGenes));
        }
        return new GeneEvidenceAnalyzer(genes);
    }

    @NotNull
    private static ActionableGene fromLineGene(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableGene.builder()
                .gene(values[0])
                .type(values[1])
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
