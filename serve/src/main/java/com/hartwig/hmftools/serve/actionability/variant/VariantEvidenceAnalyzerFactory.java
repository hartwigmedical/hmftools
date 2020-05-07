package com.hartwig.hmftools.serve.actionability.variant;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class VariantEvidenceAnalyzerFactory {

    private static final String DELIMITER = "\t";

    private VariantEvidenceAnalyzerFactory(){

    }

    @NotNull
    public static VariantEvidenceAnalyzer loadFromFileVariant(@NotNull String actionableVariantTsv) throws IOException {
        final List<ActionableVariant> variants = Lists.newArrayList();
        final List<String> lineVariant = Files.readAllLines(new File(actionableVariantTsv).toPath());

        // Skip header line for variant
        for (String lineVariants : lineVariant.subList(1, lineVariant.size())) {
            variants.add(fromLineVariant(lineVariants));
        }
        return new VariantEvidenceAnalyzer(variants);
    }

    @NotNull
    private static ActionableVariant fromLineVariant(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableVariant.builder()
                .gene(values[0])
                .chromosome(values[1])
                .position(values[2])
                .ref(values[3])
                .alt(values[4])
                .source(values[5])
                .drug(values[6])
                .drugType(values[7])
                .cancerType(values[8])
                .level(values[9])
                .direction(values[10])
                .link(values[11])
                .build();
    }
}
