package com.hartwig.hmftools.serve.actionability.variant;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class ActionableVariantFile {

    private static final String DELIMITER = "\t";
    private static final String ACTIONABLE_VARIANT_TSV = "actionableVariants.tsv";

    private ActionableVariantFile() {
    }

    @NotNull
    public static String actionableVariantTsvFilePath(@NotNull String serveActionabilityDir) {
        return serveActionabilityDir + File.separator + ACTIONABLE_VARIANT_TSV;
    }

    public static void writeToActionableVariantsTsv(@NotNull List<ActionableVariant> actionableVariants) {
        // TODO Implement
    }

    @NotNull
    public static List<ActionableVariant> loadFromActionableVariantTsv(@NotNull String actionableVariantTsv) throws IOException {
        List<ActionableVariant> actionableVariants = Lists.newArrayList();
        List<String> lines = Files.readAllLines(new File(actionableVariantTsv).toPath());

        // Skip header line for variant
        for (String line : lines.subList(1, lines.size())) {
            actionableVariants.add(fromLine(line));
        }
        return actionableVariants;
    }

    @NotNull
    private static ActionableVariant fromLine(@NotNull String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableActionableVariant.builder()
                .gene(values[0])
                .chromosome(values[1])
                .position(values[2])
                .ref(values[3])
                .alt(values[4])
                .source(values[5])
                .drug(values[6])
                .cancerType(values[7])
                .doid(values[8])
                .level(values[9])
                .direction(values[10])
                .build();
    }
}
