package com.hartwig.hmftools.serve.actionability.signature;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class ActionableSignatureFile {

    private static final String DELIMITER = "\t";
    private static final String ACTIONABLE_SIGNATURE_TSV = "actionableSignatures.tsv";

    private ActionableSignatureFile(){
    }

    @NotNull
    public static String actionableSignatureTsvFilePath(@NotNull String serveActionabilityDir) {
        return serveActionabilityDir + File.separator + ACTIONABLE_SIGNATURE_TSV;
    }

    public static void writeToActionableSignatureTsv(@NotNull List<ActionableSignature> actionableSignatures) {
        // TODO Implement
    }

    @NotNull
    public static List<ActionableSignature> loadFromActionableSignatureTsv(@NotNull String actionableSignatureTsv) throws IOException {
         List<ActionableSignature> actionableSignatures = Lists.newArrayList();
         List<String> lines = Files.readAllLines(new File(actionableSignatureTsv).toPath());

        // Skip header line for signatures
        for (String line : lines.subList(1, lines.size())) {
            actionableSignatures.add(fromLine(line));
        }
        return actionableSignatures;
    }

    @NotNull
    private static ActionableSignature fromLine(@NotNull String line) {
         String[] values = line.split(DELIMITER);
        return ImmutableActionableSignature.builder()
                .signature(values[0])
                .source(values[1])
                .treatment(values[2])
                .cancerType(values[3])
                .doid(values[4])
                .level(values[5])
                .direction(values[6])
                .build();
    }
}
