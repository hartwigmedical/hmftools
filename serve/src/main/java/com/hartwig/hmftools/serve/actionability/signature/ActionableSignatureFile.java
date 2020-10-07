package com.hartwig.hmftools.serve.actionability.signature;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.actionability.ActionableEventFactory;

import org.jetbrains.annotations.NotNull;

public final class ActionableSignatureFile {

    private static final String DELIMITER = "\t";
    private static final String ACTIONABLE_SIGNATURE_TSV = "actionableSignatures.tsv";

    private ActionableSignatureFile() {
    }

    @NotNull
    public static String actionableSignatureTsvPath(@NotNull String serveActionabilityDir) {
        return serveActionabilityDir + File.separator + ACTIONABLE_SIGNATURE_TSV;
    }

    public static void write(@NotNull String actionableSignatureTsv, @NotNull List<ActionableSignature> actionableSignatures)
            throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(actionableSignatures));
        Files.write(new File(actionableSignatureTsv).toPath(), lines);
    }

    @NotNull
    public static List<ActionableSignature> read(@NotNull String actionableSignatureTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(actionableSignatureTsv).toPath());
        // Skip header
        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("signature")
                .add("source")
                .add("treatment")
                .add("cancerType")
                .add("doid")
                .add("level")
                .add("direction")
                .toString();
    }

    @NotNull
    @VisibleForTesting
    static List<ActionableSignature> fromLines(@NotNull List<String> lines) {
        List<ActionableSignature> actionableSignatures = Lists.newArrayList();
        for (String line : lines) {
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
                .direction(ActionableEventFactory.directionFromFileValue(values[6]))
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull List<ActionableSignature> actionableSignatures) {
        List<String> lines = Lists.newArrayList();
        for (ActionableSignature actionableSignature : actionableSignatures) {
            lines.add(toLine(actionableSignature));
        }
        return lines;
    }

    @NotNull
    private static String toLine(@NotNull ActionableSignature signature) {
        return new StringJoiner(DELIMITER).add(signature.signature())
                .add(signature.source())
                .add(signature.treatment())
                .add(signature.cancerType())
                .add(signature.doid())
                .add(signature.level())
                .add(signature.direction().display())
                .toString();
    }
}
