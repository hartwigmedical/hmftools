package com.hartwig.hmftools.serve.actionability.characteristic;

import static com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions.FIELD_DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions;
import com.hartwig.hmftools.serve.extraction.characteristic.SignatureName;
import com.hartwig.hmftools.serve.util.RefGenomeVersion;

import org.jetbrains.annotations.NotNull;

public final class ActionableSignatureFile {

    private static final String ACTIONABLE_SIGNATURE_TSV = "ActionableSignatures.tsv";

    private ActionableSignatureFile() {
    }

    @NotNull
    public static String actionableSignatureTsvPath(@NotNull String serveActionabilityDir, @NotNull RefGenomeVersion refGenomeVersion) {
        return refGenomeVersion.addVersionToFilePath(serveActionabilityDir + File.separator + ACTIONABLE_SIGNATURE_TSV);
    }

    public static void write(@NotNull String actionableSignatureTsv, @NotNull Iterable<ActionableSignature> actionableSignatures)
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
        return new StringJoiner(FIELD_DELIMITER).add("name").add(ActionableFileFunctions.header()).toString();
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
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableActionableSignature.builder()
                .from(ActionableFileFunctions.fromLine(values, 1))
                .name(SignatureName.valueOf(values[0]))
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull Iterable<ActionableSignature> actionableSignatures) {
        List<String> lines = Lists.newArrayList();
        for (ActionableSignature actionableSignature : sort(actionableSignatures)) {
            lines.add(toLine(actionableSignature));
        }
        return lines;
    }

    @NotNull
    private static List<ActionableSignature> sort(@NotNull Iterable<ActionableSignature> actionableSignatures) {
        // Need to make a copy since the input may be immutable and cannot be sorted!
        List<ActionableSignature> sorted = Lists.newArrayList(actionableSignatures);
        sorted.sort(new ActionableSignatureComparator());

        return sorted;
    }

    @NotNull
    private static String toLine(@NotNull ActionableSignature signature) {
        return new StringJoiner(FIELD_DELIMITER).add(signature.name().toString()).add(ActionableFileFunctions.toLine(signature)).toString();
    }
}
