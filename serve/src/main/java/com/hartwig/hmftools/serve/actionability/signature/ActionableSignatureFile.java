package com.hartwig.hmftools.serve.actionability.signature;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.RefGenomeVersion;
import com.hartwig.hmftools.serve.extraction.signature.SignatureName;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ActionableSignatureFile {

    private static final String ACTIONABLE_SIGNATURE_TSV = "ActionableSignatures.tsv";
    private static final String DELIMITER = "\t";

    private ActionableSignatureFile() {
    }

    @NotNull
    public static String actionableSignatureTsvPath(@NotNull String serveActionabilityDir, @NotNull RefGenomeVersion refGenomeVersion) {
        return refGenomeVersion.makeVersioned(serveActionabilityDir + File.separator + ACTIONABLE_SIGNATURE_TSV);
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
        return new StringJoiner(DELIMITER).add("name")
                .add("source")
                .add("treatment")
                .add("cancerType")
                .add("doid")
                .add("level")
                .add("direction")
                .add("url")
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
        String url = values.length > 7 ? values[7] : Strings.EMPTY;

        return ImmutableActionableSignature.builder()
                .name(SignatureName.valueOf(values[0]))
                .source(Knowledgebase.valueOf(values[1]))
                .treatment(values[2])
                .cancerType(values[3])
                .doid(values[4])
                .level(EvidenceLevel.valueOf(values[5]))
                .direction(EvidenceDirection.valueOf(values[6]))
                .url(url)
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
        return new StringJoiner(DELIMITER).add(signature.name().toString())
                .add(signature.source().toString())
                .add(signature.treatment())
                .add(signature.cancerType())
                .add(signature.doid())
                .add(signature.level().toString())
                .add(signature.direction().toString())
                .add(signature.url())
                .toString();
    }
}
