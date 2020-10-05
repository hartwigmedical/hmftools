package com.hartwig.hmftools.serve.actionability.variant;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class ActionableVariantFile {

    private static final String DELIMITER = "\t";
    private static final String ACTIONABLE_VARIANT_TSV = "actionableVariants.tsv";

    private ActionableVariantFile() {
    }

    @NotNull
    public static String actionableVariantTsvPath(@NotNull String serveActionabilityDir) {
        return serveActionabilityDir + File.separator + ACTIONABLE_VARIANT_TSV;
    }

    public static void write(@NotNull String actionableVariantsTsv, @NotNull List<ActionableVariant> actionableVariants)
            throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(actionableVariants));
        Files.write(new File(actionableVariantsTsv).toPath(), lines);
    }

    @NotNull
    public static List<ActionableVariant> loadFromActionableVariantTsv(@NotNull String actionableVariantTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(actionableVariantTsv).toPath());
        // Skip header
        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("gene")
                .add("chromosome")
                .add("position")
                .add("ref")
                .add("alt")
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
    static List<ActionableVariant> fromLines(@NotNull List<String> lines) {
        List<ActionableVariant> actionableVariants = Lists.newArrayList();
        for (String line : lines) {
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
                .treatment(values[6])
                .cancerType(values[7])
                .doid(values[8])
                .level(values[9])
                .direction(values[10])
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull List<ActionableVariant> actionableVariants) {
        List<String> lines = Lists.newArrayList();
        for (ActionableVariant actionableVariant : actionableVariants) {
            lines.add(toLine(actionableVariant));
        }
        return lines;
    }

    @NotNull
    private static String toLine(@NotNull ActionableVariant variant) {
        return new StringJoiner(DELIMITER).add(variant.gene())
                .add(variant.chromosome())
                .add(variant.position())
                .add(variant.ref())
                .add(variant.alt())
                .add(variant.source())
                .add(variant.treatment())
                .add(variant.cancerType())
                .add(variant.doid())
                .add(variant.level())
                .add(variant.direction())
                .toString();
    }
}
