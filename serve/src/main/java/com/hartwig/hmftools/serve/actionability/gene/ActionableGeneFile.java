package com.hartwig.hmftools.serve.actionability.gene;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.EvidenceLevel;
import com.hartwig.hmftools.serve.actionability.ActionableEventFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ActionableGeneFile {

    private static final String DELIMITER = "\t";
    private static final String ACTIONABLE_GENE_TSV = "actionableGenes.tsv";

    private ActionableGeneFile() {
    }

    @NotNull
    public static String actionableGeneTsvPath(@NotNull String serveActionabilityDir) {
        return serveActionabilityDir + File.separator + ACTIONABLE_GENE_TSV;
    }

    public static void write(@NotNull String actionableGeneTsv, @NotNull List<ActionableGene> actionableGenes) throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(actionableGenes));
        Files.write(new File(actionableGeneTsv).toPath(), lines);
    }

    @NotNull
    public static List<ActionableGene> read(@NotNull String actionableGeneTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(actionableGeneTsv).toPath());
        // Skip header
        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER).add("gene")
                .add("event")
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
    static List<ActionableGene> fromLines(@NotNull List<String> lines) {
        List<ActionableGene> actionableGenes = Lists.newArrayList();
        for (String line : lines) {
            actionableGenes.add(fromLine(line));
        }
        return actionableGenes;
    }

    @NotNull
    private static ActionableGene fromLine(@NotNull String line) {
        String[] values = line.split(DELIMITER);
        String url = values.length > 8 ? values[8] : Strings.EMPTY;

        return ImmutableActionableGene.builder()
                .gene(values[0])
                .event(GeneLevelEvent.valueOf(values[1]))
                .source(ActionableEventFactory.sourceFromFileValue(values[2]))
                .treatment(values[3])
                .cancerType(values[4])
                .doid(values[5])
                .level(EvidenceLevel.valueOf(values[6]))
                .direction(ActionableEventFactory.directionFromFileValue(values[7]))
                .url(url)
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull List<ActionableGene> actionableGenes) {
        List<String> lines = Lists.newArrayList();
        for (ActionableGene actionableGene : actionableGenes) {
            lines.add(toLine(actionableGene));
        }
        return lines;
    }

    @NotNull
    private static String toLine(@NotNull ActionableGene gene) {
        return new StringJoiner(DELIMITER).add(gene.gene())
                .add(gene.event().toString())
                .add(gene.source().display())
                .add(gene.treatment())
                .add(gene.cancerType())
                .add(gene.doid())
                .add(gene.level().toString())
                .add(gene.direction().display())
                .add(gene.url())
                .toString();
    }
}
