package com.hartwig.hmftools.common.serve.actionability.gene;

import static com.hartwig.hmftools.common.serve.actionability.util.ActionableFileFunctions.FIELD_DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.actionability.util.ActionableFileFunctions;
import com.hartwig.hmftools.common.serve.datamodel.gene.GeneLevelEvent;

import org.jetbrains.annotations.NotNull;

public final class ActionableGeneFile {

    private static final String ACTIONABLE_GENE_TSV = "ActionableGenes.tsv";

    private ActionableGeneFile() {
    }

    @NotNull
    public static String actionableGeneTsvPath(@NotNull String serveActionabilityDir, @NotNull RefGenomeVersion refGenomeVersion) {
        return refGenomeVersion.addVersionToFilePath(serveActionabilityDir + File.separator + ACTIONABLE_GENE_TSV);
    }

    public static void write(@NotNull String actionableGeneTsv, @NotNull Iterable<ActionableGene> actionableGenes) throws IOException {
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
        return new StringJoiner(FIELD_DELIMITER).add("gene").add("event").add(ActionableFileFunctions.header()).toString();
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
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableActionableGene.builder()
                .from(ActionableFileFunctions.fromLine(values, 2))
                .gene(values[0])
                .event(GeneLevelEvent.valueOf(values[1]))
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull Iterable<ActionableGene> actionableGenes) {
        List<String> lines = Lists.newArrayList();
        for (ActionableGene actionableGene : sort(actionableGenes)) {
            lines.add(toLine(actionableGene));
        }
        return lines;
    }

    @NotNull
    private static List<ActionableGene> sort(@NotNull Iterable<ActionableGene> actionableGenes) {
        // Need to make a copy since the input may be immutable and cannot be sorted!
        List<ActionableGene> sorted = Lists.newArrayList(actionableGenes);
        sorted.sort(new ActionableGeneComparator());

        return sorted;
    }

    @NotNull
    private static String toLine(@NotNull ActionableGene gene) {
        return new StringJoiner(FIELD_DELIMITER).add(gene.gene())
                .add(gene.event().toString())
                .add(ActionableFileFunctions.toLine(gene))
                .toString();
    }
}
