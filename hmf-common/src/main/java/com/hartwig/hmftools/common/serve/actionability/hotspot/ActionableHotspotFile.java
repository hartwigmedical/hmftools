package com.hartwig.hmftools.common.serve.actionability.hotspot;

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

import org.jetbrains.annotations.NotNull;

public final class ActionableHotspotFile {

    private static final String ACTIONABLE_HOTSPOT_TSV = "ActionableHotspots.tsv";

    private ActionableHotspotFile() {
    }

    @NotNull
    public static String actionableHotspotTsvPath(@NotNull String serveActionabilityDir, @NotNull RefGenomeVersion refGenomeVersion) {
        return refGenomeVersion.addVersionToFilePath(serveActionabilityDir + File.separator + ACTIONABLE_HOTSPOT_TSV);
    }

    public static void write(@NotNull String actionableHotspotTsv, @NotNull Iterable<ActionableHotspot> actionableHotspots)
            throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(actionableHotspots));
        Files.write(new File(actionableHotspotTsv).toPath(), lines);
    }

    @NotNull
    public static List<ActionableHotspot> read(@NotNull String actionableHotspotTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(actionableHotspotTsv).toPath());
        // Skip header
        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    private static String header() {
        return new StringJoiner(FIELD_DELIMITER).add("chromosome")
                .add("position")
                .add("ref")
                .add("alt")
                .add(ActionableFileFunctions.header())
                .toString();
    }

    @NotNull
    @VisibleForTesting
    static List<ActionableHotspot> fromLines(@NotNull List<String> lines) {
        List<ActionableHotspot> actionableHotspots = Lists.newArrayList();
        for (String line : lines) {
            actionableHotspots.add(fromLine(line));
        }
        return actionableHotspots;
    }

    @NotNull
    private static ActionableHotspot fromLine(@NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableActionableHotspot.builder()
                .from(ActionableFileFunctions.fromLine(values, 4))
                .chromosome(values[0])
                .position(Integer.parseInt(values[1]))
                .ref(values[2])
                .alt(values[3])
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull Iterable<ActionableHotspot> actionableHotspots) {
        List<String> lines = Lists.newArrayList();
        for (ActionableHotspot actionableHotspot : sort(actionableHotspots)) {
            lines.add(toLine(actionableHotspot));
        }
        return lines;
    }

    @NotNull
    private static List<ActionableHotspot> sort(@NotNull Iterable<ActionableHotspot> actionableHotspots) {
        // Need to make a copy since the input may be immutable and cannot be sorted!
        List<ActionableHotspot> sorted = Lists.newArrayList(actionableHotspots);
        sorted.sort(new ActionableHotspotComparator());

        return sorted;
    }

    @NotNull
    private static String toLine(@NotNull ActionableHotspot variant) {
        return new StringJoiner(FIELD_DELIMITER).add(variant.chromosome())
                .add(Long.toString(variant.position()))
                .add(variant.ref())
                .add(variant.alt())
                .add(ActionableFileFunctions.toLine(variant))
                .toString();
    }
}
