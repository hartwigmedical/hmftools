package com.hartwig.hmftools.serve.actionability.hotspot;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.actionability.ActionableEventFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ActionableHotspotFile {

    private static final String DELIMITER = "\t";
    private static final String ACTIONABLE_HOTSPOT_TSV = "actionableHotspots.tsv";

    private ActionableHotspotFile() {
    }

    @NotNull
    public static String actionableHotspotTsvPath(@NotNull String serveActionabilityDir) {
        return serveActionabilityDir + File.separator + ACTIONABLE_HOTSPOT_TSV;
    }

    public static void write(@NotNull String actionableHotspotTsv, @NotNull List<ActionableHotspot> actionableHotspots) throws IOException {
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
        return new StringJoiner(DELIMITER).add("chromosome")
                .add("position")
                .add("ref")
                .add("alt")
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
    static List<ActionableHotspot> fromLines(@NotNull List<String> lines) {
        List<ActionableHotspot> actionableHotspots = Lists.newArrayList();
        for (String line : lines) {
            actionableHotspots.add(fromLine(line));
        }
        return actionableHotspots;
    }

    @NotNull
    private static ActionableHotspot fromLine(@NotNull String line) {
        String[] values = line.split(DELIMITER);
        String url = values.length > 10 ? values[10] : Strings.EMPTY;

        return ImmutableActionableHotspot.builder()
                .chromosome(values[0])
                .position(Long.parseLong(values[1]))
                .ref(values[2])
                .alt(values[3])
                .source(ActionableEventFactory.sourceFromFileValue(values[4]))
                .treatment(values[5])
                .cancerType(values[6])
                .doid(values[7])
                .level(values[8])
                .direction(ActionableEventFactory.directionFromFileValue(values[9]))
                .url(url)
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull List<ActionableHotspot> actionableHotspots) {
        List<String> lines = Lists.newArrayList();
        for (ActionableHotspot actionableHotspot : actionableHotspots) {
            lines.add(toLine(actionableHotspot));
        }
        return lines;
    }

    @NotNull
    private static String toLine(@NotNull ActionableHotspot variant) {
        return new StringJoiner(DELIMITER).add(variant.chromosome())
                .add(Long.toString(variant.position()))
                .add(variant.ref())
                .add(variant.alt())
                .add(variant.source().display())
                .add(variant.treatment())
                .add(variant.cancerType())
                .add(variant.doid())
                .add(variant.level())
                .add(variant.direction().display())
                .add(variant.url())
                .toString();
    }
}
