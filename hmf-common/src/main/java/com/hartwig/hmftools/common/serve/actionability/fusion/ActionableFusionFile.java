package com.hartwig.hmftools.common.serve.actionability.fusion;

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

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ActionableFusionFile {

    private static final String ACTIONABLE_FUSION_TSV = "ActionableFusions.tsv";

    private ActionableFusionFile() {
    }

    @NotNull
    public static String actionableFusionTsvPath(@NotNull String serveActionabilityDir, @NotNull RefGenomeVersion refGenomeVersion) {
        return refGenomeVersion.addVersionToFilePath(serveActionabilityDir + File.separator + ACTIONABLE_FUSION_TSV);
    }

    public static void write(@NotNull String actionableFusionTsv, @NotNull Iterable<ActionableFusion> actionableFusions)
            throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(actionableFusions));
        Files.write(new File(actionableFusionTsv).toPath(), lines);
    }

    @NotNull
    public static List<ActionableFusion> read(@NotNull String actionableFusionTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(actionableFusionTsv).toPath());
        // Skip header
        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    private static String header() {
        return new StringJoiner(FIELD_DELIMITER).add("geneUp")
                .add("minExonUp")
                .add("maxExonUp")
                .add("geneDown")
                .add("minExonDown")
                .add("maxExonDown")
                .add(ActionableFileFunctions.header())
                .toString();
    }

    @NotNull
    @VisibleForTesting
    static List<ActionableFusion> fromLines(@NotNull List<String> lines) {
        List<ActionableFusion> actionableFusions = Lists.newArrayList();
        for (String line : lines) {
            actionableFusions.add(fromLine(line));
        }
        return actionableFusions;
    }

    @NotNull
    private static ActionableFusion fromLine(@NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableActionableFusion.builder()
                .from(ActionableFileFunctions.fromLine(values, 6))
                .geneUp(values[0])
                .minExonUp(optionalInteger(values[1]))
                .maxExonUp(optionalInteger(values[2]))
                .geneDown(values[3])
                .minExonDown(optionalInteger(values[4]))
                .maxExonDown(optionalInteger(values[5]))
                .build();
    }

    @Nullable
    private static Integer optionalInteger(@NotNull String value) {
        if (value.isEmpty()) {
            return null;
        }

        return Integer.parseInt(value);
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull Iterable<ActionableFusion> actionableFusions) {
        List<String> lines = Lists.newArrayList();
        for (ActionableFusion actionableFusion : sort(actionableFusions)) {
            lines.add(toLine(actionableFusion));
        }
        return lines;
    }

    @NotNull
    private static List<ActionableFusion> sort(@NotNull Iterable<ActionableFusion> actionableFusions) {
        // Need to make a copy since the input may be immutable and cannot be sorted!
        List<ActionableFusion> sorted = Lists.newArrayList(actionableFusions);
        sorted.sort(new ActionableFusionComparator());

        return sorted;
    }

    @NotNull
    private static String toLine(@NotNull ActionableFusion fusion) {
        return new StringJoiner(FIELD_DELIMITER).add(fusion.geneUp())
                .add(fromOptionalInteger(fusion.minExonUp()))
                .add(fromOptionalInteger(fusion.maxExonUp()))
                .add(fusion.geneDown())
                .add(fromOptionalInteger(fusion.minExonDown()))
                .add(fromOptionalInteger(fusion.maxExonDown()))
                .add(ActionableFileFunctions.toLine(fusion))
                .toString();
    }

    @NotNull
    private static String fromOptionalInteger(@Nullable Integer value) {
        if (value == null) {
            return Strings.EMPTY;
        }
        return Integer.toString(value);
    }
}
