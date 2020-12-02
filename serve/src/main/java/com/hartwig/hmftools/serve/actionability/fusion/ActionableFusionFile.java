package com.hartwig.hmftools.serve.actionability.fusion;

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

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ActionableFusionFile {

    private static final String ACTIONABLE_FUSION_TSV = "ActionableFusions.tsv";
    private static final String DELIMITER = "\t";

    private ActionableFusionFile() {
    }

    @NotNull
    public static String actionableFusionTsvPath(@NotNull String serveActionabilityDir, @NotNull RefGenomeVersion refGenomeVersion) {
        return refGenomeVersion.makeVersioned(serveActionabilityDir + File.separator + ACTIONABLE_FUSION_TSV);
    }

    public static void write(@NotNull String actionableFusionTsv, @NotNull List<ActionableFusion> actionableFusions) throws IOException {
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
        return new StringJoiner(DELIMITER).add("geneUp")
                .add("minExonUp")
                .add("maxExonUp")
                .add("geneDown")
                .add("minExonDown")
                .add("maxExonDown")
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
    static List<ActionableFusion> fromLines(@NotNull List<String> lines) {
        List<ActionableFusion> actionableFusions = Lists.newArrayList();
        for (String line : lines) {
            actionableFusions.add(fromLine(line));
        }
        return actionableFusions;
    }

    @NotNull
    private static ActionableFusion fromLine(@NotNull String line) {
        String[] values = line.split(DELIMITER);
        String url = values.length > 12 ? values[12] : Strings.EMPTY;

        return ImmutableActionableFusion.builder()
                .geneUp(values[0])
                .minExonUp(optionalInteger(values[1]))
                .maxExonDown(optionalInteger(values[2]))
                .geneDown(values[3])
                .minExonDown(optionalInteger(values[4]))
                .maxExonDown(optionalInteger(values[5]))
                .source(Knowledgebase.valueOf(values[6]))
                .treatment(values[7])
                .cancerType(values[8])
                .doid(values[9])
                .level(EvidenceLevel.valueOf(values[10]))
                .direction(EvidenceDirection.valueOf(values[11]))
                .url(url)
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
    static List<String> toLines(@NotNull List<ActionableFusion> actionableFusions) {
        List<String> lines = Lists.newArrayList();
        for (ActionableFusion actionableFusion : actionableFusions) {
            lines.add(toLine(actionableFusion));
        }
        return lines;
    }

    @NotNull
    private static String toLine(@NotNull ActionableFusion fusion) {
        return new StringJoiner(DELIMITER).add(fusion.geneUp())
                .add(fromOptionalInteger(fusion.minExonUp()))
                .add(fromOptionalInteger(fusion.maxExonUp()))
                .add(fusion.geneDown())
                .add(fromOptionalInteger(fusion.minExonDown()))
                .add(fromOptionalInteger(fusion.maxExonDown()))
                .add(fusion.source().toString())
                .add(fusion.treatment())
                .add(fusion.cancerType())
                .add(fusion.doid())
                .add(fusion.level().toString())
                .add(fusion.direction().toString())
                .add(fusion.url())
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
