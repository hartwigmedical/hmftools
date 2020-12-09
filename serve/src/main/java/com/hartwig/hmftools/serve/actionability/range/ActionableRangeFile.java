package com.hartwig.hmftools.serve.actionability.range;

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
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ActionableRangeFile {

    private static final String ACTIONABLE_RANGE_TSV = "ActionableRanges.tsv";
    private static final String DELIMITER = "\t";

    private ActionableRangeFile() {
    }

    @NotNull
    public static String actionableRangeTsvPath(@NotNull String serveActionabilityDir, @NotNull RefGenomeVersion refGenomeVersion) {
        return refGenomeVersion.addVersionToFilePath(serveActionabilityDir + File.separator + ACTIONABLE_RANGE_TSV);
    }

    public static void write(@NotNull String actionableRangeTsv, @NotNull Iterable<ActionableRange> actionableRanges) throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(actionableRanges));
        Files.write(new File(actionableRangeTsv).toPath(), lines);
    }

    @NotNull
    public static List<ActionableRange> read(@NotNull String actionableRangeTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(actionableRangeTsv).toPath());
        // Skip header
        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER).add("gene")
                .add("chromosome")
                .add("start")
                .add("end")
                .add("mutationType")
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
    static List<ActionableRange> fromLines(@NotNull List<String> lines) {
        List<ActionableRange> actionableRanges = Lists.newArrayList();
        for (String line : lines) {
            actionableRanges.add(fromLine(line));
        }
        return actionableRanges;
    }

    @NotNull
    private static ActionableRange fromLine(@NotNull String line) {
        String[] values = line.split(DELIMITER);
        String url = values.length > 11 ? values[11] : Strings.EMPTY;

        return ImmutableActionableRange.builder()
                .gene(values[0])
                .chromosome(values[1])
                .start(Long.parseLong(values[2]))
                .end(Long.parseLong(values[3]))
                .mutationType(MutationTypeFilter.valueOf(values[4]))
                .source(Knowledgebase.valueOf(values[5]))
                .treatment(values[6])
                .cancerType(values[7])
                .doid(values[8])
                .level(EvidenceLevel.valueOf(values[9]))
                .direction(EvidenceDirection.valueOf(values[10]))
                .url(url)
                .build();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull Iterable<ActionableRange> actionableRanges) {
        List<String> lines = Lists.newArrayList();
        for (ActionableRange actionableRange : sort(actionableRanges)) {
            lines.add(toLine(actionableRange));
        }
        return lines;
    }

    @NotNull
    private static List<ActionableRange> sort(@NotNull Iterable<ActionableRange> actionableRanges) {
        // Need to make a copy since the input may be immutable and cannot be sorted!
        List<ActionableRange> sorted = Lists.newArrayList(actionableRanges);
        sorted.sort(new ActionableRangeComparator());

        return sorted;
    }

    @NotNull
    private static String toLine(@NotNull ActionableRange range) {
        return new StringJoiner(DELIMITER).add(range.gene())
                .add(range.chromosome())
                .add(Long.toString(range.start()))
                .add(Long.toString(range.end()))
                .add(range.mutationType().toString())
                .add(range.source().toString())
                .add(range.treatment())
                .add(range.cancerType())
                .add(range.doid())
                .add(range.level().toString())
                .add(range.direction().toString())
                .add(range.url())
                .toString();
    }
}
